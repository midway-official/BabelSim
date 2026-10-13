#include "babelsim/parallel.h"
#include "babelsim/mesh_io.h"
#include "babelsim/mpi_support.h"
#include "internal/mesh_access.h"
#include "internal/mesh_partition.h"

#include <algorithm>
#include <deque>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace babelsim {
namespace {

void broadcastPatch(std::vector<PatchSpec>& patches, const ParallelContext& parallel) {
    Index count = parallel.rank == 0 ? static_cast<Index>(patches.size()) : 0;
    detail::checkMpi(MPI_Bcast(&count, 1, MPI_INT, 0, parallel.communicator),
                     "MPI_Bcast(patch count)");
    if (count <= 0 || count > 65536)
        throw std::runtime_error("mesh patch count is invalid");
    if (parallel.rank != 0)
        patches.resize(static_cast<std::size_t>(count));
    for (Index patch = 0; patch < count; ++patch) {
        PatchSpec& specification = patches[static_cast<std::size_t>(patch)];
        int kind = parallel.rank == 0 ? static_cast<int>(specification.kind) : 0;
        int length = parallel.rank == 0 ? static_cast<int>(specification.name.size()) : 0;
        detail::checkMpi(MPI_Bcast(&kind, 1, MPI_INT, 0, parallel.communicator),
                         "MPI_Bcast(patch kind)");
        detail::checkMpi(MPI_Bcast(&length, 1, MPI_INT, 0, parallel.communicator),
                         "MPI_Bcast(patch name length)");
        if (length <= 0 || length > 4096 || kind < static_cast<int>(PatchKind::Generic) ||
            kind > static_cast<int>(PatchKind::Processor)) {
            throw std::runtime_error("mesh patch metadata is invalid");
        }
        if (parallel.rank != 0)
            specification.name.resize(static_cast<std::size_t>(length));
        detail::checkMpi(
            MPI_Bcast(specification.name.data(), length, MPI_CHAR, 0, parallel.communicator),
            "MPI_Bcast(patch name)");
        specification.kind = static_cast<PatchKind>(kind);
    }
}

// Each rank receives only its owned cells and the stencil-required halo.
// The complete mesh and the global graph partition exist on rank zero only.
struct PartitionInput {
    std::vector<Vec3> vertices;
    std::vector<PolyhedralFaceSpec> faces;
    std::vector<PatchSpec> patches;
    std::vector<Index> cell_ids, cell_owners, cell_depths, face_ids, face_owners;

    Mesh finalize(Index global_cells, Index layers, Index rank) {
        Mesh local = Mesh::polyhedral(std::move(vertices), std::move(faces), std::move(patches));
        detail::MeshAccess::setPartition(local, global_cells, layers, std::move(cell_ids),
                                         std::move(cell_owners), std::move(cell_depths),
                                         std::move(face_ids), std::move(face_owners), rank);
        local.validate();
        return local;
    }
};

void collectiveMeshError(const std::string& local_error, const ParallelContext& parallel) {
    int candidate = local_error.empty() ? parallel.size : parallel.rank;
    int source = parallel.size;
    detail::checkMpi(MPI_Allreduce(&candidate, &source, 1, MPI_INT, MPI_MIN, parallel.communicator),
                     "MPI_Allreduce(mesh failure rank)");
    if (source == parallel.size)
        return;
    std::string message = parallel.rank == source ? local_error.substr(0, 16384) : std::string{};
    int length = static_cast<int>(message.size());
    detail::checkMpi(MPI_Bcast(&length, 1, MPI_INT, source, parallel.communicator),
                     "MPI_Bcast(mesh failure length)");
    message.resize(static_cast<std::size_t>(length));
    detail::checkMpi(MPI_Bcast(message.data(), length, MPI_CHAR, source, parallel.communicator),
                     "MPI_Bcast(mesh failure message)");
    throw std::runtime_error("distributed mesh setup failed on rank " + std::to_string(source) +
                             ": " + message);
}

struct PackedPartition {
    // Wire format: five counts followed by packed coordinates and topology.
    enum CountSlot { Vertices, Faces, FacePoints, Cells, GlobalCells, CountSlots };
    Index counts[CountSlots]{};
    std::vector<double> coordinates;
    std::vector<Index> topology;
    std::vector<PatchSpec> patches;

    void allocate() {
        if (counts[PackedPartition::Vertices] <= 0 || counts[PackedPartition::Faces] <= 0 ||
            counts[PackedPartition::FacePoints] < 3 || counts[PackedPartition::Cells] <= 0 ||
            counts[PackedPartition::GlobalCells] < counts[PackedPartition::Cells])
            throw std::runtime_error("local partition counts are invalid");
        const std::size_t coordinate_count =
            static_cast<std::size_t>(counts[PackedPartition::Vertices]) * 3U;
        const std::size_t topology_count =
            static_cast<std::size_t>(counts[PackedPartition::Faces]) * 6U + 1U +
            static_cast<std::size_t>(counts[PackedPartition::FacePoints]) +
            static_cast<std::size_t>(counts[PackedPartition::Cells]) * 3U;
        detail::mpiCount(coordinate_count, "local coordinates");
        detail::mpiCount(topology_count, "local topology");
        coordinates.resize(coordinate_count);
        topology.resize(topology_count);
    }

    PackedPartition() = default;
    PackedPartition(const PartitionInput& input, Index global_cells) {
        counts[PackedPartition::Vertices] =
            detail::mpiCount(input.vertices.size(), "local vertices");
        counts[PackedPartition::Faces] = detail::mpiCount(input.faces.size(), "local faces");
        std::size_t points = 0;
        for (const auto& face : input.faces)
            points += face.vertices.size();
        counts[PackedPartition::FacePoints] = detail::mpiCount(points, "local face points");
        counts[PackedPartition::Cells] = detail::mpiCount(input.cell_ids.size(), "local cells");
        counts[PackedPartition::GlobalCells] = global_cells;
        allocate();
        patches = input.patches;
        std::size_t i = 0;
        for (const Vec3& point : input.vertices) {
            coordinates[i++] = point.x;
            coordinates[i++] = point.y;
            coordinates[i++] = point.z;
        }
        // Topology layout (pack and unpack must agree): face-point offsets,
        // point IDs, owner/neighbour/patch, face IDs/ranks, cell IDs/ranks/depths.
        i = 0;
        Index offset = 0;
        topology[i++] = offset;
        for (const auto& face : input.faces) {
            offset += static_cast<Index>(face.vertices.size());
            topology[i++] = offset;
        }
        for (const auto& face : input.faces)
            for (Index vertex : face.vertices)
                topology[i++] = vertex;
        for (const auto& face : input.faces)
            topology[i++] = face.owner;
        for (const auto& face : input.faces)
            topology[i++] = face.neighbour;
        for (const auto& face : input.faces)
            topology[i++] = face.patch;
        for (const auto* values : {&input.face_ids, &input.face_owners, &input.cell_ids,
                                   &input.cell_owners, &input.cell_depths})
            for (Index value : *values)
                topology[i++] = value;
        if (i != topology.size())
            throw std::logic_error("packed partition layout disagrees");
    }

    PartitionInput unpack() {
        PartitionInput input;
        input.vertices.resize(static_cast<std::size_t>(counts[PackedPartition::Vertices]));
        for (std::size_t i = 0; i < input.vertices.size(); ++i)
            input.vertices[i] = {coordinates[3U * i], coordinates[3U * i + 1U],
                                 coordinates[3U * i + 2U]};
        input.faces.resize(static_cast<std::size_t>(counts[PackedPartition::Faces]));
        std::size_t i = static_cast<std::size_t>(counts[PackedPartition::Faces]) + 1U;
        for (std::size_t f = 0; f < input.faces.size(); ++f)
            input.faces[f].vertices.assign(topology.begin() + i + topology[f],
                                           topology.begin() + i + topology[f + 1U]);
        i += static_cast<std::size_t>(counts[PackedPartition::FacePoints]);
        for (auto& face : input.faces)
            face.owner = topology[i++];
        for (auto& face : input.faces)
            face.neighbour = topology[i++];
        for (auto& face : input.faces)
            face.patch = topology[i++];
        for (auto* values : {&input.face_ids, &input.face_owners, &input.cell_ids,
                             &input.cell_owners, &input.cell_depths}) {
            const std::size_t count =
                (values == &input.face_ids || values == &input.face_owners)
                    ? static_cast<std::size_t>(counts[PackedPartition::Faces])
                    : static_cast<std::size_t>(counts[PackedPartition::Cells]);
            values->assign(topology.begin() + i, topology.begin() + i + count);
            i += count;
        }
        input.patches = std::move(patches);
        return input;
    }
};

PartitionInput partitionInput(const Mesh& global, const std::vector<Index>& owners, int rank,
                              int size, Index ghost_layers) {
    if (size <= 0 || rank < 0 || rank >= size || ghost_layers < 3) {
        throw std::invalid_argument("mesh partition rank, size, or halo width is invalid");
    }
    if (detail::ownedCellCount(global) != global.cellCount()) {
        throw std::invalid_argument("domain decomposition requires a complete mesh");
    }
    if (global.globalCellCount() < size)
        throw std::invalid_argument("MPI has more ranks than mesh cells");

    std::vector<Index> depth(static_cast<std::size_t>(global.cellCount()), invalid_index);
    std::deque<Index> frontier;
    for (Index cell = 0; cell < global.cellCount(); ++cell) {
        const Index owner = owners[static_cast<std::size_t>(cell)];
        if (owner == rank) {
            depth[static_cast<std::size_t>(cell)] = 0;
            frontier.push_back(cell);
        }
    }
    if (frontier.empty())
        throw std::logic_error("partition owns no global cells");
    while (!frontier.empty()) {
        const Index cell = frontier.front();
        frontier.pop_front();
        const Index current_depth = depth[static_cast<std::size_t>(cell)];
        if (current_depth >= ghost_layers)
            continue;
        for (Index neighbour : global.cellNeighbours(cell)) {
            if (neighbour == invalid_index ||
                depth[static_cast<std::size_t>(neighbour)] != invalid_index)
                continue;
            depth[static_cast<std::size_t>(neighbour)] = current_depth + 1;
            frontier.push_back(neighbour);
        }
    }

    std::vector<Index> selected;
    std::vector<Index> source_to_local(static_cast<std::size_t>(global.cellCount()), invalid_index);
    for (Index cell = 0; cell < global.cellCount(); ++cell) {
        if (depth[static_cast<std::size_t>(cell)] != invalid_index) {
            source_to_local[static_cast<std::size_t>(cell)] = static_cast<Index>(selected.size());
            selected.push_back(cell);
        }
    }
    std::vector<Index> vertex_to_local(static_cast<std::size_t>(global.vertexCount()),
                                       invalid_index);
    std::vector<Vec3> vertices;
    std::vector<Index> cell_ids;
    std::vector<Index> cell_owners;
    std::vector<Index> cell_depths;
    vertices.reserve(selected.size() * 4U);
    for (Index source : selected) {
        cell_ids.push_back(detail::globalCellId(global, source));
        cell_owners.push_back(owners[static_cast<std::size_t>(source)]);
        cell_depths.push_back(depth[static_cast<std::size_t>(source)]);
    }

    std::vector<PatchSpec> patches;
    patches.reserve(static_cast<std::size_t>(global.patchCount()) + static_cast<std::size_t>(size));
    for (Index patch = 0; patch < global.patchCount(); ++patch) {
        patches.push_back({global.patchName(patch), global.patchKind(patch)});
    }
    std::vector<Index> processor_patch(static_cast<std::size_t>(size), invalid_index);
    std::vector<Index> face_patch(static_cast<std::size_t>(global.faceCount()), invalid_index);
    for (Index face = 0; face < global.faceCount(); ++face) {
        const Index owner = global.owner(face);
        const Index neighbour = global.neighbour(face);
        const bool owner_inside = source_to_local[static_cast<std::size_t>(owner)] != invalid_index;
        const bool neighbour_inside =
            neighbour != invalid_index &&
            source_to_local[static_cast<std::size_t>(neighbour)] != invalid_index;
        if (owner_inside == neighbour_inside)
            continue;
        if (neighbour == invalid_index) {
            face_patch[static_cast<std::size_t>(face)] = detail::meshData(global).face_patch[face];
        } else {
            const Index outside = owner_inside ? neighbour : owner;
            const Index remote = owners[static_cast<std::size_t>(outside)];
            Index& patch = processor_patch[static_cast<std::size_t>(remote)];
            if (patch == invalid_index) {
                patch = static_cast<Index>(patches.size());
                patches.push_back({"processor_" + std::to_string(remote), PatchKind::Processor});
            }
            face_patch[static_cast<std::size_t>(face)] = patch;
        }
    }
    std::vector<PolyhedralFaceSpec> local_faces;
    std::vector<Index> local_global_face_ids;
    std::vector<Index> local_face_owner_ranks;
    for (Index face = 0; face < global.faceCount(); ++face) {
        const Index patch = face_patch[static_cast<std::size_t>(face)];
        const Index owner = global.owner(face);
        const Index neighbour = global.neighbour(face);
        const bool attached =
            source_to_local[static_cast<std::size_t>(owner)] != invalid_index ||
            (neighbour != invalid_index &&
             source_to_local[static_cast<std::size_t>(neighbour)] != invalid_index);
        if (!attached)
            continue;
        const bool owner_inside = source_to_local[static_cast<std::size_t>(owner)] != invalid_index;
        const Index local_owner = owner_inside
                                      ? source_to_local[static_cast<std::size_t>(owner)]
                                      : source_to_local[static_cast<std::size_t>(neighbour)];
        const bool neighbour_inside =
            neighbour != invalid_index &&
            source_to_local[static_cast<std::size_t>(neighbour)] != invalid_index;
        PolyhedralFaceSpec specification;
        specification.owner = local_owner;
        const bool both_inside = owner_inside && neighbour_inside;
        specification.neighbour =
            both_inside ? source_to_local[static_cast<std::size_t>(neighbour)] : invalid_index;
        specification.patch = both_inside ? invalid_index : patch;
        if (specification.neighbour == invalid_index && specification.patch == invalid_index) {
            throw std::logic_error("local face " + std::to_string(face) +
                                   " lost boundary/processor patch owner=" + std::to_string(owner) +
                                   " neighbour=" + std::to_string(neighbour) +
                                   " owner_inside=" + std::to_string(owner_inside) +
                                   " neighbour_inside=" + std::to_string(neighbour_inside));
        }
        for (Index original : global.facePoints(face)) {
            Index& mapped = vertex_to_local[static_cast<std::size_t>(original)];
            if (mapped == invalid_index) {
                mapped = static_cast<Index>(vertices.size());
                vertices.push_back(global.vertex(original));
            }
            specification.vertices.push_back(mapped);
        }
        if (!owner_inside && neighbour_inside) {
            std::reverse(specification.vertices.begin() + 1, specification.vertices.end());
        }
        local_faces.push_back(std::move(specification));
        local_global_face_ids.push_back(detail::globalFaceId(global, face));
        // Communication ownership is a global face property.  A local copy
        // may reverse its geometric owner when the global owner is outside
        // this rank, but requests must still go to the rank owning the global
        // owner cell so every copy shares one authoritative face record.
        local_face_owner_ranks.push_back(owners[static_cast<std::size_t>(owner)]);
    }
    return {std::move(vertices),
            std::move(local_faces),
            std::move(patches),
            std::move(cell_ids),
            std::move(cell_owners),
            std::move(cell_depths),
            std::move(local_global_face_ids),
            std::move(local_face_owner_ranks)};
}

} // namespace

Mesh decompose(const Mesh& global, const ParallelContext& parallel, Index ghost_layers) {
    parallel.validate();
    global.validate();
    if (detail::ownedCellCount(global) != global.cellCount())
        throw std::invalid_argument("domain decomposition requires a complete mesh");
    if (ghost_layers < 3)
        throw std::invalid_argument("ghostLayers must be >= 3");
    if (parallel.size == 1)
        return global;
    const auto owners = detail::partitionOwners(global, parallel.size);
    return partitionInput(global, owners, parallel.rank, parallel.size, ghost_layers)
        .finalize(global.globalCellCount(), ghost_layers, parallel.rank);
}

Mesh readDistributedMesh(const std::filesystem::path& path, const ParallelContext& parallel,
                         Index ghost_layers) {
    parallel.validate();
    if (ghost_layers < 3)
        throw std::invalid_argument("ghostLayers must be >= 3");
    if (!parallel.distributed())
        return readMeshFile(path);
    Mesh global;
    std::string error;
    std::vector<Index> owners;
    if (parallel.rank == 0) {
        try {
            detail::MeshAccess::replace(global, readMeshFile(path));
            global.validate();
            owners = detail::partitionOwners(global, parallel.size);
        } catch (const std::exception& exception) {
            error = exception.what();
        }
    }
    collectiveMeshError(error, parallel);
    PackedPartition received;
    for (int peer = 1; peer < parallel.size; ++peer) {
        PackedPartition outgoing;
        error.clear();
        if (parallel.rank == 0) {
            try {
                outgoing = PackedPartition(
                    partitionInput(global, owners, peer, parallel.size, ghost_layers),
                    global.globalCellCount());
            } catch (const std::exception& exception) {
                error = exception.what();
            }
        }
        collectiveMeshError(error, parallel);
        Index counts[PackedPartition::CountSlots]{};
        if (parallel.rank == 0)
            std::copy(std::begin(outgoing.counts), std::end(outgoing.counts), counts);
        detail::checkMpi(
            MPI_Bcast(counts, PackedPartition::CountSlots, MPI_INT, 0, parallel.communicator),
            "MPI_Bcast(local partition counts)");
        std::vector<PatchSpec> patches =
            parallel.rank == 0 ? outgoing.patches : std::vector<PatchSpec>{};
        broadcastPatch(patches, parallel);
        error.clear();
        if (parallel.rank == peer) {
            try {
                std::copy(std::begin(counts), std::end(counts), received.counts);
                received.patches = std::move(patches);
                received.allocate();
            } catch (const std::exception& exception) {
                error = exception.what();
            }
        }
        collectiveMeshError(error, parallel);
        if (parallel.rank == 0) {
            detail::checkMpi(
                MPI_Send(outgoing.coordinates.data(),
                         detail::mpiCount(outgoing.coordinates.size(), "local coordinates"),
                         MPI_DOUBLE, peer, 170, parallel.communicator),
                "MPI_Send(local coordinates)");
            detail::checkMpi(MPI_Send(outgoing.topology.data(),
                                      detail::mpiCount(outgoing.topology.size(), "local topology"),
                                      MPI_INT, peer, 171, parallel.communicator),
                             "MPI_Send(local topology)");
        } else if (parallel.rank == peer) {
            detail::checkMpi(
                MPI_Recv(received.coordinates.data(),
                         detail::mpiCount(received.coordinates.size(), "local coordinates"),
                         MPI_DOUBLE, 0, 170, parallel.communicator, MPI_STATUS_IGNORE),
                "MPI_Recv(local coordinates)");
            detail::checkMpi(MPI_Recv(received.topology.data(),
                                      detail::mpiCount(received.topology.size(), "local topology"),
                                      MPI_INT, 0, 171, parallel.communicator, MPI_STATUS_IGNORE),
                             "MPI_Recv(local topology)");
        }
    }
    Mesh local;
    error.clear();
    try {
        const Index global_cells = parallel.rank == 0
                                       ? global.globalCellCount()
                                       : received.counts[PackedPartition::GlobalCells];
        PartitionInput input = parallel.rank == 0
                                   ? partitionInput(global, owners, 0, parallel.size, ghost_layers)
                                   : received.unpack();
        // Release the full root mesh and receive buffers before constructing
        // local geometry. No global-size allocation survives on worker ranks.
        detail::MeshAccess::replace(global, Mesh{});
        received = PackedPartition{};
        detail::MeshAccess::replace(local,
                                    input.finalize(global_cells, ghost_layers, parallel.rank));
    } catch (const std::exception& exception) {
        error = exception.what();
    }
    collectiveMeshError(error, parallel);
    return local;
}

} // namespace babelsim
