#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/parallel.h"

#include "babelsim/mesh_io.h"
#include "babelsim/mpi_support.h"

#include <algorithm>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace babelsim {
namespace {

using Quad = std::array<Index, 4>;

Quad canonical(Quad vertices) {
    std::sort(vertices.begin(), vertices.end());
    return vertices;
}

void requireCommunicator(MPI_Comm communicator, const char* operation) {
    if (communicator == MPI_COMM_NULL) {
        throw std::invalid_argument(std::string(operation) + " received MPI_COMM_NULL");
    }
    detail::requireMpiReady(operation);
}

void requireCollectiveContext(const ParallelContext& parallel, const char* operation) {
    if (parallel.size <= 0 || parallel.rank < 0 || parallel.rank >= parallel.size) {
        throw std::invalid_argument("parallel rank or size is invalid");
    }
    if (parallel.communicator == MPI_COMM_NULL) {
        if (parallel.size != 1 || parallel.rank != 0) {
            throw std::invalid_argument("distributed collective has MPI_COMM_NULL");
        }
        return;
    }
    detail::requireMpiReady(operation);
}

void broadcastPatch(std::vector<PatchSpec>& patches, const ParallelContext& parallel) {
    Index count = parallel.rank == 0 ? static_cast<Index>(patches.size()) : 0;
    detail::checkMpi(MPI_Bcast(&count, 1, MPI_INT, 0, parallel.communicator),
                     "MPI_Bcast(patch count)");
    if (count <= 0 || count > 65536) throw std::runtime_error("mesh patch count is invalid");
    if (parallel.rank != 0) patches.resize(static_cast<std::size_t>(count));
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
        if (parallel.rank != 0) specification.name.resize(static_cast<std::size_t>(length));
        detail::checkMpi(MPI_Bcast(specification.name.data(), length, MPI_CHAR, 0,
                                   parallel.communicator), "MPI_Bcast(patch name)");
        specification.kind = static_cast<PatchKind>(kind);
    }
}

void broadcastMesh(Mesh& mesh, const ParallelContext& parallel) {
    Index counts[3]{};
    std::vector<PatchSpec> patches;
    std::vector<BoundaryFaceSpec> boundaries;
    if (parallel.rank == 0) {
        counts[0] = mesh.vertexCount();
        counts[1] = mesh.cellCount();
        for (Index patch = 0; patch < mesh.patchCount(); ++patch) {
            patches.push_back({mesh.patchName(patch), mesh.patchKind(patch)});
        }
        for (Index face = 0; face < mesh.faceCount(); ++face) {
            if (!mesh.boundaryFace(face)) continue;
            boundaries.push_back({mesh.faceVertices(face), detail::meshData(mesh).face_patch[face]});
        }
        counts[2] = static_cast<Index>(boundaries.size());
    }
    detail::checkMpi(MPI_Bcast(counts, 3, MPI_INT, 0, parallel.communicator),
                     "MPI_Bcast(mesh counts)");
    if (counts[0] <= 0 || counts[1] <= 0 || counts[2] <= 0) {
        throw std::runtime_error("distributed mesh counts are invalid");
    }
    broadcastPatch(patches, parallel);

    std::vector<double> coordinates(static_cast<std::size_t>(counts[0]) * 3U);
    std::vector<Index> connectivity(static_cast<std::size_t>(counts[1]) * 8U);
    std::vector<Index> boundary_vertices(static_cast<std::size_t>(counts[2]) * 4U);
    std::vector<Index> boundary_patches(static_cast<std::size_t>(counts[2]));
    if (parallel.rank == 0) {
        for (Index vertex = 0; vertex < mesh.vertexCount(); ++vertex) {
            const Vec3& point = mesh.vertex(vertex);
            coordinates[3U * static_cast<std::size_t>(vertex)] = point.x;
            coordinates[3U * static_cast<std::size_t>(vertex) + 1U] = point.y;
            coordinates[3U * static_cast<std::size_t>(vertex) + 2U] = point.z;
        }
        for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
            for (Index local = 0; local < 8; ++local) {
                connectivity[8U * static_cast<std::size_t>(cell) + static_cast<std::size_t>(local)] =
                    mesh.cellVertices(cell)[static_cast<std::size_t>(local)];
            }
        }
        for (std::size_t face = 0; face < boundaries.size(); ++face) {
            for (std::size_t local = 0; local < 4; ++local) {
                boundary_vertices[4U * face + local] = boundaries[face].vertices[local];
            }
            boundary_patches[face] = boundaries[face].patch;
        }
    }
    detail::checkMpi(MPI_Bcast(coordinates.data(), detail::mpiCount(coordinates.size(), "mesh coordinates"),
                               MPI_DOUBLE, 0, parallel.communicator), "MPI_Bcast(mesh coordinates)");
    detail::checkMpi(MPI_Bcast(connectivity.data(), detail::mpiCount(connectivity.size(), "mesh cells"),
                               MPI_INT, 0, parallel.communicator), "MPI_Bcast(mesh cells)");
    detail::checkMpi(MPI_Bcast(boundary_vertices.data(),
                               detail::mpiCount(boundary_vertices.size(), "mesh boundary vertices"),
                               MPI_INT, 0, parallel.communicator), "MPI_Bcast(mesh boundary vertices)");
    detail::checkMpi(MPI_Bcast(boundary_patches.data(),
                               detail::mpiCount(boundary_patches.size(), "mesh boundary patches"),
                               MPI_INT, 0, parallel.communicator), "MPI_Bcast(mesh boundary patches)");
    if (parallel.rank != 0) {
        std::vector<Vec3> vertices(static_cast<std::size_t>(counts[0]));
        std::vector<std::array<Index, 8>> cells(static_cast<std::size_t>(counts[1]));
        boundaries.resize(static_cast<std::size_t>(counts[2]));
        for (Index vertex = 0; vertex < counts[0]; ++vertex) {
            vertices[static_cast<std::size_t>(vertex)] = {
                coordinates[3U * static_cast<std::size_t>(vertex)],
                coordinates[3U * static_cast<std::size_t>(vertex) + 1U],
                coordinates[3U * static_cast<std::size_t>(vertex) + 2U]};
        }
        for (Index cell = 0; cell < counts[1]; ++cell) {
            for (Index local = 0; local < 8; ++local) {
                cells[static_cast<std::size_t>(cell)][static_cast<std::size_t>(local)] =
                    connectivity[8U * static_cast<std::size_t>(cell) + static_cast<std::size_t>(local)];
            }
        }
        for (std::size_t face = 0; face < boundaries.size(); ++face) {
            for (std::size_t local = 0; local < 4; ++local) {
                boundaries[face].vertices[local] = boundary_vertices[4U * face + local];
            }
            boundaries[face].patch = boundary_patches[face];
        }
        detail::MeshAccess::replace(mesh, Mesh::unstructured(std::move(vertices), std::move(cells),
                                                               std::move(patches), std::move(boundaries)));
    }
}

std::vector<Index> graphPartitionOwners(const Mesh& mesh, int partitions) {
    const Index cells = mesh.cellCount();
    if (partitions <= 0 || cells < partitions) {
        throw std::invalid_argument("graph partition has an invalid number of parts");
    }
    const auto& neighbours = detail::meshData(mesh).cell_neighbours;
    std::vector<Index> capacities(static_cast<std::size_t>(partitions));
    std::vector<Index> filled(static_cast<std::size_t>(partitions), 0);
    for (int part = 0; part < partitions; ++part) {
        capacities[static_cast<std::size_t>(part)] = cells / partitions +
            (part < cells % partitions ? 1 : 0);
    }

    // Select dispersed deterministic seeds using only graph distance.  The cell ID
    // decides ties, but never defines the partition boundary.
    const Index unreachable = std::numeric_limits<Index>::max();
    std::vector<Index> nearest(static_cast<std::size_t>(cells), unreachable);
    std::vector<Index> seeds;
    seeds.reserve(static_cast<std::size_t>(partitions));
    for (int part = 0; part < partitions; ++part) {
        Index seed = invalid_index;
        for (Index cell = 0; cell < cells; ++cell) {
            if (std::find(seeds.begin(), seeds.end(), cell) != seeds.end()) continue;
            if (seed == invalid_index ||
                nearest[static_cast<std::size_t>(cell)] > nearest[static_cast<std::size_t>(seed)] ||
                (nearest[static_cast<std::size_t>(cell)] == nearest[static_cast<std::size_t>(seed)] &&
                 detail::globalCellId(mesh, cell) < detail::globalCellId(mesh, seed))) {
                seed = cell;
            }
        }
        if (seed == invalid_index) throw std::logic_error("graph partition seed selection failed");
        seeds.push_back(seed);

        std::vector<Index> distance(static_cast<std::size_t>(cells), invalid_index);
        std::deque<Index> frontier{seed};
        distance[static_cast<std::size_t>(seed)] = 0;
        while (!frontier.empty()) {
            const Index cell = frontier.front();
            frontier.pop_front();
            for (Index neighbour : neighbours[static_cast<std::size_t>(cell)]) {
                if (neighbour == invalid_index ||
                    distance[static_cast<std::size_t>(neighbour)] != invalid_index) continue;
                distance[static_cast<std::size_t>(neighbour)] =
                    distance[static_cast<std::size_t>(cell)] + 1;
                frontier.push_back(neighbour);
            }
        }
        for (Index cell = 0; cell < cells; ++cell) {
            const Index path = distance[static_cast<std::size_t>(cell)];
            if (path != invalid_index) {
                nearest[static_cast<std::size_t>(cell)] = std::min(
                    nearest[static_cast<std::size_t>(cell)], path);
            }
        }
    }

    std::vector<Index> owners(static_cast<std::size_t>(cells), invalid_index);
    std::vector<std::deque<Index>> frontiers(static_cast<std::size_t>(partitions));
    for (int part = 0; part < partitions; ++part) {
        const Index seed = seeds[static_cast<std::size_t>(part)];
        owners[static_cast<std::size_t>(seed)] = part;
        ++filled[static_cast<std::size_t>(part)];
        frontiers[static_cast<std::size_t>(part)].push_back(seed);
    }

    Index remaining = cells - partitions;
    while (remaining > 0) {
        bool grew = false;
        for (int part = 0; part < partitions; ++part) {
            if (filled[static_cast<std::size_t>(part)] >= capacities[static_cast<std::size_t>(part)]) continue;
            bool assigned = false;
            std::deque<Index>& frontier = frontiers[static_cast<std::size_t>(part)];
            while (!frontier.empty() && !assigned) {
                const Index cell = frontier.front();
                frontier.pop_front();
                for (Index neighbour : neighbours[static_cast<std::size_t>(cell)]) {
                    if (neighbour == invalid_index ||
                        owners[static_cast<std::size_t>(neighbour)] != invalid_index) continue;
                    owners[static_cast<std::size_t>(neighbour)] = part;
                    ++filled[static_cast<std::size_t>(part)];
                    --remaining;
                    frontier.push_back(neighbour);
                    assigned = true;
                    grew = true;
                    break;
                }
            }
        }
        if (grew) continue;

        // Disconnected components have no frontier edge; seed the next component
        // in a non-full part, still without referring to any geometric axes.
        int part = 0;
        while (part < partitions &&
               filled[static_cast<std::size_t>(part)] >= capacities[static_cast<std::size_t>(part)]) {
            ++part;
        }
        Index seed = invalid_index;
        for (Index cell = 0; cell < cells; ++cell) {
            if (owners[static_cast<std::size_t>(cell)] == invalid_index &&
                (seed == invalid_index || detail::globalCellId(mesh, cell) < detail::globalCellId(mesh, seed))) {
                seed = cell;
            }
        }
        if (part == partitions || seed == invalid_index) {
            throw std::logic_error("graph partition growth failed");
        }
        owners[static_cast<std::size_t>(seed)] = part;
        ++filled[static_cast<std::size_t>(part)];
        --remaining;
        frontiers[static_cast<std::size_t>(part)].push_back(seed);
    }
    return owners;
}

Mesh partitionMesh(const Mesh& global, int rank, int size, Index ghost_layers) {
    if (size <= 0 || rank < 0 || rank >= size || ghost_layers < 3) {
        throw std::invalid_argument("mesh partition rank, size, or halo width is invalid");
    }
    global.validate();
    if (detail::ownedCellCount(global) != global.cellCount()) {
        throw std::invalid_argument("domain decomposition requires a complete mesh");
    }
    if (size == 1) return global;
    const Index global_cells = global.globalCellCount();
    if (global_cells < size) throw std::invalid_argument("MPI has more ranks than mesh cells");

    const std::vector<Index> owners = graphPartitionOwners(global, size);
    std::vector<Index> depth(static_cast<std::size_t>(global.cellCount()), invalid_index);
    std::deque<Index> frontier;
    for (Index cell = 0; cell < global.cellCount(); ++cell) {
        const Index owner = owners[static_cast<std::size_t>(cell)];
        if (owner == rank) {
            depth[static_cast<std::size_t>(cell)] = 0;
            frontier.push_back(cell);
        }
    }
    if (frontier.empty()) throw std::logic_error("partition owns no global cells");
    while (!frontier.empty()) {
        const Index cell = frontier.front();
        frontier.pop_front();
        const Index current_depth = depth[static_cast<std::size_t>(cell)];
        if (current_depth >= ghost_layers) continue;
        for (Index neighbour : detail::meshData(global).cell_neighbours[static_cast<std::size_t>(cell)]) {
            if (neighbour == invalid_index || depth[static_cast<std::size_t>(neighbour)] != invalid_index) continue;
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
    std::vector<Index> vertex_to_local(static_cast<std::size_t>(global.vertexCount()), invalid_index);
    std::vector<Vec3> vertices;
    std::vector<std::array<Index, 8>> cells;
    std::vector<Index> cell_ids;
    std::vector<Index> cell_owners;
    std::vector<Index> cell_depths;
    vertices.reserve(selected.size() * 4U);
    cells.reserve(selected.size());
    for (Index source : selected) {
        std::array<Index, 8> local_vertices{};
        for (Index local = 0; local < 8; ++local) {
            const Index original = global.cellVertices(source)[static_cast<std::size_t>(local)];
            Index& mapped = vertex_to_local[static_cast<std::size_t>(original)];
            if (mapped == invalid_index) {
                mapped = static_cast<Index>(vertices.size());
                vertices.push_back(global.vertex(original));
            }
            local_vertices[static_cast<std::size_t>(local)] = mapped;
        }
        cells.push_back(local_vertices);
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
        const bool neighbour_inside = neighbour != invalid_index &&
            source_to_local[static_cast<std::size_t>(neighbour)] != invalid_index;
        if (owner_inside == neighbour_inside) continue;
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
    std::vector<BoundaryFaceSpec> boundaries;
    std::map<Quad, Index> global_face_by_local_vertices;
    for (Index face = 0; face < global.faceCount(); ++face) {
        const Index patch = face_patch[static_cast<std::size_t>(face)];
        const Index owner = global.owner(face);
        const Index neighbour = global.neighbour(face);
        const bool attached = source_to_local[static_cast<std::size_t>(owner)] != invalid_index ||
            (neighbour != invalid_index && source_to_local[static_cast<std::size_t>(neighbour)] != invalid_index);
        if (!attached) continue;
        Quad local_vertices{};
        for (std::size_t local = 0; local < local_vertices.size(); ++local) {
            const Index original = global.faceVertices(face)[local];
            local_vertices[local] = vertex_to_local[static_cast<std::size_t>(original)];
            if (local_vertices[local] == invalid_index) {
                throw std::logic_error("local mesh omitted a face vertex");
            }
        }
        global_face_by_local_vertices.emplace(canonical(local_vertices), face);
        if (patch != invalid_index) boundaries.push_back({local_vertices, patch});
    }
    Mesh local = Mesh::unstructured(std::move(vertices), std::move(cells), std::move(patches),
                                    std::move(boundaries));
    std::vector<Index> face_ids(static_cast<std::size_t>(local.faceCount()));
    std::vector<Index> face_owners(static_cast<std::size_t>(local.faceCount()));
    for (Index face = 0; face < local.faceCount(); ++face) {
        const auto original = global_face_by_local_vertices.find(canonical(local.faceVertices(face)));
        if (original == global_face_by_local_vertices.end()) {
            throw std::logic_error("local face has no source global face");
        }
        face_ids[static_cast<std::size_t>(face)] = detail::globalFaceId(global, original->second);
        face_owners[static_cast<std::size_t>(face)] = owners[static_cast<std::size_t>(global.owner(original->second))];
    }
    detail::MeshAccess::setPartition(local, global_cells, ghost_layers, std::move(cell_ids),
                                     std::move(cell_owners), std::move(cell_depths),
                                     std::move(face_ids), std::move(face_owners), rank);
    local.validate();
    return local;
}

}  // namespace

ParallelContext ParallelContext::world(MPI_Comm communicator_value) {
    requireCommunicator(communicator_value, "ParallelContext::world");
    ParallelContext result;
    result.communicator = communicator_value;
    detail::checkMpi(MPI_Comm_rank(communicator_value, &result.rank), "MPI_Comm_rank");
    detail::checkMpi(MPI_Comm_size(communicator_value, &result.size), "MPI_Comm_size");
    result.validate();
    return result;
}

void ParallelContext::validate() const {
    if (size <= 0 || rank < 0 || rank >= size) {
        throw std::invalid_argument("parallel communicator, rank, or size is invalid");
    }
    if (communicator == MPI_COMM_NULL) {
        if (size != 1 || rank != 0) {
            throw std::invalid_argument("distributed context cannot use MPI_COMM_NULL");
        }
        return;
    }
    requireCommunicator(communicator, "ParallelContext::validate");
    int actual_rank = -1;
    int actual_size = 0;
    detail::checkMpi(MPI_Comm_rank(communicator, &actual_rank), "MPI_Comm_rank");
    detail::checkMpi(MPI_Comm_size(communicator, &actual_size), "MPI_Comm_size");
    if (actual_rank != rank || actual_size != size) {
        throw std::logic_error("parallel context rank/size no longer matches communicator");
    }
}

void ParallelContext::sum(const double* local, double* global, int count) const {
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel sum buffer is invalid");
    }
    requireCollectiveContext(*this, "ParallelContext::sum");
    if (count == 0) return;
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    detail::checkMpi(MPI_Allreduce(local, global, count, MPI_DOUBLE, MPI_SUM, communicator),
                     "MPI_Allreduce(sum)");
}

void ParallelContext::sum(const int* local, int* global, int count) const {
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel integer sum buffer is invalid");
    }
    requireCollectiveContext(*this, "ParallelContext::sum int array");
    if (count == 0) return;
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    detail::checkMpi(MPI_Allreduce(local, global, count, MPI_INT, MPI_SUM, communicator),
                     "MPI_Allreduce(sum int array)");
}

void ParallelContext::maximum(const double* local, double* global, int count) const {
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel maximum buffer is invalid");
    }
    requireCollectiveContext(*this, "ParallelContext::maximum");
    if (count == 0) return;
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    detail::checkMpi(MPI_Allreduce(local, global, count, MPI_DOUBLE, MPI_MAX, communicator),
                     "MPI_Allreduce(maximum)");
}

int ParallelContext::sum(int local) const {
    int global = 0;
    sum(&local, &global, 1);
    return global;
}

int ParallelContext::maximum(int local) const {
    requireCollectiveContext(*this, "ParallelContext::maximum int");
    if (!distributed()) return local;
    int global = 0;
    detail::checkMpi(MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, communicator),
                     "MPI_Allreduce(maximum int)");
    return global;
}

void ParallelContext::barrier() const {
    validate();
    if (distributed()) detail::checkMpi(MPI_Barrier(communicator), "MPI_Barrier");
}

Mesh decompose(const Mesh& global, const ParallelContext& parallel, Index ghost_layers) {
    parallel.validate();
    return partitionMesh(global, parallel.rank, parallel.size, ghost_layers);
}

Mesh readDistributedMesh(
    const std::filesystem::path& path,
    const ParallelContext& parallel,
    Index ghost_layers)
{
    parallel.validate();
    if (ghost_layers < 3) throw std::invalid_argument("ghostLayers must be >= 3");
    if (!parallel.distributed()) return readMeshFile(path);
    Mesh global;
    std::string error;
    int ok = 1;
    if (parallel.rank == 0) {
        try {
            detail::MeshAccess::replace(global, readMeshFile(path));
            global.validate();
        } catch (const std::exception& exception) {
            ok = 0;
            error = exception.what();
        }
    }
    detail::checkMpi(MPI_Bcast(&ok, 1, MPI_INT, 0, parallel.communicator),
                     "MPI_Bcast(mesh read status)");
    int error_length = parallel.rank == 0 ? static_cast<int>(error.size()) : 0;
    detail::checkMpi(MPI_Bcast(&error_length, 1, MPI_INT, 0, parallel.communicator),
                     "MPI_Bcast(mesh read error length)");
    if (error_length < 0 || error_length > 16384) {
        throw std::runtime_error("distributed mesh read error message is invalid");
    }
    if (parallel.rank != 0) error.resize(static_cast<std::size_t>(error_length));
    if (error_length > 0) {
        detail::checkMpi(MPI_Bcast(error.data(), error_length, MPI_CHAR, 0, parallel.communicator),
                         "MPI_Bcast(mesh read error)");
    }
    if (ok == 0) throw std::runtime_error("distributed mesh read failed: " + error);
    broadcastMesh(global, parallel);
    Mesh local = partitionMesh(global, parallel.rank, parallel.size, ghost_layers);
    local.validate();
    return local;
}

HaloExchange::HaloExchange(const Mesh& mesh, ParallelContext parallel)
    : m_mesh(&mesh), m_parallel(parallel)
{
    m_parallel.validate();
    mesh.validate();
    if (!m_parallel.distributed()) return;
    const auto build = [&](ExchangePlan& plan, bool faces, bool first_layer) {
        std::vector<std::vector<Index>> requested_ids(static_cast<std::size_t>(m_parallel.size));
        std::vector<std::vector<Index>> requested_indices(static_cast<std::size_t>(m_parallel.size));
        const Index count = faces ? mesh.faceCount() : mesh.cellCount();
        for (Index entity = 0; entity < count; ++entity) {
            const Index owner = faces ? detail::faceOwnerRank(mesh, entity) :
                detail::cellOwnerRank(mesh, entity);
            if (owner == m_parallel.rank) continue;
            if (!faces && first_layer &&
                detail::meshData(mesh).cell_ghost_depths[static_cast<std::size_t>(entity)] != 1) continue;
            if (!faces && !first_layer && detail::isOwned(mesh, entity)) continue;
            if (owner < 0 || owner >= m_parallel.size) {
                throw std::invalid_argument("halo entity has an invalid owner rank");
            }
            requested_ids[static_cast<std::size_t>(owner)].push_back(
                faces ? detail::globalFaceId(mesh, entity) : detail::globalCellId(mesh, entity));
            requested_indices[static_cast<std::size_t>(owner)].push_back(entity);
        }
        plan.receive_counts.assign(static_cast<std::size_t>(m_parallel.size), 0);
        plan.receive_offsets.assign(static_cast<std::size_t>(m_parallel.size), 0);
        std::vector<Index> outgoing_ids;
        for (int peer = 0; peer < m_parallel.size; ++peer) {
            plan.receive_offsets[static_cast<std::size_t>(peer)] =
                detail::mpiCount(outgoing_ids.size(), "halo request offset");
            plan.receive_counts[static_cast<std::size_t>(peer)] = detail::mpiCount(
                requested_ids[static_cast<std::size_t>(peer)].size(), "halo request count");
            outgoing_ids.insert(outgoing_ids.end(), requested_ids[static_cast<std::size_t>(peer)].begin(),
                                requested_ids[static_cast<std::size_t>(peer)].end());
            plan.receive_indices.insert(plan.receive_indices.end(),
                                        requested_indices[static_cast<std::size_t>(peer)].begin(),
                                        requested_indices[static_cast<std::size_t>(peer)].end());
        }
        plan.send_counts.assign(static_cast<std::size_t>(m_parallel.size), 0);
        detail::checkMpi(MPI_Alltoall(plan.receive_counts.data(), 1, MPI_INT, plan.send_counts.data(), 1,
                                      MPI_INT, m_parallel.communicator), "MPI_Alltoall(halo request counts)");
        plan.send_offsets.assign(static_cast<std::size_t>(m_parallel.size), 0);
        std::size_t incoming_size = 0;
        for (int peer = 0; peer < m_parallel.size; ++peer) {
            plan.send_offsets[static_cast<std::size_t>(peer)] = detail::mpiCount(incoming_size, "halo request offset");
            incoming_size += static_cast<std::size_t>(plan.send_counts[static_cast<std::size_t>(peer)]);
        }
        std::vector<Index> incoming_ids(incoming_size);
        Index dummy = 0;
        detail::checkMpi(MPI_Alltoallv(
            outgoing_ids.empty() ? &dummy : outgoing_ids.data(), plan.receive_counts.data(), plan.receive_offsets.data(),
            MPI_INT, incoming_ids.empty() ? &dummy : incoming_ids.data(), plan.send_counts.data(), plan.send_offsets.data(),
            MPI_INT, m_parallel.communicator), "MPI_Alltoallv(halo requests)");
        std::map<Index, Index> owned;
        for (Index entity = 0; entity < count; ++entity) {
            const bool owner = faces ? detail::faceOwnerRank(mesh, entity) == m_parallel.rank :
                detail::isOwned(mesh, entity);
            if (owner) owned.emplace(faces ? detail::globalFaceId(mesh, entity) :
                                      detail::globalCellId(mesh, entity), entity);
        }
        for (Index id : incoming_ids) {
            const auto local = owned.find(id);
            if (local == owned.end()) throw std::runtime_error("halo request is not owned by this rank");
            plan.send_indices.push_back(local->second);
        }
    };
    build(m_cells, false, false);
    build(m_first_layer_cells, false, true);
    build(m_faces, true, false);
}

void HaloExchange::exchange(double* values, std::size_t components, ExchangePlan& plan) {
    if (!m_parallel.distributed()) return;
    if (values == nullptr || components == 0) throw std::invalid_argument("halo exchange values are invalid");
    plan.send_buffer.resize(plan.send_indices.size() * components);
    plan.receive_buffer.resize(plan.receive_indices.size() * components);
    for (std::size_t index = 0; index < plan.send_indices.size(); ++index) {
        for (std::size_t component = 0; component < components; ++component) {
            plan.send_buffer[index * components + component] =
                values[static_cast<std::size_t>(plan.send_indices[index]) * components + component];
        }
    }
    std::vector<int> send_counts(static_cast<std::size_t>(m_parallel.size));
    std::vector<int> send_offsets(static_cast<std::size_t>(m_parallel.size));
    std::vector<int> receive_counts(static_cast<std::size_t>(m_parallel.size));
    std::vector<int> receive_offsets(static_cast<std::size_t>(m_parallel.size));
    for (int peer = 0; peer < m_parallel.size; ++peer) {
        const std::size_t p = static_cast<std::size_t>(peer);
        send_counts[p] = detail::mpiCount(static_cast<std::size_t>(plan.send_counts[p]) * components,
                                           "halo send values");
        send_offsets[p] = detail::mpiCount(static_cast<std::size_t>(plan.send_offsets[p]) * components,
                                            "halo send offset");
        receive_counts[p] = detail::mpiCount(static_cast<std::size_t>(plan.receive_counts[p]) * components,
                                              "halo receive values");
        receive_offsets[p] = detail::mpiCount(static_cast<std::size_t>(plan.receive_offsets[p]) * components,
                                               "halo receive offset");
    }
    double dummy = 0.0;
    detail::checkMpi(MPI_Alltoallv(
        plan.send_buffer.empty() ? &dummy : plan.send_buffer.data(), send_counts.data(), send_offsets.data(), MPI_DOUBLE,
        plan.receive_buffer.empty() ? &dummy : plan.receive_buffer.data(), receive_counts.data(), receive_offsets.data(),
        MPI_DOUBLE, m_parallel.communicator), "MPI_Alltoallv(halo values)");
    for (std::size_t index = 0; index < plan.receive_indices.size(); ++index) {
        for (std::size_t component = 0; component < components; ++component) {
            values[static_cast<std::size_t>(plan.receive_indices[index]) * components + component] =
                plan.receive_buffer[index * components + component];
        }
    }
}

void HaloExchange::exchange(double* values, std::size_t components) {
    exchange(values, components, m_cells);
}

void HaloExchange::exchangeFaces(double* values, std::size_t components) {
    exchange(values, components, m_faces);
}

void HaloExchange::exchange(std::vector<double>& values) {
    if (m_mesh == nullptr || values.size() != static_cast<std::size_t>(m_mesh->cellCount())) {
        throw std::invalid_argument("raw halo field has the wrong size");
    }
    exchange(values.data(), 1);
}

void HaloExchange::exchangeFirstLayer(std::vector<double>& values) {
    if (m_mesh == nullptr || values.size() != static_cast<std::size_t>(m_mesh->cellCount())) {
        throw std::invalid_argument("raw first-layer halo field has the wrong size");
    }
    exchange(values.data(), 1, m_first_layer_cells);
}

void HaloExchange::exchange(ScalarField& field) {
    if (&field.mesh() != m_mesh) throw std::invalid_argument("scalar halo field is incompatible");
    field.validateStorage();
    if (field.location() == FieldLocation::Cell) exchange(detail::fieldData(field), 1);
    else if (field.location() == FieldLocation::Face) exchangeFaces(detail::fieldData(field), 1);
    else throw std::invalid_argument("vertex scalar halo exchange is not supported");
}

void HaloExchange::exchange(VectorField& field) {
    if (&field.mesh() != m_mesh) throw std::invalid_argument("vector halo field is incompatible");
    field.validateStorage();
    if (field.location() == FieldLocation::Cell) exchange(&detail::fieldData(field)->x, 3);
    else if (field.location() == FieldLocation::Face) exchangeFaces(&detail::fieldData(field)->x, 3);
    else throw std::invalid_argument("vertex vector halo exchange is not supported");
}

void HaloExchange::exchange(TensorField& field) {
    if (&field.mesh() != m_mesh) throw std::invalid_argument("tensor halo field is incompatible");
    field.validateStorage();
    if (field.location() == FieldLocation::Cell) exchange(&detail::fieldData(field)->rows[0].x, 9);
    else if (field.location() == FieldLocation::Face) exchangeFaces(&detail::fieldData(field)->rows[0].x, 9);
    else throw std::invalid_argument("vertex tensor halo exchange is not supported");
}

}  // namespace babelsim
