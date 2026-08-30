#include "babelsim/parallel.h"

#include "babelsim/mesh_io.h"
#include "babelsim/mpi_support.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <stdexcept>
#include <utility>

namespace babelsim {
namespace {

std::array<PatchSpec, 6> temporaryPatches() {
    auto patches = defaultPatches();
    for (std::size_t side = 0; side < patches.size(); ++side) {
        patches[side].name = "local_side_" + std::to_string(side);
    }
    return patches;
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
    // 热路径只检查 MPI 生命周期和空通信器；rank/size 与通信器的一致性在
    // world()/Solver/HaloExchange 构造时完整校验，避免每个点积重复查询。
    detail::requireMpiReady(operation);
}

void assignPartitionPatches(
    Mesh& local,
    const std::array<PatchSpec, 6>& specifications,
    int rank,
    int size)
{
    std::vector<BoundaryPatch> patches;
    patches.reserve(8);
    for (const PatchSpec& specification : specifications) {
        patches.push_back({specification.name, specification.kind, {}});
    }
    const Index left_processor_patch = rank == 0
        ? invalid_index : static_cast<Index>(patches.size());
    if (left_processor_patch != invalid_index) {
        patches.push_back({"processor_left", PatchKind::Processor, {}});
    }
    const Index right_processor_patch = rank + 1 == size
        ? invalid_index : static_cast<Index>(patches.size());
    if (right_processor_patch != invalid_index) {
        patches.push_back({"processor_right", PatchKind::Processor, {}});
    }
    for (Index face = 0; face < local.faceCount(); ++face) {
        if (!local.boundaryFace(face)) continue;
        const auto f = static_cast<std::size_t>(face);
        const Side side = static_cast<Side>(local.face_patch[f]);
        Index patch = static_cast<Index>(side);
        if (side == Side::XMin && left_processor_patch != invalid_index) {
            patch = left_processor_patch;
        } else if (side == Side::XMax && right_processor_patch != invalid_index) {
            patch = right_processor_patch;
        }
        if (patch < 0 || static_cast<std::size_t>(patch) >= patches.size()) {
            throw std::logic_error("decomposed boundary face has no patch");
        }
        local.face_patch[f] = patch;
    }
    local.setPatches(std::move(patches));
    for (Index face = 0; face < local.faceCount(); ++face) {
        if (local.boundaryFace(face)) {
            const Index patch = local.face_patch[static_cast<std::size_t>(face)];
            local.addPatchFace(patch, face);
        }
    }
}

Mesh partitionMesh(const Mesh& global, int rank, int size, Index ghost_layers) {
    if (size <= 0 || rank < 0 || rank >= size || ghost_layers < 1) {
        throw std::invalid_argument("mesh partition rank, size, or halo width is invalid");
    }
    if (global.ownedCellCount() != global.cellCount()) {
        throw std::invalid_argument("domain decomposition requires a global mesh");
    }
    if (size == 1) return global;
    const Index global_nx = global.dimensions[0];
    if (static_cast<std::int64_t>(global_nx) <
        static_cast<std::int64_t>(ghost_layers) * size) {
        throw std::invalid_argument("each MPI partition needs at least ghost_layers owned x cells");
    }
    std::vector<Index> widths(static_cast<std::size_t>(size));
    Index remaining = global_nx;
    for (int part = 0; part < size; ++part) {
        widths[static_cast<std::size_t>(part)] =
            remaining / static_cast<Index>(size - part);
        remaining -= widths[static_cast<std::size_t>(part)];
    }
    if (std::any_of(widths.begin(), widths.end(), [ghost_layers](Index width) {
            return width < ghost_layers;
        })) {
        throw std::invalid_argument("an MPI partition is narrower than its halo");
    }
    Index owned_global_begin = 0;
    for (int part = 0; part < rank; ++part) {
        owned_global_begin += widths[static_cast<std::size_t>(part)];
    }
    const Index owned_width = widths[static_cast<std::size_t>(rank)];
    const Index left_ghost = rank == 0 ? 0 : ghost_layers;
    const Index right_ghost = rank + 1 == size ? 0 : ghost_layers;
    const Index global_offset = owned_global_begin - left_ghost;
    const Index local_nx = left_ghost + owned_width + right_ghost;
    std::vector<Vec3> points;
    points.reserve(
        static_cast<std::size_t>(local_nx + 1) *
        static_cast<std::size_t>(global.dimensions[1] + 1) *
        static_cast<std::size_t>(global.dimensions[2] + 1));
    for (Index k = 0; k <= global.dimensions[2]; ++k) {
        for (Index j = 0; j <= global.dimensions[1]; ++j) {
            for (Index i = 0; i <= local_nx; ++i) {
                points.push_back(global.vertices[static_cast<std::size_t>(
                    global.vertexId(global_offset + i, j, k))]);
            }
        }
    }
    std::array<PatchSpec, 6> specifications{};
    if (global.patches.size() != specifications.size()) {
        throw std::invalid_argument("global structured mesh must have six logical patches");
    }
    for (std::size_t side = 0; side < specifications.size(); ++side) {
        specifications[side] = {
            global.patches[side].name, global.patches[side].kind};
    }
    Mesh local = Mesh::structured(
        {local_nx, global.dimensions[1], global.dimensions[2]},
        std::move(points), temporaryPatches());
    local.setOwnership(
        global.dimensions, global_offset, left_ghost,
        left_ghost + owned_width, ghost_layers);
    assignPartitionPatches(local, specifications, rank, size);
    local.validate();
    return local;
}

void broadcastPatchSpecifications(
    std::array<PatchSpec, 6>& specifications,
    const ParallelContext& parallel)
{
    for (PatchSpec& specification : specifications) {
        int kind = static_cast<int>(specification.kind);
        if (parallel.rank == 0 &&
            (specification.name.empty() || specification.name.size() > 4096)) {
            throw std::runtime_error("distributed mesh patch name is too long");
        }
        int name_length = static_cast<int>(specification.name.size());
        detail::checkMpi(
            MPI_Bcast(&kind, 1, MPI_INT, 0, parallel.communicator), "MPI_Bcast(patch kind)");
        detail::checkMpi(
            MPI_Bcast(&name_length, 1, MPI_INT, 0, parallel.communicator),
            "MPI_Bcast(patch name length)");
        if (name_length <= 0 || name_length > 4096) {
            throw std::runtime_error("distributed mesh contains an invalid patch name length");
        }
        if (parallel.rank != 0) specification.name.resize(static_cast<std::size_t>(name_length));
        detail::checkMpi(
            MPI_Bcast(specification.name.data(), name_length, MPI_CHAR, 0, parallel.communicator),
            "MPI_Bcast(patch name)");
        specification.kind = static_cast<PatchKind>(kind);
    }
}

void sendPartition(const Mesh& local, int destination, MPI_Comm communicator) {
    const int metadata[5] = {
        local.dimensions[0], local.global_i_offset, local.owned_i_begin,
        local.owned_i_end, local.ghost_layers};
    detail::checkMpi(
        MPI_Send(metadata, 5, MPI_INT, destination, 710, communicator),
        "MPI_Send(mesh metadata)");
    std::vector<double> coordinates;
    coordinates.reserve(local.vertices.size() * 3U);
    for (const Vec3& point : local.vertices) {
        coordinates.push_back(point.x);
        coordinates.push_back(point.y);
        coordinates.push_back(point.z);
    }
    detail::checkMpi(
        MPI_Send(
            coordinates.data(),
            detail::mpiCount(coordinates.size(), "distributed mesh coordinates"),
            MPI_DOUBLE, destination, 711, communicator),
        "MPI_Send(mesh coordinates)");
}

Mesh receivePartition(
    const std::array<Index, 3>& global_dimensions,
    const std::array<PatchSpec, 6>& specifications,
    int source,
    const ParallelContext& parallel)
{
    int metadata[5]{};
    detail::checkMpi(
        MPI_Recv(metadata, 5, MPI_INT, source, 710, parallel.communicator,
                 MPI_STATUS_IGNORE),
        "MPI_Recv(mesh metadata)");
    const Index local_nx = metadata[0];
    const Index global_offset = metadata[1];
    const Index owned_begin = metadata[2];
    const Index owned_end = metadata[3];
    const Index ghost_layers = metadata[4];
    if (local_nx <= 0 || global_offset < 0 || owned_begin < 0 ||
        owned_end <= owned_begin || owned_end > local_nx || ghost_layers < 1 ||
        global_offset + local_nx > global_dimensions[0]) {
        throw std::runtime_error("received distributed mesh metadata is invalid");
    }
    const std::size_t point_count =
        static_cast<std::size_t>(local_nx + 1) *
        static_cast<std::size_t>(global_dimensions[1] + 1) *
        static_cast<std::size_t>(global_dimensions[2] + 1);
    std::vector<double> coordinates(point_count * 3U);
    detail::checkMpi(
        MPI_Recv(
            coordinates.data(), detail::mpiCount(coordinates.size(), "distributed mesh coordinates"),
            MPI_DOUBLE, source, 711, parallel.communicator, MPI_STATUS_IGNORE),
        "MPI_Recv(mesh coordinates)");
    std::vector<Vec3> points(point_count);
    for (std::size_t point = 0; point < point_count; ++point) {
        points[point] = {
            coordinates[3U * point], coordinates[3U * point + 1U],
            coordinates[3U * point + 2U]};
    }
    Mesh local = Mesh::structured(
        {local_nx, global_dimensions[1], global_dimensions[2]},
        std::move(points), temporaryPatches());
    local.setOwnership(
        global_dimensions, global_offset, owned_begin, owned_end, ghost_layers);
    assignPartitionPatches(local, specifications, parallel.rank, parallel.size);
    local.validate();
    return local;
}

}  // 匿名命名空间

ParallelContext ParallelContext::world(MPI_Comm communicator_value) {
    requireCommunicator(communicator_value, "ParallelContext::world");
    ParallelContext result;
    result.communicator = communicator_value;
    detail::checkMpi(
        MPI_Comm_rank(communicator_value, &result.rank), "MPI_Comm_rank");
    detail::checkMpi(
        MPI_Comm_size(communicator_value, &result.size), "MPI_Comm_size");
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

void ParallelContext::sum(
    const double* local,
    double* global,
    int count) const
{
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel sum buffer is invalid");
    }
    requireCollectiveContext(*this, "ParallelContext::sum");
    if (count == 0) {
        return;
    }
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    detail::checkMpi(
        MPI_Allreduce(local, global, count, MPI_DOUBLE, MPI_SUM, communicator),
        "MPI_Allreduce(sum)");
}

void ParallelContext::sum(
    const int* local,
    int* global,
    int count) const
{
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel integer sum buffer is invalid");
    }
    requireCollectiveContext(*this, "ParallelContext::sum int array");
    if (count == 0) {
        return;
    }
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    detail::checkMpi(
        MPI_Allreduce(local, global, count, MPI_INT, MPI_SUM, communicator),
        "MPI_Allreduce(sum int array)");
}

void ParallelContext::maximum(
    const double* local,
    double* global,
    int count) const
{
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel maximum buffer is invalid");
    }
    requireCollectiveContext(*this, "ParallelContext::maximum");
    if (count == 0) {
        return;
    }
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    detail::checkMpi(
        MPI_Allreduce(local, global, count, MPI_DOUBLE, MPI_MAX, communicator),
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
    detail::checkMpi(
        MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, communicator),
        "MPI_Allreduce(maximum int)");
    return global;
}

void ParallelContext::barrier() const {
    validate();
    if (distributed()) {
        detail::checkMpi(MPI_Barrier(communicator), "MPI_Barrier");
    }
}

Mesh decompose(
    const Mesh& global,
    const ParallelContext& parallel,
    Index ghost_layers)
{
    parallel.validate();
    global.validate();
    return partitionMesh(global, parallel.rank, parallel.size, ghost_layers);
}

Mesh readDistributedMesh(
    const std::filesystem::path& path,
    const ParallelContext& parallel,
    Index ghost_layers)
{
    parallel.validate();
    if (!parallel.distributed()) {
        return readMeshFile(path);
    }

    Mesh global;
    std::string error;
    int ok = 1;
    if (parallel.rank == 0) {
        try {
            global = readMeshFile(path);
            global.validate();
            if (global.patches.size() != 6) {
                throw std::invalid_argument("distributed mesh requires six logical patches");
            }
            for (const BoundaryPatch& patch : global.patches) {
                if (patch.name.empty() || patch.name.size() > 4096) {
                    throw std::invalid_argument("distributed mesh patch name is too long");
                }
            }
        } catch (const std::exception& exception) {
            ok = 0;
            error = exception.what();
        }
    }
    detail::checkMpi(
        MPI_Bcast(&ok, 1, MPI_INT, 0, parallel.communicator),
        "MPI_Bcast(mesh read status)");
    int error_length = parallel.rank == 0 ? static_cast<int>(error.size()) : 0;
    detail::checkMpi(
        MPI_Bcast(&error_length, 1, MPI_INT, 0, parallel.communicator),
        "MPI_Bcast(mesh read error length)");
    if (error_length < 0 || error_length > 16384) {
        throw std::runtime_error("distributed mesh read error message is invalid");
    }
    if (parallel.rank != 0) error.resize(static_cast<std::size_t>(error_length));
    if (error_length > 0) {
        detail::checkMpi(
            MPI_Bcast(error.data(), error_length, MPI_CHAR, 0, parallel.communicator),
            "MPI_Bcast(mesh read error)");
    }
    if (ok == 0) {
        throw std::runtime_error("distributed mesh read failed: " + error);
    }

    std::array<Index, 3> dimensions{};
    std::array<PatchSpec, 6> specifications{};
    if (parallel.rank == 0) {
        dimensions = global.dimensions;
        for (std::size_t side = 0; side < specifications.size(); ++side) {
            specifications[side] = {
                global.patches[side].name, global.patches[side].kind};
        }
    }
    detail::checkMpi(
        MPI_Bcast(dimensions.data(), 3, MPI_INT, 0, parallel.communicator),
        "MPI_Bcast(mesh dimensions)");
    broadcastPatchSpecifications(specifications, parallel);

    Mesh local;
    if (parallel.rank == 0) {
        local = partitionMesh(global, 0, parallel.size, ghost_layers);
        for (int destination = 1; destination < parallel.size; ++destination) {
            sendPartition(
                partitionMesh(global, destination, parallel.size, ghost_layers),
                destination, parallel.communicator);
        }
    } else {
        local = receivePartition(dimensions, specifications, 0, parallel);
    }
    local.validate();
    return local;
}

HaloExchange::HaloExchange(const Mesh& mesh, ParallelContext parallel)
    : m_mesh(&mesh), m_parallel(parallel)
{
    m_parallel.validate();
    mesh.validate();
    if (!m_parallel.distributed()) {
        return;
    }
    if (mesh.ghost_layers < 1 ||
        mesh.owned_i_begin < (m_parallel.rank == 0 ? 0 : mesh.ghost_layers) ||
        mesh.dimensions[0] - mesh.owned_i_end <
            (m_parallel.rank + 1 == m_parallel.size ? 0 : mesh.ghost_layers)) {
        throw std::invalid_argument("mesh halo layout does not match MPI partition");
    }
    m_left = m_parallel.rank == 0 ? MPI_PROC_NULL : m_parallel.rank - 1;
    m_right = m_parallel.rank + 1 == m_parallel.size
        ? MPI_PROC_NULL : m_parallel.rank + 1;
    const auto appendPlane = [&](std::vector<Index>& indices, Index begin) {
        for (Index k = 0; k < mesh.dimensions[2]; ++k) {
            for (Index j = 0; j < mesh.dimensions[1]; ++j) {
                for (Index layer = 0; layer < mesh.ghost_layers; ++layer) {
                    indices.push_back(mesh.cellId(begin + layer, j, k));
                }
            }
        }
    };
    if (m_left != MPI_PROC_NULL) {
        appendPlane(m_send_left, mesh.owned_i_begin);
        appendPlane(m_receive_left, mesh.owned_i_begin - mesh.ghost_layers);
    }
    if (m_right != MPI_PROC_NULL) {
        appendPlane(m_send_right, mesh.owned_i_end - mesh.ghost_layers);
        appendPlane(m_receive_right, mesh.owned_i_end);
    }
    const auto appendFirstPlane = [&](std::vector<Index>& indices, Index begin) {
        for (Index k = 0; k < mesh.dimensions[2]; ++k) {
            for (Index j = 0; j < mesh.dimensions[1]; ++j) {
                indices.push_back(mesh.cellId(begin, j, k));
            }
        }
    };
    if (m_left != MPI_PROC_NULL) {
        appendFirstPlane(m_send_left_first, mesh.owned_i_begin);
        appendFirstPlane(m_receive_left_first, mesh.owned_i_begin - 1);
    }
    if (m_right != MPI_PROC_NULL) {
        appendFirstPlane(m_send_right_first, mesh.owned_i_end - 1);
        appendFirstPlane(m_receive_right_first, mesh.owned_i_end);
    }
    const auto appendInterface = [&](std::vector<Index>& indices, Index cell_i, Side side) {
        for (Index k = 0; k < mesh.dimensions[2]; ++k) {
            for (Index j = 0; j < mesh.dimensions[1]; ++j) {
                const Index cell = mesh.cellId(cell_i, j, k);
                indices.push_back(mesh.cell_faces[static_cast<std::size_t>(cell)]
                    [static_cast<std::size_t>(side)]);
            }
        }
    };
    if (m_left != MPI_PROC_NULL) {
        // 分区界面由较小 rank（左侧 owned cell）作为唯一发布者。
        appendInterface(m_receive_face_left, mesh.owned_i_begin, Side::XMin);
    }
    if (m_right != MPI_PROC_NULL) {
        appendInterface(m_send_face_right, mesh.owned_i_end - 1, Side::XMax);
    }
    // 常用场最多包含 9 个 double（Tensor3）。预留一次后，后续时间步/外迭代
    // 的打包和接收不会再次触发堆分配。
    m_send_buffer_left.reserve(m_send_left.size() * 9U);
    m_send_buffer_right.reserve(m_send_right.size() * 9U);
    m_receive_buffer_left.reserve(m_receive_left.size() * 9U);
    m_receive_buffer_right.reserve(m_receive_right.size() * 9U);
    m_send_face_buffer_right.reserve(m_send_face_right.size() * 9U);
    m_receive_face_buffer_left.reserve(m_receive_face_left.size() * 9U);
}

void HaloExchange::exchange(double* values, std::size_t components) {
    exchangeCells(
        values, components, m_send_left, m_send_right, m_receive_left, m_receive_right);
}

void HaloExchange::exchangeCells(
    double* values,
    std::size_t components,
    const std::vector<Index>& send_left,
    const std::vector<Index>& send_right,
    const std::vector<Index>& receive_left,
    const std::vector<Index>& receive_right)
{
    if (!m_parallel.distributed()) {
        return;
    }
    if (values == nullptr || components == 0) {
        throw std::invalid_argument("halo exchange values are invalid");
    }
    const auto pack = [values, components](
                          const std::vector<Index>& indices,
                          std::vector<double>& buffer) {
        buffer.resize(indices.size() * components);
        std::size_t output = 0;
        for (Index cell : indices) {
            const std::size_t begin = static_cast<std::size_t>(cell) * components;
            for (std::size_t component = 0; component < components; ++component) {
                buffer[output++] = values[begin + component];
            }
        }
    };
    const auto unpack = [values, components](
                            const std::vector<Index>& indices,
                            const std::vector<double>& buffer) {
        std::size_t input = 0;
        for (Index cell : indices) {
            const std::size_t begin = static_cast<std::size_t>(cell) * components;
            for (std::size_t component = 0; component < components; ++component) {
                values[begin + component] = buffer[input++];
            }
        }
    };

    pack(send_left, m_send_buffer_left);
    pack(send_right, m_send_buffer_right);
    m_receive_buffer_right.resize(receive_right.size() * components);
    m_receive_buffer_left.resize(receive_left.size() * components);
    double dummy = 0.0;
    detail::checkMpi(MPI_Sendrecv(
        m_send_buffer_left.empty() ? &dummy : m_send_buffer_left.data(),
        detail::mpiCount(m_send_buffer_left.size(), "left halo buffer"), MPI_DOUBLE, m_left, 101,
        m_receive_buffer_right.empty() ? &dummy : m_receive_buffer_right.data(),
        detail::mpiCount(m_receive_buffer_right.size(), "right halo buffer"), MPI_DOUBLE, m_right, 101,
        m_parallel.communicator, MPI_STATUS_IGNORE), "MPI_Sendrecv(left halo)");
    detail::checkMpi(MPI_Sendrecv(
        m_send_buffer_right.empty() ? &dummy : m_send_buffer_right.data(),
        detail::mpiCount(m_send_buffer_right.size(), "right halo buffer"), MPI_DOUBLE, m_right, 102,
        m_receive_buffer_left.empty() ? &dummy : m_receive_buffer_left.data(),
        detail::mpiCount(m_receive_buffer_left.size(), "left halo buffer"), MPI_DOUBLE, m_left, 102,
        m_parallel.communicator, MPI_STATUS_IGNORE), "MPI_Sendrecv(right halo)");
    unpack(receive_right, m_receive_buffer_right);
    unpack(receive_left, m_receive_buffer_left);
}

void HaloExchange::exchangeFaces(double* values, std::size_t components) {
    if (!m_parallel.distributed()) {
        return;
    }
    if (values == nullptr || components == 0) {
        throw std::invalid_argument("face halo exchange values are invalid");
    }
    const auto pack = [values, components](
                          const std::vector<Index>& indices,
                          std::vector<double>& buffer) {
        buffer.resize(indices.size() * components);
        std::size_t output = 0;
        for (Index face : indices) {
            const std::size_t begin = static_cast<std::size_t>(face) * components;
            for (std::size_t component = 0; component < components; ++component) {
                buffer[output++] = values[begin + component];
            }
        }
    };
    const auto unpack = [values, components](
                            const std::vector<Index>& indices,
                            const std::vector<double>& buffer) {
        std::size_t input = 0;
        for (Index face : indices) {
            const std::size_t begin = static_cast<std::size_t>(face) * components;
            for (std::size_t component = 0; component < components; ++component) {
                values[begin + component] = buffer[input++];
            }
        }
    };
    pack(m_send_face_right, m_send_face_buffer_right);
    m_receive_face_buffer_left.resize(m_receive_face_left.size() * components);
    double dummy = 0.0;
    // 每个 rank 只向右侧发布自己的 XMax 界面；同时从左侧接收 XMin
    // 界面。这样一个物理面始终由较小 rank 决定，不会出现交换振荡。
    detail::checkMpi(MPI_Sendrecv(
        m_send_face_buffer_right.empty() ? &dummy : m_send_face_buffer_right.data(),
        detail::mpiCount(m_send_face_buffer_right.size(), "right face halo buffer"), MPI_DOUBLE, m_right, 201,
        m_receive_face_buffer_left.empty() ? &dummy : m_receive_face_buffer_left.data(),
        detail::mpiCount(m_receive_face_buffer_left.size(), "left face halo buffer"), MPI_DOUBLE, m_left, 201,
        m_parallel.communicator, MPI_STATUS_IGNORE), "MPI_Sendrecv(face owner halo)");
    unpack(m_receive_face_left, m_receive_face_buffer_left);
}

void HaloExchange::exchange(std::vector<double>& values) {
    if (m_mesh == nullptr) throw std::logic_error("halo exchange has no mesh");
    if (values.size() != static_cast<std::size_t>(m_mesh->cellCount())) {
        throw std::invalid_argument("raw halo field has the wrong size");
    }
    exchange(values.data(), 1);
}

void HaloExchange::exchangeFirstLayer(std::vector<double>& values) {
    if (m_mesh == nullptr) throw std::logic_error("halo exchange has no mesh");
    if (values.size() != static_cast<std::size_t>(m_mesh->cellCount())) {
        throw std::invalid_argument("raw first-layer halo field has the wrong size");
    }
    exchangeCells(
        values.data(), 1, m_send_left_first, m_send_right_first,
        m_receive_left_first, m_receive_right_first);
}

void HaloExchange::exchange(ScalarField& field) {
    if (&field.mesh() != m_mesh) {
        throw std::invalid_argument("scalar halo field is incompatible");
    }
    field.validateStorage();
    if (field.location() == FieldLocation::Cell) {
        exchange(field.mutableData(), 1);
    } else if (field.location() == FieldLocation::Face) {
        exchangeFaces(field.mutableData(), 1);
    } else {
        throw std::invalid_argument("vertex scalar halo exchange is not supported");
    }
}

void HaloExchange::exchange(VectorField& field) {
    if (&field.mesh() != m_mesh) {
        throw std::invalid_argument("vector halo field is incompatible");
    }
    field.validateStorage();
    if (field.location() == FieldLocation::Cell) {
        exchange(&field.mutableData()->x, 3);
    } else if (field.location() == FieldLocation::Face) {
        exchangeFaces(&field.mutableData()->x, 3);
    } else {
        throw std::invalid_argument("vertex vector halo exchange is not supported");
    }
}

void HaloExchange::exchange(TensorField& field) {
    if (&field.mesh() != m_mesh) {
        throw std::invalid_argument("tensor halo field is incompatible");
    }
    field.validateStorage();
    if (field.location() == FieldLocation::Cell) {
        exchange(&field.mutableData()->rows[0].x, 9);
    } else if (field.location() == FieldLocation::Face) {
        exchangeFaces(&field.mutableData()->rows[0].x, 9);
    } else {
        throw std::invalid_argument("vertex tensor halo exchange is not supported");
    }
}

}  // babelsim 命名空间
