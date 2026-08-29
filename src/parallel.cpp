#include "babelsim/parallel.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <utility>

namespace babelsim {
namespace {

void requireMpiInitialized() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized == 0) {
        throw std::logic_error("MPI must be initialized before creating a world context");
    }
}

int mpiCount(std::size_t count) {
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("MPI halo buffer exceeds the MPI int count limit");
    }
    return static_cast<int>(count);
}

std::array<PatchSpec, 6> temporaryPatches() {
    auto patches = defaultPatches();
    for (std::size_t side = 0; side < patches.size(); ++side) {
        patches[side].name = "local_side_" + std::to_string(side);
    }
    return patches;
}

}  // namespace

ParallelContext ParallelContext::world(MPI_Comm communicator_value) {
    requireMpiInitialized();
    ParallelContext result;
    result.communicator = communicator_value;
    MPI_Comm_rank(communicator_value, &result.rank);
    MPI_Comm_size(communicator_value, &result.size);
    result.validate();
    return result;
}

void ParallelContext::validate() const {
    if (size <= 0 || rank < 0 || rank >= size ||
        (size > 1 && communicator == MPI_COMM_NULL)) {
        throw std::invalid_argument("parallel communicator, rank, or size is invalid");
    }
}

void ParallelContext::sum(
    const double* local,
    double* global,
    int count) const
{
    validate();
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel sum buffer is invalid");
    }
    if (count == 0) {
        return;
    }
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    MPI_Allreduce(local, global, count, MPI_DOUBLE, MPI_SUM, communicator);
}

void ParallelContext::maximum(
    const double* local,
    double* global,
    int count) const
{
    validate();
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel maximum buffer is invalid");
    }
    if (count == 0) {
        return;
    }
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    MPI_Allreduce(local, global, count, MPI_DOUBLE, MPI_MAX, communicator);
}

int ParallelContext::sum(int local) const {
    validate();
    if (!distributed()) {
        return local;
    }
    int global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, communicator);
    return global;
}

int ParallelContext::maximum(int local) const {
    validate();
    if (!distributed()) {
        return local;
    }
    int global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, communicator);
    return global;
}

void ParallelContext::barrier() const {
    validate();
    if (distributed()) {
        MPI_Barrier(communicator);
    }
}

Mesh decompose(
    const Mesh& global,
    const ParallelContext& parallel,
    Index ghost_layers)
{
    parallel.validate();
    global.validate();
    if (global.ownedCellCount() != global.cellCount() || ghost_layers < 1) {
        throw std::invalid_argument("domain decomposition requires a global mesh");
    }
    if (parallel.size == 1) {
        return global;
    }
    const Index global_nx = global.dimensions[0];
    if (global_nx < ghost_layers * parallel.size) {
        throw std::invalid_argument("each MPI partition needs at least ghost_layers owned x cells");
    }

    std::vector<Index> widths(static_cast<std::size_t>(parallel.size));
    Index remaining = global_nx;
    for (int part = 0; part < parallel.size; ++part) {
        widths[static_cast<std::size_t>(part)] =
            remaining / static_cast<Index>(parallel.size - part);
        remaining -= widths[static_cast<std::size_t>(part)];
    }
    if (std::any_of(widths.begin(), widths.end(), [ghost_layers](Index width) {
            return width < ghost_layers;
        })) {
        throw std::invalid_argument("an MPI partition is narrower than its halo");
    }

    Index owned_global_begin = 0;
    for (int part = 0; part < parallel.rank; ++part) {
        owned_global_begin += widths[static_cast<std::size_t>(part)];
    }
    const Index owned_width = widths[static_cast<std::size_t>(parallel.rank)];
    const Index left_ghost = parallel.rank == 0 ? 0 : ghost_layers;
    const Index right_ghost = parallel.rank + 1 == parallel.size ? 0 : ghost_layers;
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
    Mesh local = Mesh::structured(
        {local_nx, global.dimensions[1], global.dimensions[2]},
        std::move(points), temporaryPatches());
    local.setOwnership(
        global.dimensions, global_offset, left_ghost,
        left_ghost + owned_width, ghost_layers);

    std::vector<BoundaryPatch> patches = global.patches;
    for (auto& patch : patches) {
        patch.faces.clear();
    }
    const Index left_processor_patch = left_ghost == 0
        ? invalid_index : static_cast<Index>(patches.size());
    if (left_processor_patch != invalid_index) {
        patches.push_back({"processor_left", PatchKind::Processor, {}});
    }
    const Index right_processor_patch = right_ghost == 0
        ? invalid_index : static_cast<Index>(patches.size());
    if (right_processor_patch != invalid_index) {
        patches.push_back({"processor_right", PatchKind::Processor, {}});
    }

    for (Index face = 0; face < local.faceCount(); ++face) {
        if (!local.boundaryFace(face)) {
            continue;
        }
        const auto f = static_cast<std::size_t>(face);
        const Side side = static_cast<Side>(local.face_patch[f]);
        Index patch = invalid_index;
        if (side == Side::XMin && left_processor_patch != invalid_index) {
            patch = left_processor_patch;
        } else if (side == Side::XMax && right_processor_patch != invalid_index) {
            patch = right_processor_patch;
        } else {
            const Index owner = local.face_owner[f];
            const Index global_cell = local.globalCellId(owner);
            const Index global_i = global_cell % global.dimensions[0];
            const Index global_j =
                (global_cell / global.dimensions[0]) % global.dimensions[1];
            const Index global_k = global_cell /
                (global.dimensions[0] * global.dimensions[1]);
            const Index source_cell = global.cellId(global_i, global_j, global_k);
            const Index source_face = global.cell_faces[static_cast<std::size_t>(source_cell)]
                [static_cast<std::size_t>(side)];
            patch = global.face_patch[static_cast<std::size_t>(source_face)];
        }
        if (patch < 0 || static_cast<std::size_t>(patch) >= patches.size()) {
            throw std::logic_error("decomposed boundary face has no patch");
        }
        local.face_patch[f] = patch;
        patches[static_cast<std::size_t>(patch)].faces.push_back(face);
    }
    local.patches = std::move(patches);
    local.validate();
    return local;
}

HaloExchange::HaloExchange(const Mesh& mesh, ParallelContext parallel)
    : mesh_(&mesh), parallel_(parallel)
{
    parallel_.validate();
    if (!parallel_.distributed()) {
        return;
    }
    if (mesh.ghost_layers < 1 ||
        mesh.owned_i_begin < (parallel_.rank == 0 ? 0 : mesh.ghost_layers) ||
        mesh.dimensions[0] - mesh.owned_i_end <
            (parallel_.rank + 1 == parallel_.size ? 0 : mesh.ghost_layers)) {
        throw std::invalid_argument("mesh halo layout does not match MPI partition");
    }
    left_ = parallel_.rank == 0 ? MPI_PROC_NULL : parallel_.rank - 1;
    right_ = parallel_.rank + 1 == parallel_.size
        ? MPI_PROC_NULL : parallel_.rank + 1;
    const auto appendPlane = [&](std::vector<Index>& indices, Index begin) {
        for (Index k = 0; k < mesh.dimensions[2]; ++k) {
            for (Index j = 0; j < mesh.dimensions[1]; ++j) {
                for (Index layer = 0; layer < mesh.ghost_layers; ++layer) {
                    indices.push_back(mesh.cellId(begin + layer, j, k));
                }
            }
        }
    };
    if (left_ != MPI_PROC_NULL) {
        appendPlane(send_left_, mesh.owned_i_begin);
        appendPlane(receive_left_, mesh.owned_i_begin - mesh.ghost_layers);
    }
    if (right_ != MPI_PROC_NULL) {
        appendPlane(send_right_, mesh.owned_i_end - mesh.ghost_layers);
        appendPlane(receive_right_, mesh.owned_i_end);
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
    if (left_ != MPI_PROC_NULL) {
        appendInterface(send_face_left_, mesh.owned_i_begin, Side::XMin);
        appendInterface(receive_face_left_, mesh.owned_i_begin, Side::XMin);
    }
    if (right_ != MPI_PROC_NULL) {
        appendInterface(send_face_right_, mesh.owned_i_end - 1, Side::XMax);
        appendInterface(receive_face_right_, mesh.owned_i_end - 1, Side::XMax);
    }
}

void HaloExchange::exchange(double* values, std::size_t components) {
    if (!parallel_.distributed()) {
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

    pack(send_left_, send_buffer_left_);
    pack(send_right_, send_buffer_right_);
    receive_buffer_right_.resize(receive_right_.size() * components);
    receive_buffer_left_.resize(receive_left_.size() * components);
    double dummy = 0.0;
    MPI_Sendrecv(
        send_buffer_left_.empty() ? &dummy : send_buffer_left_.data(),
        mpiCount(send_buffer_left_.size()), MPI_DOUBLE, left_, 101,
        receive_buffer_right_.empty() ? &dummy : receive_buffer_right_.data(),
        mpiCount(receive_buffer_right_.size()), MPI_DOUBLE, right_, 101,
        parallel_.communicator, MPI_STATUS_IGNORE);
    MPI_Sendrecv(
        send_buffer_right_.empty() ? &dummy : send_buffer_right_.data(),
        mpiCount(send_buffer_right_.size()), MPI_DOUBLE, right_, 102,
        receive_buffer_left_.empty() ? &dummy : receive_buffer_left_.data(),
        mpiCount(receive_buffer_left_.size()), MPI_DOUBLE, left_, 102,
        parallel_.communicator, MPI_STATUS_IGNORE);
    unpack(receive_right_, receive_buffer_right_);
    unpack(receive_left_, receive_buffer_left_);
}

void HaloExchange::exchangeFaces(double* values, std::size_t components) {
    if (!parallel_.distributed()) {
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
    pack(send_face_left_, send_face_buffer_left_);
    pack(send_face_right_, send_face_buffer_right_);
    receive_face_buffer_right_.resize(receive_face_right_.size() * components);
    receive_face_buffer_left_.resize(receive_face_left_.size() * components);
    double dummy = 0.0;
    MPI_Sendrecv(
        send_face_buffer_left_.empty() ? &dummy : send_face_buffer_left_.data(),
        mpiCount(send_face_buffer_left_.size()), MPI_DOUBLE, left_, 201,
        receive_face_buffer_right_.empty() ? &dummy : receive_face_buffer_right_.data(),
        mpiCount(receive_face_buffer_right_.size()), MPI_DOUBLE, right_, 201,
        parallel_.communicator, MPI_STATUS_IGNORE);
    MPI_Sendrecv(
        send_face_buffer_right_.empty() ? &dummy : send_face_buffer_right_.data(),
        mpiCount(send_face_buffer_right_.size()), MPI_DOUBLE, right_, 202,
        receive_face_buffer_left_.empty() ? &dummy : receive_face_buffer_left_.data(),
        mpiCount(receive_face_buffer_left_.size()), MPI_DOUBLE, left_, 202,
        parallel_.communicator, MPI_STATUS_IGNORE);
    unpack(receive_face_right_, receive_face_buffer_right_);
    unpack(receive_face_left_, receive_face_buffer_left_);
}

void HaloExchange::exchange(std::vector<double>& values) {
    if (values.size() != static_cast<std::size_t>(mesh_->cellCount())) {
        throw std::invalid_argument("raw halo field has the wrong size");
    }
    exchange(values.data(), 1);
}

void HaloExchange::exchange(ScalarField& field) {
    if (&field.mesh() != mesh_) {
        throw std::invalid_argument("scalar halo field is incompatible");
    }
    if (field.location() == FieldLocation::Cell) {
        exchange(field.data(), 1);
    } else if (field.location() == FieldLocation::Face) {
        exchangeFaces(field.data(), 1);
    } else {
        throw std::invalid_argument("vertex scalar halo exchange is not supported");
    }
}

void HaloExchange::exchange(VectorField& field) {
    if (&field.mesh() != mesh_) {
        throw std::invalid_argument("vector halo field is incompatible");
    }
    if (field.location() == FieldLocation::Cell) {
        exchange(&field.data()->x, 3);
    } else if (field.location() == FieldLocation::Face) {
        exchangeFaces(&field.data()->x, 3);
    } else {
        throw std::invalid_argument("vertex vector halo exchange is not supported");
    }
}

void HaloExchange::exchange(TensorField& field) {
    if (&field.mesh() != mesh_) {
        throw std::invalid_argument("tensor halo field is incompatible");
    }
    if (field.location() == FieldLocation::Cell) {
        exchange(&field.data()->rows[0].x, 9);
    } else if (field.location() == FieldLocation::Face) {
        exchangeFaces(&field.data()->rows[0].x, 9);
    } else {
        throw std::invalid_argument("vertex tensor halo exchange is not supported");
    }
}

void writeOwnedCsv(
    const std::filesystem::path& directory,
    const VectorField& velocity,
    const ScalarField& pressure,
    const ParallelContext& parallel,
    const std::string& prefix)
{
    const Mesh& mesh = velocity.mesh();
    if (&pressure.mesh() != &mesh ||
        velocity.location() != FieldLocation::Cell ||
        pressure.location() != FieldLocation::Cell || prefix.empty()) {
        throw std::invalid_argument("parallel output fields are incompatible");
    }
    int directory_error = 0;
    if (parallel.rank == 0) {
        std::error_code error;
        std::filesystem::create_directories(directory, error);
        directory_error = error ? 1 : 0;
    }
    if (parallel.maximum(directory_error) != 0) {
        throw std::runtime_error("cannot create parallel output directory: " +
                                 directory.string());
    }
    parallel.barrier();
    const auto path = directory /
        (prefix + "_" + std::to_string(parallel.rank) + ".csv");
    std::ofstream output(path);
    if (parallel.maximum(output ? 0 : 1) != 0) {
        throw std::runtime_error("cannot create parallel field output: " + path.string());
    }
    output << "global_id,x,y,z,u,v,w,p\n" << std::setprecision(17);
    for (Index cell : mesh.owned_cells) {
        const Vec3& point = mesh.cell_centres[static_cast<std::size_t>(cell)];
        const Vec3& value = velocity[cell];
        output << mesh.globalCellId(cell) << ','
               << point.x << ',' << point.y << ',' << point.z << ','
               << value.x << ',' << value.y << ',' << value.z << ','
               << pressure[cell] << '\n';
    }
    output.close();
    if (parallel.maximum(output ? 0 : 1) != 0) {
        throw std::runtime_error("cannot flush parallel field output: " + path.string());
    }
    parallel.barrier();
}

}  // namespace babelsim
