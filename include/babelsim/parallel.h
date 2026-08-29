#pragma once

#include "babelsim/field.h"

#include <mpi.h>

#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace babelsim {

// A serial default keeps the mathematical core usable without MPI_Init.
// world() is the explicit entry point for MPI applications.
struct ParallelContext {
    MPI_Comm communicator = MPI_COMM_NULL;
    int rank = 0;
    int size = 1;

    static ParallelContext world(MPI_Comm communicator = MPI_COMM_WORLD);
    bool distributed() const { return size > 1; }
    void validate() const;
    void sum(const double* local, double* global, int count) const;
    void maximum(const double* local, double* global, int count) const;
    int sum(int local) const;
    int maximum(int local) const;
    void barrier() const;
};

// Structured x decomposition. Two layers are the default because corrected
// face diffusion may consume a gradient reconstructed in the first ghost cell.
Mesh decompose(
    const Mesh& global,
    const ParallelContext& parallel,
    Index ghost_layers = 2);

template <typename T>
void copyBoundaryConditions(const Field<T>& global, Field<T>& local) {
    if (global.location() != FieldLocation::Cell ||
        local.location() != FieldLocation::Cell ||
        global.mesh().patches.size() > local.mesh().patches.size()) {
        throw std::invalid_argument("boundary-condition copy fields are incompatible");
    }
    for (Index patch = 0;
         patch < static_cast<Index>(global.mesh().patches.size()); ++patch) {
        local.setBoundary(patch, global.boundary(patch));
    }
}

// Packs strided x-normal planes into persistent buffers. Field storage remains
// contiguous and independent of MPI data types.
class HaloExchange {
public:
    HaloExchange(const Mesh& mesh, ParallelContext parallel);

    void exchange(std::vector<double>& values);
    void exchange(ScalarField& field);
    void exchange(VectorField& field);
    void exchange(TensorField& field);

private:
    void exchange(double* values, std::size_t components);
    void exchangeFaces(double* values, std::size_t components);

    const Mesh* mesh_;
    ParallelContext parallel_;
    int left_ = MPI_PROC_NULL;
    int right_ = MPI_PROC_NULL;
    std::vector<Index> send_left_;
    std::vector<Index> send_right_;
    std::vector<Index> receive_left_;
    std::vector<Index> receive_right_;
    std::vector<Index> send_face_left_;
    std::vector<Index> send_face_right_;
    std::vector<Index> receive_face_left_;
    std::vector<Index> receive_face_right_;
    std::vector<double> send_buffer_left_;
    std::vector<double> send_buffer_right_;
    std::vector<double> receive_buffer_left_;
    std::vector<double> receive_buffer_right_;
    std::vector<double> send_face_buffer_left_;
    std::vector<double> send_face_buffer_right_;
    std::vector<double> receive_face_buffer_left_;
    std::vector<double> receive_face_buffer_right_;
};

void writeOwnedCsv(
    const std::filesystem::path& directory,
    const VectorField& velocity,
    const ScalarField& pressure,
    const ParallelContext& parallel,
    const std::string& prefix = "fields");

}  // namespace babelsim
