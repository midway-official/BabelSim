#pragma once

#include "babelsim/field.h"

#include <mpi.h>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace babelsim {

// 串行默认上下文使数学核心无需 MPI_Init 也可使用；world() 是 MPI 应用的显式入口。
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

// 结构化 x 向分区。默认两层 ghost，因为修正面扩散可能读取第一层 ghost cell 中重构的梯度。
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

// 将非连续的 x 法向平面打包到持久缓冲区。Field 存储保持连续，且不依赖 MPI 数据类型。
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

}  // babelsim 命名空间
