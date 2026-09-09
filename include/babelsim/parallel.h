#pragma once

#include "babelsim/field.h"

#include <mpi.h>

#include <cstddef>
#include <filesystem>
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
    void sum(const int* local, int* global, int count) const;
    void maximum(const double* local, double* global, int count) const;
    int sum(int local) const;
    int maximum(int local) const;
    void barrier() const;
};

// 基于单元邻接图的分区。默认两层 ghost，因为修正面扩散可能读取第一层 ghost
// cell 中重构的梯度。
Mesh decompose(
    const Mesh& global,
    const ParallelContext& parallel,
    Index ghost_layers = 2);

// 并行读取原生网格：rank 0 负责磁盘读取，分区细节完全留在 Parallel 层；
// 返回值始终是当前 rank 的局部 Mesh。串行时退化为 readMeshFile。
Mesh readDistributedMesh(
    const std::filesystem::path& path,
    const ParallelContext& parallel,
    Index ghost_layers = 2);

template <typename T>
void copyBoundaryConditions(const Field<T>& global, Field<T>& local) {
    if (global.location() != FieldLocation::Cell ||
        local.location() != FieldLocation::Cell ||
        global.mesh().patchCount() > local.mesh().patchCount()) {
        throw std::invalid_argument("boundary-condition copy fields are incompatible");
    }
    for (Index patch = 0;
         patch < global.mesh().patchCount(); ++patch) {
        local.setBoundary(patch, global.boundary(patch));
    }
}

// 将任意拓扑的 cell/face 邻居打包到持久缓冲区。Field 存储保持连续，且不依赖
// MPI 数据类型。值的发布者由对应实体的 owner rank 唯一确定。
class HaloExchange {
public:
    HaloExchange(const Mesh& mesh, ParallelContext parallel);

    void exchange(std::vector<double>& values);
    // 分布式稀疏矩阵乘只访问接口两侧第一层 ghost；该入口避免传输非正交
    // 重构所需的第二层 ghost，从而减少 Krylov 热路径的通信量。
    void exchangeFirstLayer(std::vector<double>& values);
    void exchange(ScalarField& field);
    void exchange(VectorField& field);
    void exchange(TensorField& field);

private:
    struct ExchangePlan {
        std::vector<int> send_counts;
        std::vector<int> send_offsets;
        std::vector<int> receive_counts;
        std::vector<int> receive_offsets;
        std::vector<Index> send_indices;
        std::vector<Index> receive_indices;
        std::vector<double> send_buffer;
        std::vector<double> receive_buffer;
    };

    void exchange(double* values, std::size_t components);
    void exchange(double* values, std::size_t components, ExchangePlan& plan);
    void exchangeFaces(double* values, std::size_t components);

    const Mesh* m_mesh;
    ParallelContext m_parallel;
    ExchangePlan m_cells;
    ExchangePlan m_first_layer_cells;
    ExchangePlan m_faces;
};

}  // babelsim 命名空间
