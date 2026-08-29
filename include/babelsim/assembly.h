#pragma once

#include "babelsim/equation.h"

#include <Eigen/Sparse>

#include <array>
#include <vector>

namespace babelsim {

struct LinearSystem {
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
};

// 预计算依赖网格的压缩稀疏结构及每个 LDU 系数在值数组中的直接位置。分区网格仅
// 生成自有行；跨进程的面系数由分布式求解器的 halo 矩阵向量乘施加。因此 update()
// 不遍历拓扑，不分配三元组，不排序，也不插入稀疏项。
class SparseAssembly {
public:
    explicit SparseAssembly(const Mesh& mesh);
    SparseAssembly(const Mesh&&) = delete;

    void update(const ScalarEquation& equation);
    void update(const VectorEquation& equation);
    const Eigen::SparseMatrix<double>& matrix() const { return matrix_; }

private:
    void update(
        const Mesh* equation_mesh,
        const std::vector<double>& diagonal,
        const std::vector<double>& upper,
        const std::vector<double>& lower);

    const Mesh* mesh_;
    Eigen::SparseMatrix<double> matrix_;
    std::vector<Eigen::Index> diagonal_positions_;
    std::vector<Eigen::Index> upper_positions_;
    std::vector<Eigen::Index> lower_positions_;
};

Eigen::SparseMatrix<double> assembleMatrix(const ScalarEquation& equation);
Eigen::SparseMatrix<double> assembleMatrix(const VectorEquation& equation);
LinearSystem assemble(const ScalarEquation& equation);
void assembleSource(const ScalarEquation& equation, Eigen::VectorXd& result);
void assembleSource(
    const VectorEquation& equation,
    std::array<Eigen::VectorXd, 3>& result);
std::array<Eigen::VectorXd, 3> assembleSource(const VectorEquation& equation);

}  // babelsim 命名空间
