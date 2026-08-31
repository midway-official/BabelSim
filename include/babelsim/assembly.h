#pragma once

#include "babelsim/discrete_equation.h"

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

    void update(const ScalarDiscreteEquation& equation);
    void update(const VectorDiscreteEquation& equation);
    const Eigen::SparseMatrix<double>& matrix() const { return m_matrix; }

private:
    void update(
        const Mesh* equation_mesh,
        const std::vector<double>& diagonal,
        const std::vector<double>& upper,
        const std::vector<double>& lower);

    const Mesh* m_mesh;
    Eigen::SparseMatrix<double> m_matrix;
    std::vector<Eigen::Index> m_diagonal_positions;
    std::vector<Eigen::Index> m_upper_positions;
    std::vector<Eigen::Index> m_lower_positions;
    std::vector<Index> m_coupled_faces;
};

Eigen::SparseMatrix<double> assembleMatrix(const ScalarDiscreteEquation& equation);
Eigen::SparseMatrix<double> assembleMatrix(const VectorDiscreteEquation& equation);
LinearSystem assemble(const ScalarDiscreteEquation& equation);
void assembleSource(const ScalarDiscreteEquation& equation, Eigen::VectorXd& result);
void assembleSource(
    const VectorDiscreteEquation& equation,
    std::array<Eigen::VectorXd, 3>& result);
std::array<Eigen::VectorXd, 3> assembleSource(const VectorDiscreteEquation& equation);

}  // babelsim 命名空间
