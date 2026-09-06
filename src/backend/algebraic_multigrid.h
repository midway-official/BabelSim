#pragma once

#include "babelsim/solver_control.h"

#include <Eigen/Sparse>

#include <cstdint>
#include <memory>

namespace babelsim::detail {

// 默认计算后端的串行聚合 AMG。它只接收代数矩阵，不知道 Mesh、Field 或 MPI；
// 分布式求解器在代数层实现全局粗空间，不复用这个串行层级。
class AlgebraicMultigrid {
public:
    explicit AlgebraicMultigrid(LinearSolverConfig config);
    ~AlgebraicMultigrid();
    AlgebraicMultigrid(AlgebraicMultigrid&&) noexcept;
    AlgebraicMultigrid& operator=(AlgebraicMultigrid&&) noexcept;
    AlgebraicMultigrid(const AlgebraicMultigrid&) = delete;
    AlgebraicMultigrid& operator=(const AlgebraicMultigrid&) = delete;

    // compute 建立固定聚合和粗网格模式；factorize 复用该模式，只更新系数与末层分解。
    void compute(const Eigen::SparseMatrix<double>& matrix);
    void factorize(const Eigen::SparseMatrix<double>& matrix);
    bool ready() const;

    // 一个 V-cycle，只作为 CG/BiCGSTAB 的预条件器。输出与输入不允许别名。
    bool apply(const Eigen::VectorXd& input, Eigen::VectorXd& output);
    std::uint64_t lastSparseMatvecs() const;
    double lastSparseMatvecSeconds() const;

private:
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
};

}  // babelsim::detail 命名空间
