#pragma once

#include "babelsim/solver_control.h"

#include <Eigen/Sparse>

#include <memory>

namespace babelsim::detail {

// 默认计算后端的轻量聚合 AMG。它只接收代数矩阵，不知道 Mesh、Field 或 MPI；分布式
// 求解器将其作为每个子域的可复用块预条件器，并继续负责跨分区矩阵向量乘和全局归约。
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

    // 一个 V-cycle，可作为 CG/BiCGSTAB/GMRES 的预条件器。输出与输入不允许别名。
    bool apply(const Eigen::VectorXd& input, Eigen::VectorXd& output);
    // 以重复 V-cycle 直接求解；用于显式选择 amg none 的线性系统。
    SolveResult solve(const Eigen::VectorXd& right_hand_side, Eigen::VectorXd& solution);

private:
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
};

}  // babelsim::detail 命名空间
