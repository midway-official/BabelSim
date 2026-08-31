#pragma once

#include "babelsim/solver_control.h"

#include <Eigen/Sparse>

#include <memory>

namespace babelsim {

// 将矩阵准备与求解分离，使多个右端项复用排序和预条件器工作，同时不向物理层暴露 Eigen。
class PreparedLinearSolver {
public:
    explicit PreparedLinearSolver(LinearSolverConfig config = {});
    ~PreparedLinearSolver();
    PreparedLinearSolver(PreparedLinearSolver&&) noexcept;
    PreparedLinearSolver& operator=(PreparedLinearSolver&&) noexcept;
    PreparedLinearSolver(const PreparedLinearSolver&) = delete;
    PreparedLinearSolver& operator=(const PreparedLinearSolver&) = delete;

    // 复用 compute() 的稀疏结构分析；A 必须具有相同维度和非零模式。
    void factorize(const Eigen::SparseMatrix<double>& A);
    void compute(const Eigen::SparseMatrix<double>& A);
    SolveResult solve(const Eigen::VectorXd& b, Eigen::VectorXd& x);

private:
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
};

SolveResult solve(
    const Eigen::SparseMatrix<double>& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& x,
    const LinearSolverConfig& config = {});

}  // babelsim 命名空间
