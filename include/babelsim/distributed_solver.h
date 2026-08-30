#pragma once

#include "babelsim/linear_solver.h"
#include "babelsim/parallel.h"

#include <memory>

namespace babelsim {

// 带各 rank 局部 ILUT/IC 块预条件器的分布式 Krylov 求解器。跨 rank 的 LDU 系数在
// halo 交换后施加；所有点积和残差范数均在给定通信器上归约。
class DistributedLinearSolver {
public:
    DistributedLinearSolver(
        const Mesh& mesh,
        ParallelContext parallel,
        LinearSolverConfig config = {});
    DistributedLinearSolver(
        const Mesh&&,
        ParallelContext,
        LinearSolverConfig = {}) = delete;
    ~DistributedLinearSolver();
    DistributedLinearSolver(DistributedLinearSolver&&) noexcept;
    DistributedLinearSolver& operator=(DistributedLinearSolver&&) noexcept;
    DistributedLinearSolver(const DistributedLinearSolver&) = delete;
    DistributedLinearSolver& operator=(const DistributedLinearSolver&) = delete;

    void compute(
        const Eigen::SparseMatrix<double>& local_matrix,
        const ScalarEquation& equation);
    void compute(
        const Eigen::SparseMatrix<double>& local_matrix,
        const VectorEquation& equation);
    void factorize(
        const Eigen::SparseMatrix<double>& local_matrix,
        const ScalarEquation& equation);
    void factorize(
        const Eigen::SparseMatrix<double>& local_matrix,
        const VectorEquation& equation);
    SolveResult solve(const Eigen::VectorXd& b, Eigen::VectorXd& x);

private:
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
};

}  // babelsim 命名空间
