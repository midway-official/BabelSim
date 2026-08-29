#pragma once

#include "babelsim/linear_solver.h"
#include "babelsim/parallel.h"

#include <memory>

namespace babelsim {

// 带 rank-local ILUT/IC 块预条件器的分布式 Krylov 求解器。跨 rank LDU 系数在
// halo exchange 后施加；所有点积和残差范数均在给定 communicator 上归约。
class DistributedLinearSolver {
public:
    DistributedLinearSolver(
        const Mesh& mesh,
        ParallelContext parallel,
        LinearSolverConfig config = {});
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
    std::unique_ptr<Implementation> implementation_;
};

}  // babelsim 命名空间
