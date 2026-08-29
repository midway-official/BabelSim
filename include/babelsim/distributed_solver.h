#pragma once

#include "babelsim/linear_solver.h"
#include "babelsim/parallel.h"

#include <memory>

namespace babelsim {

// Distributed Krylov solver with a rank-local ILUT/IC block preconditioner.
// Cross-rank LDU coefficients are applied after a halo exchange, while every
// dot product and residual norm is reduced over the supplied communicator.
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

}  // namespace babelsim
