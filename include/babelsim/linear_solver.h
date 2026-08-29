#pragma once

#include "babelsim/assembly.h"

#include <memory>

namespace babelsim {

enum class LinearSolverType {
    ConjugateGradient,
    BiCGSTAB,
};

enum class PreconditionerType {
    IncompleteCholesky,
    ILUT,
};

enum class SolveStatus {
    Converged,
    MaxIterations,
    NumericalFailure,
};

struct LinearSolverConfig {
    LinearSolverType solver = LinearSolverType::BiCGSTAB;
    PreconditionerType preconditioner = PreconditionerType::ILUT;
    double absolute_tolerance = 1e-12;
    double relative_tolerance = 1e-8;
    int max_iterations = 1000;
    bool warm_start = false;

    void validate() const;
};

struct SolveResult {
    SolveStatus status = SolveStatus::NumericalFailure;
    int iterations = 0;
    double initial_residual = 0.0;
    double final_residual = 0.0;
    double relative_residual = 0.0;

    bool converged() const { return status == SolveStatus::Converged; }
};

// Separating matrix setup from solve permits several right-hand sides to
// share ordering and preconditioner work without exposing Eigen to physics.
class PreparedLinearSolver {
public:
    explicit PreparedLinearSolver(LinearSolverConfig config = {});
    ~PreparedLinearSolver();
    PreparedLinearSolver(PreparedLinearSolver&&) noexcept;
    PreparedLinearSolver& operator=(PreparedLinearSolver&&) noexcept;
    PreparedLinearSolver(const PreparedLinearSolver&) = delete;
    PreparedLinearSolver& operator=(const PreparedLinearSolver&) = delete;

    // Reuses the sparsity analysis from compute(); A must have the same
    // dimensions and nonzero pattern.
    void factorize(const Eigen::SparseMatrix<double>& A);
    void compute(const Eigen::SparseMatrix<double>& A);
    SolveResult solve(const Eigen::VectorXd& b, Eigen::VectorXd& x);

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

SolveResult solve(
    const Eigen::SparseMatrix<double>& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& x,
    const LinearSolverConfig& config = {});

}  // namespace babelsim
