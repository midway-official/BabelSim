#include "babelsim/linear_solver.h"

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>

namespace babelsim {
namespace {

template <typename Solver>
SolveResult runPreparedSolver(
    Solver& solver,
    const Eigen::SparseMatrix<double>& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& x,
    const LinearSolverConfig& config,
    double initial_residual,
    double target)
{
    solver.setTolerance(target / std::max(b.norm(), 1e-30));
    if (config.warm_start) {
        x = solver.solveWithGuess(b, x);
    } else {
        x = solver.solve(b);
    }
    const double final_residual = (b - A * x).norm();
    const double scale = std::max({initial_residual, b.norm(), 1e-30});
    SolveStatus status = SolveStatus::NumericalFailure;
    if (std::isfinite(final_residual) && final_residual <= target * (1.0 + 1e-10)) {
        status = SolveStatus::Converged;
    } else if (solver.info() == Eigen::NoConvergence && std::isfinite(final_residual)) {
        status = SolveStatus::MaxIterations;
    }
    return {
        status,
        static_cast<int>(solver.iterations()),
        initial_residual,
        final_residual,
        final_residual / scale,
    };
}

}  // 匿名命名空间

struct PreparedLinearSolver::Implementation {
    using ConjugateGradient = Eigen::ConjugateGradient<
        Eigen::SparseMatrix<double>,
        Eigen::Lower | Eigen::Upper,
        Eigen::IncompleteCholesky<double>>;
    using BiCGSTAB = Eigen::BiCGSTAB<
        Eigen::SparseMatrix<double>,
        Eigen::IncompleteLUT<double>>;

    explicit Implementation(LinearSolverConfig value)
        : config(std::move(value))
    {}

    LinearSolverConfig config;
    const Eigen::SparseMatrix<double>* matrix = nullptr;
    ConjugateGradient conjugate_gradient;
    BiCGSTAB bicgstab;
    bool pattern_analyzed = false;
    bool factorization_succeeded = false;
};

void LinearSolverConfig::validate() const {
    if (!(absolute_tolerance > 0.0) || !(relative_tolerance > 0.0) ||
        !std::isfinite(absolute_tolerance) ||
        !std::isfinite(relative_tolerance) || max_iterations <= 0) {
        throw std::invalid_argument("linear solver tolerances or iteration limit are invalid");
    }
    const bool supported =
        (solver == LinearSolverType::ConjugateGradient &&
         preconditioner == PreconditionerType::IncompleteCholesky) ||
        (solver == LinearSolverType::BiCGSTAB &&
         preconditioner == PreconditionerType::ILUT);
    if (!supported) {
        throw std::invalid_argument("unsupported linear solver/preconditioner pair");
    }
}

PreparedLinearSolver::PreparedLinearSolver(LinearSolverConfig config)
    : implementation_(std::make_unique<Implementation>(std::move(config)))
{
    implementation_->config.validate();
}

PreparedLinearSolver::~PreparedLinearSolver() = default;
PreparedLinearSolver::PreparedLinearSolver(PreparedLinearSolver&&) noexcept = default;
PreparedLinearSolver& PreparedLinearSolver::operator=(
    PreparedLinearSolver&&) noexcept = default;

void PreparedLinearSolver::compute(const Eigen::SparseMatrix<double>& A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("linear-system matrix must be square");
    }
    auto& state = *implementation_;
    state.matrix = &A;
    state.pattern_analyzed = false;
    state.factorization_succeeded = false;
    if (state.config.solver == LinearSolverType::ConjugateGradient) {
        state.conjugate_gradient.setMaxIterations(state.config.max_iterations);
        state.conjugate_gradient.compute(A);
        state.pattern_analyzed = true;
        state.factorization_succeeded =
            state.conjugate_gradient.info() == Eigen::Success;
        return;
    }
    state.bicgstab.setMaxIterations(state.config.max_iterations);
    state.bicgstab.preconditioner().setDroptol(1e-3);
    state.bicgstab.preconditioner().setFillfactor(2);
    // 对 Eigen 3.4 的 IncompleteLUT 直接调用 IterativeSolverBase::analyzePattern
    // 不安全，因为该预条件器直到 factorize 才初始化 ComputationInfo。compute() 会
    // 安全地完成两步。
    state.bicgstab.compute(A);
    state.pattern_analyzed = true;
    state.factorization_succeeded = state.bicgstab.info() == Eigen::Success;
}

void PreparedLinearSolver::factorize(
    const Eigen::SparseMatrix<double>& A)
{
    auto& state = *implementation_;
    if (!state.pattern_analyzed || A.rows() != A.cols()) {
        throw std::logic_error(
            "linear-solver pattern must be analyzed before factorization");
    }
    state.matrix = &A;
    state.factorization_succeeded = false;
    if (state.config.solver == LinearSolverType::ConjugateGradient) {
        state.conjugate_gradient.factorize(A);
        state.factorization_succeeded =
            state.conjugate_gradient.info() == Eigen::Success;
        return;
    }
    state.bicgstab.factorize(A);
    state.factorization_succeeded = state.bicgstab.info() == Eigen::Success;
}

SolveResult PreparedLinearSolver::solve(
    const Eigen::VectorXd& b,
    Eigen::VectorXd& x)
{
    auto& state = *implementation_;
    if (state.matrix == nullptr || state.matrix->rows() != b.size()) {
        throw std::invalid_argument("linear system dimensions are inconsistent");
    }
    if (!state.config.warm_start || x.size() != b.size()) {
        x = Eigen::VectorXd::Zero(b.size());
    }

    const auto& A = *state.matrix;
    const double initial_residual = (b - A * x).norm();
    const double scale = std::max({initial_residual, b.norm(), 1e-30});
    const double target = std::max(
        state.config.absolute_tolerance,
        state.config.relative_tolerance * scale);
    if (initial_residual <= target) {
        return {
            SolveStatus::Converged,
            0,
            initial_residual,
            initial_residual,
            initial_residual / scale,
        };
    }
    if (!state.factorization_succeeded) {
        return {
            SolveStatus::NumericalFailure,
            0,
            initial_residual,
            initial_residual,
            initial_residual / scale,
        };
    }

    if (state.config.solver == LinearSolverType::ConjugateGradient) {
        return runPreparedSolver(
            state.conjugate_gradient, A, b, x, state.config,
            initial_residual, target);
    }
    return runPreparedSolver(
        state.bicgstab, A, b, x, state.config, initial_residual, target);
}

SolveResult solve(
    const Eigen::SparseMatrix<double>& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& x,
    const LinearSolverConfig& config)
{
    PreparedLinearSolver solver(config);
    solver.compute(A);
    return solver.solve(b, x);
}

}  // babelsim 命名空间
