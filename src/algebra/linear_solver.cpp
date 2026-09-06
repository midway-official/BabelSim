#include "babelsim/linear_solver.h"

#include "backend/algebraic_multigrid.h"

#include <Eigen/IterativeLinearSolvers>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <limits>
#include <stdexcept>
#include <utility>

namespace babelsim {
namespace {

constexpr double breakdown_tolerance = 1e-30;
using Clock = std::chrono::steady_clock;

double secondsSince(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

bool usesAmg(const LinearSolverConfig& config) {
    return config.preconditioner == PreconditionerType::AlgebraicMultigrid;
}

struct KrylovWorkspace {
    Eigen::VectorXd residual;
    Eigen::VectorXd product;
    Eigen::VectorXd shadow;
    Eigen::VectorXd direction;
    Eigen::VectorXd preconditioned_direction;
    Eigen::VectorXd direction_product;
    Eigen::VectorXd intermediate;
    Eigen::VectorXd preconditioned_intermediate;
    Eigen::VectorXd intermediate_product;

    void resize(Eigen::Index rows) {
        residual.resize(rows);
        product.resize(rows);
        shadow.resize(rows);
        direction.resize(rows);
        preconditioned_direction.resize(rows);
        direction_product.resize(rows);
        intermediate.resize(rows);
        preconditioned_intermediate.resize(rows);
        intermediate_product.resize(rows);
    }
};

}  // 匿名命名空间

struct PreparedLinearSolver::Implementation {
    explicit Implementation(LinearSolverConfig value) : config(std::move(value)) {}

    void apply(const Eigen::VectorXd& input, Eigen::VectorXd& output) {
        const Clock::time_point start = Clock::now();
        output.noalias() = matrix * input;
        ++current_performance.sparse_matvecs;
        current_performance.sparse_matvec_seconds += secondsSince(start);
    }

    bool precondition(const Eigen::VectorXd& input, Eigen::VectorXd& output) {
        const Clock::time_point start = Clock::now();
        bool success = false;
        if (usesAmg(config)) {
            success = amg && amg->apply(input, output);
            if (amg) {
                current_performance.sparse_matvecs += amg->lastSparseMatvecs();
                current_performance.sparse_matvec_seconds +=
                    amg->lastSparseMatvecSeconds();
            }
        } else if (config.solver == LinearSolverType::ConjugateGradient) {
            output = incomplete_cholesky.solve(input);
            success = incomplete_cholesky.info() == Eigen::Success;
        } else {
            output = ilut.solve(input);
            success = ilut.info() == Eigen::Success;
        }
        ++current_performance.preconditioner_applications;
        current_performance.preconditioner_apply_seconds += secondsSince(start);
        return success && output.allFinite();
    }

    SolveResult finish(
        SolveStatus status,
        int iterations,
        double initial_residual,
        double scale,
        const Eigen::VectorXd& right_hand_side,
        const Eigen::VectorXd& solution)
    {
        apply(solution, workspace.product);
        const double final_residual = (right_hand_side - workspace.product).norm();
        const double target = std::max(
            config.absolute_tolerance, config.relative_tolerance * scale);
        // 当 b-Ax 已达到稀疏乘与向量范数的舍入误差下限时，继续迭代只会
        // 触发 Krylov breakdown。该阈值不改变正常尺度问题的用户容差。
        const double roundoff = 64.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, right_hand_side.norm() + matrix.norm() * solution.norm());
        const double attainable_target = std::max(target, roundoff);
        if (!std::isfinite(final_residual)) status = SolveStatus::NumericalFailure;
        else if (final_residual <= attainable_target * (1.0 + 1e-10))
            status = SolveStatus::Converged;
        else if (status == SolveStatus::Converged) status = SolveStatus::MaxIterations;
        SolveResult result{status, iterations, initial_residual, final_residual,
                           final_residual / scale};
        current_performance.linear_solves = 1;
        current_performance.krylov_iterations = static_cast<std::uint64_t>(iterations);
        result.performance = current_performance;
        return result;
    }

    SolveResult solvePcg(
        const Eigen::VectorXd& right_hand_side,
        Eigen::VectorXd& solution,
        double initial_residual,
        double target,
        double scale)
    {
        if (!precondition(workspace.residual, workspace.preconditioned_direction)) {
            return finish(SolveStatus::NumericalFailure, 0, initial_residual, scale,
                          right_hand_side, solution);
        }
        double rho = workspace.residual.dot(workspace.preconditioned_direction);
        workspace.direction = workspace.preconditioned_direction;
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        for (int iteration = 1; iteration <= config.max_iterations; ++iteration) {
            apply(workspace.direction, workspace.direction_product);
            const double denominator = workspace.direction.dot(workspace.direction_product);
            if (!std::isfinite(denominator) || denominator <= breakdown_tolerance ||
                !std::isfinite(rho) || std::abs(rho) <= breakdown_tolerance) break;
            const double alpha = rho / denominator;
            solution.noalias() += alpha * workspace.direction;
            workspace.residual.noalias() -= alpha * workspace.direction_product;
            iterations = iteration;
            if (workspace.residual.norm() <= target) {
                status = SolveStatus::Converged;
                break;
            }
            if (!precondition(workspace.residual, workspace.preconditioned_direction)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            const double next = workspace.residual.dot(workspace.preconditioned_direction);
            if (!std::isfinite(next)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            workspace.direction = workspace.preconditioned_direction +
                (next / rho) * workspace.direction;
            rho = next;
        }
        return finish(status, iterations, initial_residual, scale, right_hand_side, solution);
    }

    SolveResult solveBicgstab(
        const Eigen::VectorXd& right_hand_side,
        Eigen::VectorXd& solution,
        double initial_residual,
        double target,
        double scale)
    {
        workspace.shadow = workspace.residual;
        workspace.direction.setZero();
        workspace.direction_product.setZero();
        double previous_rho = 1.0;
        double alpha = 1.0;
        double omega = 1.0;
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        for (int iteration = 1; iteration <= config.max_iterations; ++iteration) {
            const double rho = workspace.shadow.dot(workspace.residual);
            if (!std::isfinite(rho) || !std::isfinite(omega) ||
                std::abs(rho) <= breakdown_tolerance ||
                std::abs(omega) <= breakdown_tolerance) break;
            const double beta = (rho / previous_rho) * (alpha / omega);
            workspace.direction = workspace.residual + beta *
                (workspace.direction - omega * workspace.direction_product);
            if (!precondition(workspace.direction, workspace.preconditioned_direction)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            apply(workspace.preconditioned_direction, workspace.direction_product);
            const double denominator = workspace.shadow.dot(workspace.direction_product);
            if (!std::isfinite(denominator) ||
                std::abs(denominator) <= breakdown_tolerance) break;
            alpha = rho / denominator;
            workspace.intermediate = workspace.residual - alpha * workspace.direction_product;
            iterations = iteration;
            if (workspace.intermediate.norm() <= target) {
                solution.noalias() += alpha * workspace.preconditioned_direction;
                workspace.residual = workspace.intermediate;
                status = SolveStatus::Converged;
                break;
            }
            if (!precondition(workspace.intermediate, workspace.preconditioned_intermediate)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            apply(workspace.preconditioned_intermediate,
                  workspace.intermediate_product);
            const double product_norm = workspace.intermediate_product.squaredNorm();
            if (!std::isfinite(product_norm) || product_norm <= breakdown_tolerance) break;
            omega = workspace.intermediate_product.dot(workspace.intermediate) / product_norm;
            if (!std::isfinite(omega) || std::abs(omega) <= breakdown_tolerance) break;
            solution.noalias() += alpha * workspace.preconditioned_direction +
                omega * workspace.preconditioned_intermediate;
            workspace.residual = workspace.intermediate - omega * workspace.intermediate_product;
            if (workspace.residual.norm() <= target) {
                status = SolveStatus::Converged;
                break;
            }
            previous_rho = rho;
        }
        return finish(status, iterations, initial_residual, scale, right_hand_side, solution);
    }

    LinearSolverConfig config;
    Eigen::SparseMatrix<double> matrix;
    Eigen::IncompleteCholesky<double> incomplete_cholesky;
    Eigen::IncompleteLUT<double> ilut;
    std::unique_ptr<detail::AlgebraicMultigrid> amg;
    KrylovWorkspace workspace;
    PerformanceCounters current_performance;
    bool pattern_analyzed = false;
    bool factorization_succeeded = false;
    int amg_updates_since_factorization = 0;
};

void LinearSolverConfig::validate() const {
    if (!(absolute_tolerance > 0.0) || !(relative_tolerance > 0.0) ||
        !std::isfinite(absolute_tolerance) || !std::isfinite(relative_tolerance) ||
        max_iterations <= 0 || amg_max_levels <= 0 || amg_coarse_size <= 0 ||
        amg_smoothing_steps <= 0 || amg_refresh_interval <= 0) {
        throw std::invalid_argument("linear solver configuration is invalid");
    }
    const bool supported =
        (solver == LinearSolverType::ConjugateGradient &&
         (preconditioner == PreconditionerType::IncompleteCholesky ||
          preconditioner == PreconditionerType::AlgebraicMultigrid)) ||
        (solver == LinearSolverType::BiCGSTAB &&
         (preconditioner == PreconditionerType::ILUT ||
          preconditioner == PreconditionerType::AlgebraicMultigrid));
    if (!supported) {
        throw std::invalid_argument("unsupported linear solver/preconditioner pair");
    }
}

PreparedLinearSolver::PreparedLinearSolver(LinearSolverConfig config)
    : m_implementation(std::make_unique<Implementation>(std::move(config)))
{
    m_implementation->config.validate();
}

PreparedLinearSolver::~PreparedLinearSolver() = default;
PreparedLinearSolver::PreparedLinearSolver(PreparedLinearSolver&&) noexcept = default;
PreparedLinearSolver& PreparedLinearSolver::operator=(PreparedLinearSolver&&) noexcept = default;

void PreparedLinearSolver::compute(const Eigen::SparseMatrix<double>& matrix) {
    if (matrix.rows() != matrix.cols()) {
        throw std::invalid_argument("linear-system matrix must be square");
    }
    if (!m_implementation) throw std::logic_error("linear solver is moved-from");
    Implementation& state = *m_implementation;
    state.matrix = matrix;
    state.workspace.resize(state.matrix.rows());
    state.pattern_analyzed = false;
    state.factorization_succeeded = false;
    if (usesAmg(state.config)) {
        state.amg = std::make_unique<detail::AlgebraicMultigrid>(state.config);
        state.amg->compute(state.matrix);
        state.pattern_analyzed = true;
        state.factorization_succeeded = state.amg->ready();
        state.amg_updates_since_factorization = 0;
        return;
    }
    if (state.config.solver == LinearSolverType::ConjugateGradient) {
        state.incomplete_cholesky.compute(state.matrix);
        state.pattern_analyzed = true;
        state.factorization_succeeded =
            state.incomplete_cholesky.info() == Eigen::Success;
        return;
    }
    state.ilut.setDroptol(1e-3);
    state.ilut.setFillfactor(2);
    state.ilut.compute(state.matrix);
    state.pattern_analyzed = true;
    state.factorization_succeeded = state.ilut.info() == Eigen::Success;
}

void PreparedLinearSolver::factorize(const Eigen::SparseMatrix<double>& matrix) {
    if (!m_implementation) throw std::logic_error("linear solver is moved-from");
    Implementation& state = *m_implementation;
    if (!state.pattern_analyzed || matrix.rows() != matrix.cols()) {
        throw std::logic_error("linear-solver pattern must be analyzed before factorization");
    }
    state.matrix = matrix;
    state.factorization_succeeded = false;
    if (usesAmg(state.config)) {
        if (++state.amg_updates_since_factorization >= state.config.amg_refresh_interval) {
            state.amg->factorize(state.matrix);
            state.amg_updates_since_factorization = 0;
        }
        state.factorization_succeeded = state.amg->ready();
        return;
    }
    if (state.config.solver == LinearSolverType::ConjugateGradient) {
        state.incomplete_cholesky.factorize(state.matrix);
        state.factorization_succeeded =
            state.incomplete_cholesky.info() == Eigen::Success;
        return;
    }
    state.ilut.factorize(state.matrix);
    state.factorization_succeeded = state.ilut.info() == Eigen::Success;
}

SolveResult PreparedLinearSolver::solve(
    const Eigen::VectorXd& right_hand_side,
    Eigen::VectorXd& solution)
{
    if (!m_implementation) throw std::logic_error("linear solver is moved-from");
    Implementation& state = *m_implementation;
    if (!state.pattern_analyzed || state.matrix.rows() != right_hand_side.size()) {
        throw std::invalid_argument("linear system dimensions are inconsistent");
    }
    if (!state.config.warm_start || solution.size() != right_hand_side.size()) {
        solution.setZero(right_hand_side.size());
    }
    state.current_performance = {};
    state.apply(solution, state.workspace.product);
    state.workspace.residual = right_hand_side - state.workspace.product;
    const double initial_residual = state.workspace.residual.norm();
    const double scale = std::max({initial_residual, right_hand_side.norm(), 1e-30});
    const double target = std::max(
        state.config.absolute_tolerance, state.config.relative_tolerance * scale);
    if (initial_residual <= target || !state.factorization_succeeded) {
        const SolveStatus status = initial_residual <= target
            ? SolveStatus::Converged : SolveStatus::NumericalFailure;
        SolveResult result{status, 0, initial_residual, initial_residual,
                           initial_residual / scale};
        state.current_performance.linear_solves = 1;
        result.performance = state.current_performance;
        return result;
    }
    return state.config.solver == LinearSolverType::ConjugateGradient
        ? state.solvePcg(right_hand_side, solution, initial_residual, target, scale)
        : state.solveBicgstab(right_hand_side, solution, initial_residual, target, scale);
}

SolveResult solve(
    const Eigen::SparseMatrix<double>& matrix,
    const Eigen::VectorXd& right_hand_side,
    Eigen::VectorXd& solution,
    const LinearSolverConfig& config)
{
    PreparedLinearSolver solver(config);
    solver.compute(matrix);
    return solver.solve(right_hand_side, solution);
}

}  // babelsim 命名空间
