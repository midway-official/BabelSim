#include "babelsim/linear_solver.h"

#include "backend/algebraic_multigrid.h"

#include <Eigen/IterativeLinearSolvers>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace babelsim {
namespace {

constexpr double breakdown_tolerance = 1e-30;

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
    if (config.warm_start) x = solver.solveWithGuess(b, x);
    else x = solver.solve(b);
    const double final_residual = (b - A * x).norm();
    const double scale = std::max({initial_residual, b.norm(), 1e-30});
    SolveStatus status = SolveStatus::NumericalFailure;
    if (std::isfinite(final_residual) && final_residual <= target * (1.0 + 1e-10)) {
        status = SolveStatus::Converged;
    } else if (solver.info() == Eigen::NoConvergence && std::isfinite(final_residual)) {
        status = SolveStatus::MaxIterations;
    }
    return {status, static_cast<int>(solver.iterations()), initial_residual,
            final_residual, final_residual / scale};
}

bool usesAmg(const LinearSolverConfig& config) {
    return config.solver == LinearSolverType::AlgebraicMultigrid ||
        config.preconditioner == PreconditionerType::AlgebraicMultigrid;
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
    std::vector<Eigen::VectorXd> basis;
    std::vector<Eigen::VectorXd> preconditioned_basis;
    Eigen::MatrixXd hessenberg;
    Eigen::VectorXd givens_cosine;
    Eigen::VectorXd givens_sine;
    Eigen::VectorXd least_squares;
    Eigen::VectorXd coefficients;

    void resize(Eigen::Index rows, int restart) {
        residual.resize(rows);
        product.resize(rows);
        shadow.resize(rows);
        direction.resize(rows);
        preconditioned_direction.resize(rows);
        direction_product.resize(rows);
        intermediate.resize(rows);
        preconditioned_intermediate.resize(rows);
        intermediate_product.resize(rows);
        basis.resize(static_cast<std::size_t>(restart + 1));
        preconditioned_basis.resize(static_cast<std::size_t>(restart));
        for (Eigen::VectorXd& value : basis) value.resize(rows);
        for (Eigen::VectorXd& value : preconditioned_basis) value.resize(rows);
        hessenberg.resize(restart + 1, restart);
        givens_cosine.resize(restart);
        givens_sine.resize(restart);
        least_squares.resize(restart + 1);
        coefficients.resize(restart);
    }
};

}  // 匿名命名空间

struct PreparedLinearSolver::Implementation {
    using ConjugateGradient = Eigen::ConjugateGradient<
        Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
        Eigen::IncompleteCholesky<double>>;
    using BiCGSTAB = Eigen::BiCGSTAB<
        Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>;

    explicit Implementation(LinearSolverConfig value)
        : config(std::move(value))
    {}

    bool precondition(const Eigen::VectorXd& input, Eigen::VectorXd& output) {
        return amg && amg->apply(input, output) && output.allFinite();
    }

    SolveResult finish(
        SolveStatus status,
        int iterations,
        double initial_residual,
        double scale,
        const Eigen::VectorXd& right_hand_side,
        const Eigen::VectorXd& solution)
    {
        workspace.product.noalias() = matrix * solution;
        const double final_residual = (right_hand_side - workspace.product).norm();
        const double target = std::max(
            config.absolute_tolerance, config.relative_tolerance * scale);
        if (!std::isfinite(final_residual)) status = SolveStatus::NumericalFailure;
        else if (final_residual <= target * (1.0 + 1e-10)) status = SolveStatus::Converged;
        else if (status == SolveStatus::Converged) status = SolveStatus::MaxIterations;
        return {status, iterations, initial_residual, final_residual,
                final_residual / scale};
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
        if (!std::isfinite(rho)) {
            return finish(SolveStatus::NumericalFailure, 0, initial_residual, scale,
                          right_hand_side, solution);
        }
        workspace.direction = workspace.preconditioned_direction;
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        for (int iteration = 1; iteration <= config.max_iterations; ++iteration) {
            workspace.direction_product.noalias() = matrix * workspace.direction;
            const double denominator = workspace.direction.dot(workspace.direction_product);
            if (!std::isfinite(denominator) || denominator <= breakdown_tolerance ||
                std::abs(rho) <= breakdown_tolerance) {
                break;
            }
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
            workspace.direction_product.noalias() = matrix * workspace.preconditioned_direction;
            const double denominator = workspace.shadow.dot(workspace.direction_product);
            if (!std::isfinite(denominator) || std::abs(denominator) <= breakdown_tolerance) break;
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
            workspace.intermediate_product.noalias() =
                matrix * workspace.preconditioned_intermediate;
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

    SolveResult solveGmres(
        const Eigen::VectorXd& right_hand_side,
        Eigen::VectorXd& solution,
        double initial_residual,
        double target,
        double scale)
    {
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        while (iterations < config.max_iterations) {
            workspace.product.noalias() = matrix * solution;
            workspace.residual = right_hand_side - workspace.product;
            const double beta = workspace.residual.norm();
            if (!std::isfinite(beta)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (beta <= target) {
                status = SolveStatus::Converged;
                break;
            }
            workspace.basis.front() = workspace.residual / beta;
            workspace.hessenberg.setZero();
            workspace.least_squares.setZero();
            workspace.least_squares[0] = beta;
            int columns = 0;
            const int cycle_size = std::min(config.gmres_restart,
                                            config.max_iterations - iterations);
            for (int column = 0; column < cycle_size; ++column) {
                if (!precondition(workspace.basis[static_cast<std::size_t>(column)],
                                  workspace.preconditioned_basis[static_cast<std::size_t>(column)])) {
                    status = SolveStatus::NumericalFailure;
                    break;
                }
                workspace.product.noalias() = matrix *
                    workspace.preconditioned_basis[static_cast<std::size_t>(column)];
                for (int row = 0; row <= column; ++row) {
                    const double value = workspace.basis[static_cast<std::size_t>(row)].dot(
                        workspace.product);
                    workspace.hessenberg(row, column) = value;
                    workspace.product.noalias() -= value *
                        workspace.basis[static_cast<std::size_t>(row)];
                }
                workspace.hessenberg(column + 1, column) = workspace.product.norm();
                if (workspace.hessenberg(column + 1, column) > breakdown_tolerance) {
                    workspace.basis[static_cast<std::size_t>(column + 1)] = workspace.product /
                        workspace.hessenberg(column + 1, column);
                }
                for (int row = 0; row < column; ++row) {
                    const double upper = workspace.hessenberg(row, column);
                    const double lower = workspace.hessenberg(row + 1, column);
                    workspace.hessenberg(row, column) =
                        workspace.givens_cosine[row] * upper +
                        workspace.givens_sine[row] * lower;
                    workspace.hessenberg(row + 1, column) =
                        -workspace.givens_sine[row] * upper +
                        workspace.givens_cosine[row] * lower;
                }
                const double upper = workspace.hessenberg(column, column);
                const double lower = workspace.hessenberg(column + 1, column);
                const double norm = std::hypot(upper, lower);
                if (!std::isfinite(norm) || norm <= breakdown_tolerance) {
                    status = SolveStatus::NumericalFailure;
                    break;
                }
                workspace.givens_cosine[column] = upper / norm;
                workspace.givens_sine[column] = lower / norm;
                workspace.hessenberg(column, column) = norm;
                workspace.hessenberg(column + 1, column) = 0.0;
                const double first = workspace.least_squares[column];
                const double second = workspace.least_squares[column + 1];
                workspace.least_squares[column] =
                    workspace.givens_cosine[column] * first +
                    workspace.givens_sine[column] * second;
                workspace.least_squares[column + 1] =
                    -workspace.givens_sine[column] * first +
                    workspace.givens_cosine[column] * second;
                ++iterations;
                columns = column + 1;
                if (std::abs(workspace.least_squares[column + 1]) <= target) break;
            }
            if (columns == 0) break;
            workspace.coefficients.head(columns) =
                workspace.least_squares.head(columns);
            for (int row = columns - 1; row >= 0; --row) {
                workspace.coefficients[row] -= workspace.hessenberg.row(row)
                    .segment(row + 1, columns - row - 1)
                    .dot(workspace.coefficients.segment(row + 1, columns - row - 1));
                workspace.coefficients[row] /= workspace.hessenberg(row, row);
            }
            for (int column = 0; column < columns; ++column) {
                solution.noalias() += workspace.coefficients[column] *
                    workspace.preconditioned_basis[static_cast<std::size_t>(column)];
            }
            if (status == SolveStatus::NumericalFailure) break;
        }
        return finish(status, iterations, initial_residual, scale, right_hand_side, solution);
    }

    LinearSolverConfig config;
    // 求解器拥有快照，调用者可以在 compute/factorize 返回后释放临时矩阵。
    Eigen::SparseMatrix<double> matrix;
    ConjugateGradient conjugate_gradient;
    BiCGSTAB bicgstab;
    std::unique_ptr<detail::AlgebraicMultigrid> amg;
    KrylovWorkspace workspace;
    bool pattern_analyzed = false;
    bool factorization_succeeded = false;
};

void LinearSolverConfig::validate() const {
    if (!(absolute_tolerance > 0.0) || !(relative_tolerance > 0.0) ||
        !std::isfinite(absolute_tolerance) || !std::isfinite(relative_tolerance) ||
        max_iterations <= 0 || gmres_restart <= 0 || amg_max_levels <= 0 ||
        amg_coarse_size <= 0 || amg_smoothing_steps <= 0) {
        throw std::invalid_argument("linear solver configuration is invalid");
    }
    const bool supported =
        (solver == LinearSolverType::ConjugateGradient &&
         preconditioner == PreconditionerType::IncompleteCholesky) ||
        (solver == LinearSolverType::BiCGSTAB &&
         preconditioner == PreconditionerType::ILUT) ||
        ((solver == LinearSolverType::ConjugateGradient ||
          solver == LinearSolverType::BiCGSTAB || solver == LinearSolverType::GMRES) &&
         preconditioner == PreconditionerType::AlgebraicMultigrid) ||
        (solver == LinearSolverType::AlgebraicMultigrid &&
         preconditioner == PreconditionerType::None);
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
PreparedLinearSolver& PreparedLinearSolver::operator=(
    PreparedLinearSolver&&) noexcept = default;

void PreparedLinearSolver::compute(const Eigen::SparseMatrix<double>& A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("linear-system matrix must be square");
    }
    if (!m_implementation) throw std::logic_error("linear solver is moved-from");
    auto& state = *m_implementation;
    state.matrix = A;
    state.pattern_analyzed = false;
    state.factorization_succeeded = false;
    if (usesAmg(state.config)) {
        state.amg = std::make_unique<detail::AlgebraicMultigrid>(state.config);
        state.amg->compute(state.matrix);
        state.workspace.resize(state.matrix.rows(), state.config.gmres_restart);
        state.pattern_analyzed = true;
        state.factorization_succeeded = state.amg->ready();
        return;
    }
    if (state.config.solver == LinearSolverType::ConjugateGradient) {
        state.conjugate_gradient.setMaxIterations(state.config.max_iterations);
        state.conjugate_gradient.compute(state.matrix);
        state.pattern_analyzed = true;
        state.factorization_succeeded = state.conjugate_gradient.info() == Eigen::Success;
        return;
    }
    state.bicgstab.setMaxIterations(state.config.max_iterations);
    state.bicgstab.preconditioner().setDroptol(1e-3);
    state.bicgstab.preconditioner().setFillfactor(2);
    state.bicgstab.compute(state.matrix);
    state.pattern_analyzed = true;
    state.factorization_succeeded = state.bicgstab.info() == Eigen::Success;
}

void PreparedLinearSolver::factorize(const Eigen::SparseMatrix<double>& A) {
    if (!m_implementation) throw std::logic_error("linear solver is moved-from");
    auto& state = *m_implementation;
    if (!state.pattern_analyzed || A.rows() != A.cols()) {
        throw std::logic_error("linear-solver pattern must be analyzed before factorization");
    }
    state.matrix = A;
    state.factorization_succeeded = false;
    if (usesAmg(state.config)) {
        state.amg->factorize(state.matrix);
        state.factorization_succeeded = state.amg->ready();
        return;
    }
    if (state.config.solver == LinearSolverType::ConjugateGradient) {
        state.conjugate_gradient.factorize(state.matrix);
        state.factorization_succeeded = state.conjugate_gradient.info() == Eigen::Success;
        return;
    }
    state.bicgstab.factorize(state.matrix);
    state.factorization_succeeded = state.bicgstab.info() == Eigen::Success;
}

SolveResult PreparedLinearSolver::solve(const Eigen::VectorXd& b, Eigen::VectorXd& x) {
    if (!m_implementation) throw std::logic_error("linear solver is moved-from");
    auto& state = *m_implementation;
    if (!state.pattern_analyzed || state.matrix.rows() != b.size()) {
        throw std::invalid_argument("linear system dimensions are inconsistent");
    }
    if (!state.config.warm_start || x.size() != b.size()) x.setZero(b.size());
    const double initial_residual = (b - state.matrix * x).norm();
    const double scale = std::max({initial_residual, b.norm(), 1e-30});
    const double target = std::max(
        state.config.absolute_tolerance, state.config.relative_tolerance * scale);
    if (initial_residual <= target) {
        return {SolveStatus::Converged, 0, initial_residual, initial_residual,
                initial_residual / scale};
    }
    if (!state.factorization_succeeded) {
        return {SolveStatus::NumericalFailure, 0, initial_residual, initial_residual,
                initial_residual / scale};
    }
    if (!usesAmg(state.config)) {
        return state.config.solver == LinearSolverType::ConjugateGradient
            ? runPreparedSolver(state.conjugate_gradient, state.matrix, b, x, state.config,
                                initial_residual, target)
            : runPreparedSolver(state.bicgstab, state.matrix, b, x, state.config,
                                initial_residual, target);
    }
    if (state.config.solver == LinearSolverType::AlgebraicMultigrid) {
        return state.amg->solve(b, x);
    }
    state.workspace.product.noalias() = state.matrix * x;
    state.workspace.residual = b - state.workspace.product;
    if (state.config.solver == LinearSolverType::ConjugateGradient) {
        return state.solvePcg(b, x, initial_residual, target, scale);
    }
    if (state.config.solver == LinearSolverType::BiCGSTAB) {
        return state.solveBicgstab(b, x, initial_residual, target, scale);
    }
    return state.solveGmres(b, x, initial_residual, target, scale);
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
