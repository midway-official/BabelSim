#include "babelsim/distributed_solver.h"

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace babelsim {
namespace {

constexpr double breakdown_tolerance = 1e-30;

struct RemoteCoupling {
    Index row = invalid_index;
    Index ghost_cell = invalid_index;
    Index face = invalid_index;
    bool upper = false;
};

bool invalid(double value) {
    return !std::isfinite(value);
}

}  // 匿名命名空间

struct DistributedLinearSolver::Implementation {
    Implementation(
        const Mesh& mesh_value,
        ParallelContext parallel_value,
        LinearSolverConfig config_value)
        : mesh(mesh_value),
          parallel(parallel_value),
          config(std::move(config_value)),
          halo(mesh, parallel),
          local_values(static_cast<std::size_t>(mesh.cellCount()), 0.0)
    {
        parallel.validate();
        config.validate();
        if (!parallel.distributed() ||
            mesh.ownedCellCount() >= mesh.cellCount()) {
            throw std::invalid_argument(
                "distributed solver requires a decomposed multi-rank mesh");
        }
        for (Index face = 0; face < mesh.faceCount(); ++face) {
            const auto f = static_cast<std::size_t>(face);
            const Index owner = mesh.face_owner[f];
            const Index neighbour = mesh.face_neighbour[f];
            if (neighbour == invalid_index ||
                mesh.isOwned(owner) == mesh.isOwned(neighbour)) {
                continue;
            }
            if (mesh.isOwned(owner)) {
                remote.push_back({mesh.ownedIndex(owner), neighbour, face, true});
            } else {
                remote.push_back({mesh.ownedIndex(neighbour), owner, face, false});
            }
        }
        const Eigen::Index rows = mesh.ownedCellCount();
        residual.resize(rows);
        matrix_product.resize(rows);
        shadow.resize(rows);
        direction.resize(rows);
        preconditioned_direction.resize(rows);
        direction_product.resize(rows);
        intermediate.resize(rows);
        preconditioned_intermediate.resize(rows);
        intermediate_product.resize(rows);
    }

    void setEquation(
        const Mesh* equation_mesh,
        const std::vector<double>& equation_upper,
        const std::vector<double>& equation_lower)
    {
        if (equation_mesh != &mesh ||
            equation_upper.size() != static_cast<std::size_t>(mesh.faceCount()) ||
            equation_lower.size() != static_cast<std::size_t>(mesh.faceCount())) {
            throw std::invalid_argument(
                "distributed equation coefficients do not match the mesh");
        }
        upper = &equation_upper;
        lower = &equation_lower;
    }

    void setMatrix(const Eigen::SparseMatrix<double>& value) {
        if (value.rows() != mesh.ownedCellCount() ||
            value.cols() != mesh.ownedCellCount()) {
            throw std::invalid_argument("distributed local matrix size is invalid");
        }
        matrix = &value;
    }

    void computePreconditioner() {
        factorization_succeeded = false;
        if (config.solver == LinearSolverType::ConjugateGradient) {
            incomplete_cholesky.compute(*matrix);
            factorization_succeeded =
                incomplete_cholesky.info() == Eigen::Success;
        } else {
            ilut.setDroptol(1e-3);
            ilut.setFillfactor(2);
            ilut.compute(*matrix);
            factorization_succeeded = ilut.info() == Eigen::Success;
        }
        factorization_succeeded =
            parallel.maximum(factorization_succeeded ? 0 : 1) == 0;
        pattern_ready = true;
    }

    void factorizePreconditioner() {
        if (!pattern_ready) {
            throw std::logic_error(
                "distributed pattern must be computed before factorization");
        }
        factorization_succeeded = false;
        if (config.solver == LinearSolverType::ConjugateGradient) {
            incomplete_cholesky.factorize(*matrix);
            factorization_succeeded =
                incomplete_cholesky.info() == Eigen::Success;
        } else {
            ilut.factorize(*matrix);
            factorization_succeeded = ilut.info() == Eigen::Success;
        }
        factorization_succeeded =
            parallel.maximum(factorization_succeeded ? 0 : 1) == 0;
    }

    double dotGlobal(
        const Eigen::VectorXd& left,
        const Eigen::VectorXd& right) const
    {
        const double local = left.dot(right);
        double global = 0.0;
        parallel.sum(&local, &global, 1);
        return global;
    }

    double normGlobal(const Eigen::VectorXd& value) const {
        return std::sqrt(std::max(dotGlobal(value, value), 0.0));
    }

    void apply(const Eigen::VectorXd& input, Eigen::VectorXd& output) {
        output.noalias() = (*matrix) * input;
        for (Index cell : mesh.owned_cells) {
            local_values[static_cast<std::size_t>(cell)] =
                input[mesh.ownedIndex(cell)];
        }
        halo.exchange(local_values);
        for (const RemoteCoupling& coupling : remote) {
            const auto f = static_cast<std::size_t>(coupling.face);
            const double coefficient = coupling.upper
                ? (*upper)[f] : (*lower)[f];
            output[coupling.row] += coefficient *
                local_values[static_cast<std::size_t>(coupling.ghost_cell)];
        }
    }

    bool precondition(
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output)
    {
        bool local_success = false;
        if (config.solver == LinearSolverType::ConjugateGradient) {
            output = incomplete_cholesky.solve(input);
            local_success = incomplete_cholesky.info() == Eigen::Success;
        } else {
            output = ilut.solve(input);
            local_success = ilut.info() == Eigen::Success;
        }
        return parallel.maximum(local_success ? 0 : 1) == 0;
    }

    SolveResult finish(
        SolveStatus status,
        int iterations,
        double initial_residual,
        double scale,
        const Eigen::VectorXd& b,
        const Eigen::VectorXd& x)
    {
        apply(x, matrix_product);
        const double final_residual = normGlobal(b - matrix_product);
        const double target = std::max(
            config.absolute_tolerance,
            config.relative_tolerance * scale);
        if (!std::isfinite(final_residual)) {
            status = SolveStatus::NumericalFailure;
        } else if (final_residual <= target * (1.0 + 1e-8)) {
            status = SolveStatus::Converged;
        } else if (status == SolveStatus::Converged) {
            // 递推 Krylov 残差可能偏离真实残差，因此周期性计算实际残差。
            status = SolveStatus::MaxIterations;
        }
        return {
            status,
            iterations,
            initial_residual,
            final_residual,
            final_residual / scale,
        };
    }

    SolveResult solvePcg(
        const Eigen::VectorXd& b,
        Eigen::VectorXd& x,
        double initial_residual,
        double target,
        double scale)
    {
        if (!precondition(residual, preconditioned_direction)) {
            return finish(
                SolveStatus::NumericalFailure, 0, initial_residual, scale, b, x);
        }
        direction = preconditioned_direction;
        double residual_preconditioned =
            dotGlobal(residual, preconditioned_direction);
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        for (int iteration = 1; iteration <= config.max_iterations; ++iteration) {
            apply(direction, direction_product);
            const double denominator = dotGlobal(direction, direction_product);
            if (invalid(denominator) || invalid(residual_preconditioned)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (denominator <= breakdown_tolerance) {
                status = SolveStatus::MaxIterations;
                break;
            }
            const double alpha = residual_preconditioned / denominator;
            x.noalias() += alpha * direction;
            residual.noalias() -= alpha * direction_product;
            iterations = iteration;
            const double residual_norm = normGlobal(residual);
            if (invalid(residual_norm)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (residual_norm <= target) {
                status = SolveStatus::Converged;
                break;
            }
            if (!precondition(residual, preconditioned_direction)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            const double next = dotGlobal(residual, preconditioned_direction);
            if (invalid(next)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (std::abs(residual_preconditioned) <= breakdown_tolerance) {
                status = SolveStatus::MaxIterations;
                break;
            }
            direction = preconditioned_direction +
                (next / residual_preconditioned) * direction;
            residual_preconditioned = next;
        }
        return finish(status, iterations, initial_residual, scale, b, x);
    }

    SolveResult solveBicgstab(
        const Eigen::VectorXd& b,
        Eigen::VectorXd& x,
        double initial_residual,
        double target,
        double scale)
    {
        shadow = residual;
        direction.setZero();
        direction_product.setZero();
        double previous_rho = 1.0;
        double alpha = 1.0;
        double omega = 1.0;
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        for (int iteration = 1; iteration <= config.max_iterations; ++iteration) {
            const double rho = dotGlobal(shadow, residual);
            if (invalid(rho) || invalid(omega)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (std::abs(rho) <= breakdown_tolerance ||
                std::abs(omega) <= breakdown_tolerance) {
                status = SolveStatus::MaxIterations;
                break;
            }
            const double beta = (rho / previous_rho) * (alpha / omega);
            direction = residual + beta * (direction - omega * direction_product);
            if (!precondition(direction, preconditioned_direction)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            apply(preconditioned_direction, direction_product);
            const double shadow_product = dotGlobal(shadow, direction_product);
            if (invalid(shadow_product)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (std::abs(shadow_product) <= breakdown_tolerance) {
                status = SolveStatus::MaxIterations;
                break;
            }
            alpha = rho / shadow_product;
            intermediate = residual - alpha * direction_product;
            const double intermediate_norm = normGlobal(intermediate);
            iterations = iteration;
            if (invalid(intermediate_norm)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (intermediate_norm <= target) {
                x.noalias() += alpha * preconditioned_direction;
                residual = intermediate;
                status = SolveStatus::Converged;
                break;
            }
            if (!precondition(intermediate, preconditioned_intermediate)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            apply(preconditioned_intermediate, intermediate_product);
            const double local_products[2] = {
                intermediate_product.dot(intermediate),
                intermediate_product.squaredNorm(),
            };
            double global_products[2]{};
            parallel.sum(local_products, global_products, 2);
            if (invalid(global_products[0]) || invalid(global_products[1])) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (global_products[1] <= breakdown_tolerance) {
                status = SolveStatus::MaxIterations;
                break;
            }
            omega = global_products[0] / global_products[1];
            if (invalid(omega)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (std::abs(omega) <= breakdown_tolerance) {
                status = SolveStatus::MaxIterations;
                break;
            }
            x.noalias() += alpha * preconditioned_direction +
                omega * preconditioned_intermediate;
            residual = intermediate - omega * intermediate_product;
            const double residual_norm = normGlobal(residual);
            if (invalid(residual_norm)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (residual_norm <= target) {
                status = SolveStatus::Converged;
                break;
            }
            previous_rho = rho;
        }
        return finish(status, iterations, initial_residual, scale, b, x);
    }

    const Mesh& mesh;
    ParallelContext parallel;
    LinearSolverConfig config;
    HaloExchange halo;
    const Eigen::SparseMatrix<double>* matrix = nullptr;
    const std::vector<double>* upper = nullptr;
    const std::vector<double>* lower = nullptr;
    std::vector<RemoteCoupling> remote;
    std::vector<double> local_values;
    Eigen::IncompleteCholesky<double> incomplete_cholesky;
    Eigen::IncompleteLUT<double> ilut;
    Eigen::VectorXd residual;
    Eigen::VectorXd matrix_product;
    Eigen::VectorXd shadow;
    Eigen::VectorXd direction;
    Eigen::VectorXd preconditioned_direction;
    Eigen::VectorXd direction_product;
    Eigen::VectorXd intermediate;
    Eigen::VectorXd preconditioned_intermediate;
    Eigen::VectorXd intermediate_product;
    bool pattern_ready = false;
    bool factorization_succeeded = false;
};

DistributedLinearSolver::DistributedLinearSolver(
    const Mesh& mesh,
    ParallelContext parallel,
    LinearSolverConfig config)
    : implementation_(std::make_unique<Implementation>(
          mesh, parallel, std::move(config)))
{}

DistributedLinearSolver::~DistributedLinearSolver() = default;
DistributedLinearSolver::DistributedLinearSolver(
    DistributedLinearSolver&&) noexcept = default;
DistributedLinearSolver& DistributedLinearSolver::operator=(
    DistributedLinearSolver&&) noexcept = default;

void DistributedLinearSolver::compute(
    const Eigen::SparseMatrix<double>& local_matrix,
    const ScalarEquation& equation)
{
    auto& state = *implementation_;
    state.setMatrix(local_matrix);
    state.setEquation(equation.mesh, equation.upper, equation.lower);
    state.computePreconditioner();
}

void DistributedLinearSolver::compute(
    const Eigen::SparseMatrix<double>& local_matrix,
    const VectorEquation& equation)
{
    auto& state = *implementation_;
    state.setMatrix(local_matrix);
    state.setEquation(equation.mesh, equation.upper, equation.lower);
    state.computePreconditioner();
}

void DistributedLinearSolver::factorize(
    const Eigen::SparseMatrix<double>& local_matrix,
    const ScalarEquation& equation)
{
    auto& state = *implementation_;
    state.setMatrix(local_matrix);
    state.setEquation(equation.mesh, equation.upper, equation.lower);
    state.factorizePreconditioner();
}

void DistributedLinearSolver::factorize(
    const Eigen::SparseMatrix<double>& local_matrix,
    const VectorEquation& equation)
{
    auto& state = *implementation_;
    state.setMatrix(local_matrix);
    state.setEquation(equation.mesh, equation.upper, equation.lower);
    state.factorizePreconditioner();
}

SolveResult DistributedLinearSolver::solve(
    const Eigen::VectorXd& b,
    Eigen::VectorXd& x)
{
    auto& state = *implementation_;
    if (state.matrix == nullptr || state.upper == nullptr ||
        b.size() != state.mesh.ownedCellCount()) {
        throw std::invalid_argument("distributed linear system is not prepared");
    }
    if (!state.config.warm_start || x.size() != b.size()) {
        x = Eigen::VectorXd::Zero(b.size());
    }
    state.apply(x, state.matrix_product);
    state.residual = b - state.matrix_product;
    const double initial_residual = state.normGlobal(state.residual);
    const double rhs_norm = state.normGlobal(b);
    const double scale = std::max({initial_residual, rhs_norm, 1e-30});
    const double target = std::max(
        state.config.absolute_tolerance,
        state.config.relative_tolerance * scale);
    if (initial_residual <= target) {
        return {
            SolveStatus::Converged, 0, initial_residual,
            initial_residual, initial_residual / scale,
        };
    }
    if (!state.factorization_succeeded) {
        return {
            SolveStatus::NumericalFailure, 0, initial_residual,
            initial_residual, initial_residual / scale,
        };
    }
    if (state.config.solver == LinearSolverType::ConjugateGradient) {
        return state.solvePcg(b, x, initial_residual, target, scale);
    }
    return state.solveBicgstab(b, x, initial_residual, target, scale);
}

}  // babelsim 命名空间
