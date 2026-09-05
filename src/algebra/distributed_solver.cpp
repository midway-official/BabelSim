#include "internal/mesh_access.h"
#include "babelsim/distributed_solver.h"

#include "backend/algebraic_multigrid.h"

#include <Eigen/IterativeLinearSolvers>

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
    double coefficient = 0.0;
};

bool invalid(double value) {
    return !std::isfinite(value);
}

bool usesAmg(const LinearSolverConfig& config) {
    return config.solver == LinearSolverType::AlgebraicMultigrid ||
        config.preconditioner == PreconditionerType::AlgebraicMultigrid;
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
        mesh.validate();
        config.validate();
        if (!parallel.distributed() ||
            detail::ownedCellCount(mesh) >= mesh.cellCount()) {
            throw std::invalid_argument(
                "distributed solver requires a decomposed multi-rank mesh");
        }
        for (Index face = 0; face < mesh.faceCount(); ++face) {
            const auto f = static_cast<std::size_t>(face);
            const Index owner = detail::meshData(mesh).face_owner[f];
            const Index neighbour = detail::meshData(mesh).face_neighbour[f];
            if (neighbour == invalid_index ||
                detail::isOwned(mesh, owner) == detail::isOwned(mesh, neighbour)) {
                continue;
            }
            if (detail::isOwned(mesh, owner)) {
                remote.push_back({detail::ownedIndex(mesh, owner), neighbour, face, true});
            } else {
                remote.push_back({detail::ownedIndex(mesh, neighbour), owner, face, false});
            }
        }
        const Eigen::Index rows = detail::ownedCellCount(mesh);
        residual.resize(rows);
        matrix_product.resize(rows);
        shadow.resize(rows);
        direction.resize(rows);
        preconditioned_direction.resize(rows);
        direction_product.resize(rows);
        intermediate.resize(rows);
        preconditioned_intermediate.resize(rows);
        intermediate_product.resize(rows);
        resizeGmresWorkspace(rows);
    }

    void resizeGmresWorkspace(Eigen::Index rows) {
        const int restart = config.gmres_restart;
        gmres_basis.resize(static_cast<std::size_t>(restart + 1));
        gmres_preconditioned_basis.resize(static_cast<std::size_t>(restart));
        for (Eigen::VectorXd& value : gmres_basis) value.resize(rows);
        for (Eigen::VectorXd& value : gmres_preconditioned_basis) value.resize(rows);
        gmres_hessenberg.resize(restart + 1, restart);
        gmres_cosine.resize(restart);
        gmres_sine.resize(restart);
        gmres_least_squares.resize(restart + 1);
        gmres_coefficients.resize(restart);
        // 最后一项存放本 rank 的预条件器失败标志，随 Arnoldi 点积一次归约。
        gmres_local_products.resize(static_cast<std::size_t>(restart + 1));
        gmres_global_products.resize(static_cast<std::size_t>(restart + 1));
    }

    void setEquation(
        const Mesh* equation_mesh,
        const std::vector<double>& equation_upper,
        const std::vector<double>& equation_lower)
    {
        if (equation_mesh == nullptr) {
            throw std::invalid_argument("distributed equation has no mesh");
        }
        if (equation_mesh != &mesh ||
            equation_upper.size() != static_cast<std::size_t>(mesh.faceCount()) ||
            equation_lower.size() != static_cast<std::size_t>(mesh.faceCount())) {
            throw std::invalid_argument(
                "distributed equation coefficients do not match the mesh");
        }
        // 只有跨分区面会在 halo 矩阵向量乘中使用系数；本地系数已经由
        // Eigen 稀疏矩阵保存，因此不再为每个外迭代复制整套 LDU 数组。
        for (RemoteCoupling& coupling : remote) {
            const std::size_t f = static_cast<std::size_t>(coupling.face);
            coupling.coefficient = coupling.upper
                ? equation_upper[f] : equation_lower[f];
        }
        equation_ready = true;
    }

    void setMatrix(const Eigen::SparseMatrix<double>& value) {
        if (value.rows() != detail::ownedCellCount(mesh) ||
            value.cols() != detail::ownedCellCount(mesh)) {
            throw std::invalid_argument("distributed local matrix size is invalid");
        }
        matrix = value;
    }

    void computePreconditioner() {
        factorization_succeeded = false;
        if (usesAmg(config)) {
            amg = std::make_unique<detail::AlgebraicMultigrid>(config);
            amg->compute(matrix);
            factorization_succeeded = amg->ready();
        } else if (config.solver == LinearSolverType::ConjugateGradient) {
            incomplete_cholesky.compute(matrix);
            factorization_succeeded =
                incomplete_cholesky.info() == Eigen::Success;
        } else {
            ilut.setDroptol(1e-3);
            ilut.setFillfactor(2);
            ilut.compute(matrix);
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
        if (usesAmg(config)) {
            amg->factorize(matrix);
            factorization_succeeded = amg->ready();
        } else if (config.solver == LinearSolverType::ConjugateGradient) {
            incomplete_cholesky.factorize(matrix);
            factorization_succeeded =
                incomplete_cholesky.info() == Eigen::Success;
        } else {
            ilut.factorize(matrix);
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

    double dotGlobalWithStatus(
        const Eigen::VectorXd& left,
        const Eigen::VectorXd& right,
        bool local_success,
        bool& global_success) const
    {
        const double local[2] = {
            local_success ? left.dot(right) : 0.0,
            local_success ? 0.0 : 1.0,
        };
        double global[2]{};
        parallel.sum(local, global, 2);
        global_success = global[1] == 0.0;
        return global[0];
    }

    void productsGlobalWithStatus(
        double local_first,
        double local_second,
        bool local_success,
        double& global_first,
        double& global_second,
        bool& global_success) const
    {
        const double local[3] = {
            local_success ? local_first : 0.0,
            local_success ? local_second : 0.0,
            local_success ? 0.0 : 1.0,
        };
        double global[3]{};
        parallel.sum(local, global, 3);
        global_first = global[0];
        global_second = global[1];
        global_success = global[2] == 0.0;
    }

    void apply(const Eigen::VectorXd& input, Eigen::VectorXd& output) {
        output.noalias() = matrix * input;
        for (std::size_t owned = 0; owned < detail::meshData(mesh).owned_cells.size(); ++owned) {
            const Index cell = detail::meshData(mesh).owned_cells[owned];
            local_values[static_cast<std::size_t>(cell)] = input[
                static_cast<Eigen::Index>(owned)];
        }
        halo.exchangeFirstLayer(local_values);
        for (const RemoteCoupling& coupling : remote) {
            output[coupling.row] += coupling.coefficient *
                local_values[static_cast<std::size_t>(coupling.ghost_cell)];
        }
    }

    bool precondition(
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output)
    {
        bool local_success = false;
        if (usesAmg(config)) {
            local_success = amg && amg->apply(input, output);
        } else if (config.solver == LinearSolverType::ConjugateGradient) {
            output = incomplete_cholesky.solve(input);
            local_success = incomplete_cholesky.info() == Eigen::Success;
        } else {
            output = ilut.solve(input);
            local_success = ilut.info() == Eigen::Success;
        }
        if (!local_success || !output.allFinite()) {
            // 失败 rank 仍需参加下一次全局归约；零向量避免把 NaN 传播给其他 rank。
            output.setZero();
            return false;
        }
        return true;
    }

    bool preconditionAll(
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output)
    {
        const bool local_success = precondition(input, output);
        const double local_failure = local_success ? 0.0 : 1.0;
        double global_failure = 0.0;
        parallel.sum(&local_failure, &global_failure, 1);
        return global_failure == 0.0;
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
        double scale,
        int iteration_limit)
    {
        const bool local_precondition_success =
            precondition(residual, preconditioned_direction);
        bool global_precondition_success = false;
        double residual_preconditioned = dotGlobalWithStatus(
            residual, preconditioned_direction, local_precondition_success,
            global_precondition_success);
        if (!global_precondition_success) {
            return finish(
                SolveStatus::NumericalFailure, 0, initial_residual, scale, b, x);
        }
        direction = preconditioned_direction;
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        for (int iteration = 1; iteration <= iteration_limit; ++iteration) {
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
            const bool local_success =
                precondition(residual, preconditioned_direction);
            bool global_success = false;
            const double next = dotGlobalWithStatus(
                residual, preconditioned_direction, local_success,
                global_success);
            if (!global_success) {
                status = SolveStatus::NumericalFailure;
                break;
            }
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
        double scale,
        int iteration_limit)
    {
        shadow = residual;
        direction.setZero();
        direction_product.setZero();
        double previous_rho = 1.0;
        double alpha = 1.0;
        double omega = 1.0;
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        for (int iteration = 1; iteration <= iteration_limit; ++iteration) {
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
            const bool local_precondition_success =
                precondition(direction, preconditioned_direction);
            apply(preconditioned_direction, direction_product);
            bool global_precondition_success = false;
            const double shadow_product = dotGlobalWithStatus(
                shadow, direction_product, local_precondition_success,
                global_precondition_success);
            if (!global_precondition_success) {
                status = SolveStatus::NumericalFailure;
                break;
            }
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
            const bool local_intermediate_precondition_success =
                precondition(intermediate, preconditioned_intermediate);
            apply(preconditioned_intermediate, intermediate_product);
            const double local_products[2] = {
                intermediate_product.dot(intermediate),
                intermediate_product.squaredNorm(),
            };
            double global_products[2]{};
            bool global_intermediate_precondition_success = false;
            productsGlobalWithStatus(
                local_products[0], local_products[1],
                local_intermediate_precondition_success,
                global_products[0], global_products[1],
                global_intermediate_precondition_success);
            if (!global_intermediate_precondition_success) {
                status = SolveStatus::NumericalFailure;
                break;
            }
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

    SolveResult solveGmres(
        const Eigen::VectorXd& b,
        Eigen::VectorXd& x,
        double initial_residual,
        double target,
        double scale)
    {
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        while (iterations < config.max_iterations) {
            apply(x, matrix_product);
            residual = b - matrix_product;
            const double beta = normGlobal(residual);
            if (invalid(beta)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (beta <= target) {
                status = SolveStatus::Converged;
                break;
            }
            gmres_basis.front() = residual / beta;
            gmres_hessenberg.setZero();
            gmres_least_squares.setZero();
            gmres_least_squares[0] = beta;
            const int cycle_size = std::min(
                config.gmres_restart, config.max_iterations - iterations);
            int columns = 0;
            for (int column = 0; column < cycle_size; ++column) {
                const bool local_success = precondition(
                    gmres_basis[static_cast<std::size_t>(column)],
                    gmres_preconditioned_basis[static_cast<std::size_t>(column)]);
                apply(gmres_preconditioned_basis[static_cast<std::size_t>(column)],
                      matrix_product);
                for (int row = 0; row <= column; ++row) {
                    gmres_local_products[static_cast<std::size_t>(row)] = local_success
                        ? gmres_basis[static_cast<std::size_t>(row)].dot(matrix_product) : 0.0;
                }
                gmres_local_products[static_cast<std::size_t>(column + 1)] =
                    local_success ? 0.0 : 1.0;
                parallel.sum(gmres_local_products.data(), gmres_global_products.data(),
                             column + 2);
                if (gmres_global_products[static_cast<std::size_t>(column + 1)] != 0.0) {
                    status = SolveStatus::NumericalFailure;
                    break;
                }
                for (int row = 0; row <= column; ++row) {
                    const double value = gmres_global_products[static_cast<std::size_t>(row)];
                    gmres_hessenberg(row, column) = value;
                    matrix_product.noalias() -= value *
                        gmres_basis[static_cast<std::size_t>(row)];
                }
                gmres_hessenberg(column + 1, column) = normGlobal(matrix_product);
                if (invalid(gmres_hessenberg(column + 1, column))) {
                    status = SolveStatus::NumericalFailure;
                    break;
                }
                if (gmres_hessenberg(column + 1, column) > breakdown_tolerance) {
                    gmres_basis[static_cast<std::size_t>(column + 1)] = matrix_product /
                        gmres_hessenberg(column + 1, column);
                }
                for (int row = 0; row < column; ++row) {
                    const double upper = gmres_hessenberg(row, column);
                    const double lower = gmres_hessenberg(row + 1, column);
                    gmres_hessenberg(row, column) = gmres_cosine[row] * upper +
                        gmres_sine[row] * lower;
                    gmres_hessenberg(row + 1, column) = -gmres_sine[row] * upper +
                        gmres_cosine[row] * lower;
                }
                const double upper = gmres_hessenberg(column, column);
                const double lower = gmres_hessenberg(column + 1, column);
                const double hessenberg_norm = std::hypot(upper, lower);
                if (invalid(hessenberg_norm) || hessenberg_norm <= breakdown_tolerance) {
                    status = SolveStatus::NumericalFailure;
                    break;
                }
                gmres_cosine[column] = upper / hessenberg_norm;
                gmres_sine[column] = lower / hessenberg_norm;
                gmres_hessenberg(column, column) = hessenberg_norm;
                gmres_hessenberg(column + 1, column) = 0.0;
                const double first = gmres_least_squares[column];
                const double second = gmres_least_squares[column + 1];
                gmres_least_squares[column] = gmres_cosine[column] * first +
                    gmres_sine[column] * second;
                gmres_least_squares[column + 1] = -gmres_sine[column] * first +
                    gmres_cosine[column] * second;
                ++iterations;
                columns = column + 1;
                if (std::abs(gmres_least_squares[column + 1]) <= target) break;
            }
            if (columns == 0) break;
            gmres_coefficients.head(columns) = gmres_least_squares.head(columns);
            for (int row = columns - 1; row >= 0; --row) {
                gmres_coefficients[row] -= gmres_hessenberg.row(row)
                    .segment(row + 1, columns - row - 1)
                    .dot(gmres_coefficients.segment(row + 1, columns - row - 1));
                gmres_coefficients[row] /= gmres_hessenberg(row, row);
            }
            for (int column = 0; column < columns; ++column) {
                x.noalias() += gmres_coefficients[column] *
                    gmres_preconditioned_basis[static_cast<std::size_t>(column)];
            }
            if (status == SolveStatus::NumericalFailure) break;
        }
        return finish(status, iterations, initial_residual, scale, b, x);
    }

    SolveResult solveAmg(
        const Eigen::VectorXd& b,
        Eigen::VectorXd& x,
        double initial_residual,
        double target,
        double scale)
    {
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        for (; iterations < config.max_iterations; ++iterations) {
            const double residual_norm = normGlobal(residual);
            if (invalid(residual_norm)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (residual_norm <= target) {
                status = SolveStatus::Converged;
                break;
            }
            if (!preconditionAll(residual, preconditioned_direction)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            x += preconditioned_direction;
            apply(x, matrix_product);
            residual = b - matrix_product;
        }
        return finish(status, iterations, initial_residual, scale, b, x);
    }

    const Mesh& mesh;
    ParallelContext parallel;
    LinearSolverConfig config;
    HaloExchange halo;
    // 系数和局部矩阵均由求解器拥有快照，Equation/Assembly 可安全地在调用后销毁。
    Eigen::SparseMatrix<double> matrix;
    std::vector<RemoteCoupling> remote;
    std::vector<double> local_values;
    Eigen::IncompleteCholesky<double> incomplete_cholesky;
    Eigen::IncompleteLUT<double> ilut;
    std::unique_ptr<detail::AlgebraicMultigrid> amg;
    Eigen::VectorXd residual;
    Eigen::VectorXd matrix_product;
    Eigen::VectorXd shadow;
    Eigen::VectorXd direction;
    Eigen::VectorXd preconditioned_direction;
    Eigen::VectorXd direction_product;
    Eigen::VectorXd intermediate;
    Eigen::VectorXd preconditioned_intermediate;
    Eigen::VectorXd intermediate_product;
    // GMRES 的基、预条件基和 Hessenberg 缓冲在构造期一次分配，Krylov 热循环不分配。
    std::vector<Eigen::VectorXd> gmres_basis;
    std::vector<Eigen::VectorXd> gmres_preconditioned_basis;
    Eigen::MatrixXd gmres_hessenberg;
    Eigen::VectorXd gmres_cosine;
    Eigen::VectorXd gmres_sine;
    Eigen::VectorXd gmres_least_squares;
    Eigen::VectorXd gmres_coefficients;
    std::vector<double> gmres_local_products;
    std::vector<double> gmres_global_products;
    bool pattern_ready = false;
    bool factorization_succeeded = false;
    bool equation_ready = false;
};

DistributedLinearSolver::DistributedLinearSolver(
    const Mesh& mesh,
    ParallelContext parallel,
    LinearSolverConfig config)
    : m_implementation(std::make_unique<Implementation>(
          mesh, parallel, std::move(config)))
{}

DistributedLinearSolver::~DistributedLinearSolver() = default;
DistributedLinearSolver::DistributedLinearSolver(
    DistributedLinearSolver&&) noexcept = default;
DistributedLinearSolver& DistributedLinearSolver::operator=(
    DistributedLinearSolver&&) noexcept = default;

void DistributedLinearSolver::compute(
    const Eigen::SparseMatrix<double>& local_matrix,
    const ScalarDiscreteEquation& equation)
{
    if (!m_implementation) throw std::logic_error("distributed solver is moved-from");
    auto& state = *m_implementation;
    equation.validateStorage();
    state.setMatrix(local_matrix);
    state.setEquation(equation.mesh, equation.upper, equation.lower);
    state.computePreconditioner();
}

void DistributedLinearSolver::compute(
    const Eigen::SparseMatrix<double>& local_matrix,
    const VectorDiscreteEquation& equation)
{
    if (!m_implementation) throw std::logic_error("distributed solver is moved-from");
    auto& state = *m_implementation;
    equation.validateStorage();
    state.setMatrix(local_matrix);
    state.setEquation(equation.mesh, equation.upper, equation.lower);
    state.computePreconditioner();
}

void DistributedLinearSolver::factorize(
    const Eigen::SparseMatrix<double>& local_matrix,
    const ScalarDiscreteEquation& equation)
{
    if (!m_implementation) throw std::logic_error("distributed solver is moved-from");
    auto& state = *m_implementation;
    equation.validateStorage();
    state.setMatrix(local_matrix);
    state.setEquation(equation.mesh, equation.upper, equation.lower);
    state.factorizePreconditioner();
}

void DistributedLinearSolver::factorize(
    const Eigen::SparseMatrix<double>& local_matrix,
    const VectorDiscreteEquation& equation)
{
    if (!m_implementation) throw std::logic_error("distributed solver is moved-from");
    auto& state = *m_implementation;
    equation.validateStorage();
    state.setMatrix(local_matrix);
    state.setEquation(equation.mesh, equation.upper, equation.lower);
    state.factorizePreconditioner();
}

SolveResult DistributedLinearSolver::solve(
    const Eigen::VectorXd& b,
    Eigen::VectorXd& x)
{
    if (!m_implementation) throw std::logic_error("distributed solver is moved-from");
    auto& state = *m_implementation;
    if (!state.pattern_ready || !state.equation_ready || state.matrix.rows() == 0 ||
        b.size() != detail::ownedCellCount(state.mesh)) {
        throw std::invalid_argument("distributed linear system is not prepared");
    }
    if (!state.config.warm_start || x.size() != b.size()) {
        x = Eigen::VectorXd::Zero(b.size());
    }
    state.apply(x, state.matrix_product);
    state.residual = b - state.matrix_product;
    const double local_norms[2] = {
        state.residual.squaredNorm(), b.squaredNorm()};
    double global_norms[2]{};
    state.parallel.sum(local_norms, global_norms, 2);
    const double initial_residual = std::sqrt(std::max(global_norms[0], 0.0));
    const double rhs_norm = std::sqrt(std::max(global_norms[1], 0.0));
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
    if (state.config.solver == LinearSolverType::AlgebraicMultigrid) {
        return state.solveAmg(b, x, initial_residual, target, scale);
    }
    if (state.config.solver == LinearSolverType::GMRES) {
        return state.solveGmres(b, x, initial_residual, target, scale);
    }
    // 有限精度下可能提前停在近似 breakdown，或递推残差与真实残差不一致。
    // 有进展且预算尚有剩余时，以真实残差重启 Krylov；不增加外迭代，不放宽容差。
    // 所有判据均来自全局范数，所有 rank 的分支及剩余预算完全一致。
    SolveResult result;
    int completed = 0;
    double previous_residual = initial_residual;
    while (completed < state.config.max_iterations) {
        const int remaining = state.config.max_iterations - completed;
        result = state.config.solver == LinearSolverType::ConjugateGradient
            ? state.solvePcg(b, x, initial_residual, target, scale, remaining)
            : state.solveBicgstab(b, x, initial_residual, target, scale, remaining);
        const int used = result.iterations;
        completed += used;
        result.iterations = completed;
        if (result.converged() || !result.healthy() || used == 0 ||
            result.final_residual >= previous_residual * (1.0 - 1e-8)) break;
        // finish() 已计算 A*x，因此不需要为重启再做一次 halo matvec。
        state.residual = b - state.matrix_product;
        previous_residual = result.final_residual;
    }
    return result;
}

}  // babelsim 命名空间
