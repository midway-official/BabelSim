#include "internal/mesh_access.h"
#include "babelsim/distributed_solver.h"

#include "babelsim/mpi_support.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace babelsim {
namespace {

constexpr double breakdown_tolerance = 1e-30;
using Clock = std::chrono::steady_clock;

double secondsSince(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct RemoteCoupling {
    Index row = invalid_index;
    Index ghost_cell = invalid_index;
    Index face = invalid_index;
    bool upper = false;
    double coefficient = 0.0;
};

// Krylov 向量只包含 owned 行。把它投影到局部 cell 布局后复用拓扑无关的
// HaloExchange；这使远程矩阵耦合只依赖 global cell ID，而非空间方向。
class KrylovHalo {
public:
    KrylovHalo(const Mesh& mesh, ParallelContext parallel)
        : m_mesh(&mesh), m_exchange(mesh, parallel),
          m_values(static_cast<std::size_t>(mesh.cellCount()), 0.0) {}

    void begin(const Eigen::VectorXd& owned_values) {
        if (m_active) throw std::logic_error("Krylov halo exchange is already active");
        if (owned_values.size() != detail::ownedCellCount(*m_mesh)) {
            throw std::invalid_argument("Krylov vector does not match owned cells");
        }
        std::fill(m_values.begin(), m_values.end(), 0.0);
        for (Index cell : detail::meshData(*m_mesh).owned_cells) {
            m_values[static_cast<std::size_t>(cell)] = owned_values[detail::ownedIndex(*m_mesh, cell)];
        }
        m_exchange.exchangeFirstLayer(m_values);
        m_active = true;
    }

    void finish() {
        if (!m_active) throw std::logic_error("Krylov halo exchange is not active");
        m_active = false;
    }

    double value(Index ghost_cell) const {
        return m_values.at(static_cast<std::size_t>(ghost_cell));
    }

private:
    const Mesh* m_mesh;
    HaloExchange m_exchange;
    std::vector<double> m_values;
    bool m_active = false;
};

bool invalid(double value) {
    return !std::isfinite(value);
}

bool usesAmg(const LinearSolverConfig& config) {
    return config.preconditioner == PreconditionerType::AlgebraicMultigrid;
}

bool hasPreconditioner(const LinearSolverConfig& config) {
    return config.preconditioner != PreconditionerType::None;
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
          krylov_halo(mesh, parallel)
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
        if (!spmv_pattern_ready) {
            boundary_rows.assign(static_cast<std::size_t>(value.rows()), 0);
            for (const RemoteCoupling& coupling : remote) {
                boundary_rows[static_cast<std::size_t>(coupling.row)] = 1;
            }
            std::vector<Eigen::Triplet<double>> interior_entries;
            std::vector<Eigen::Triplet<double>> boundary_entries;
            interior_entries.reserve(static_cast<std::size_t>(value.nonZeros()));
            boundary_entries.reserve(remote.size() * 8U);
            for (Eigen::Index column = 0; column < value.outerSize(); ++column) {
                for (Eigen::SparseMatrix<double>::InnerIterator entry(value, column);
                     entry; ++entry) {
                    auto& entries = boundary_rows[static_cast<std::size_t>(entry.row())]
                        ? boundary_entries : interior_entries;
                    entries.emplace_back(entry.row(), entry.col(), entry.value());
                }
            }
            interior_matrix.resize(value.rows(), value.cols());
            boundary_matrix.resize(value.rows(), value.cols());
            interior_matrix.setFromTriplets(interior_entries.begin(), interior_entries.end());
            boundary_matrix.setFromTriplets(boundary_entries.begin(), boundary_entries.end());
            interior_matrix.makeCompressed();
            boundary_matrix.makeCompressed();
            spmv_pattern_ready = true;
        } else {
            // 稀疏模式固定时只覆盖已有系数，不在每个外迭代重新分配 Triplet/矩阵。
            for (Eigen::Index column = 0; column < value.outerSize(); ++column) {
                for (Eigen::SparseMatrix<double>::InnerIterator entry(value, column);
                     entry; ++entry) {
                    Eigen::SparseMatrix<double>& target =
                        boundary_rows[static_cast<std::size_t>(entry.row())]
                        ? boundary_matrix : interior_matrix;
                    target.coeffRef(entry.row(), entry.col()) = entry.value();
                }
            }
        }
    }

    void computePreconditioner() {
        factorization_succeeded = false;
        if (!hasPreconditioner(config)) {
            factorization_succeeded = true;
        } else if (usesAmg(config)) {
            factorization_succeeded = updateDistributedAmg(true);
            amg_updates_since_factorization = 0;
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
        if (hasPreconditioner(config)) {
            factorization_succeeded = globallyReady(factorization_succeeded);
        }
        pattern_ready = true;
    }

    void factorizePreconditioner() {
        if (!pattern_ready) {
            throw std::logic_error(
                "distributed pattern must be computed before factorization");
        }
        factorization_succeeded = false;
        if (!hasPreconditioner(config)) {
            factorization_succeeded = true;
        } else if (usesAmg(config)) {
            const bool refresh =
                ++amg_updates_since_factorization >= config.amg_refresh_interval;
            if (refresh) {
                factorization_succeeded = updateDistributedAmg(false);
                amg_updates_since_factorization = 0;
            } else {
                factorization_succeeded = amg_ready;
            }
        } else if (config.solver == LinearSolverType::ConjugateGradient) {
            incomplete_cholesky.factorize(matrix);
            factorization_succeeded =
                incomplete_cholesky.info() == Eigen::Success;
        } else {
            ilut.factorize(matrix);
            factorization_succeeded = ilut.info() == Eigen::Success;
        }
        if (hasPreconditioner(config)) {
            factorization_succeeded = globallyReady(factorization_succeeded);
        }
    }

    bool globallyReady(bool local_ready) {
        const Clock::time_point start = Clock::now();
        const bool ready = parallel.maximum(local_ready ? 0 : 1) == 0;
        ++pending_performance.global_reductions;
        pending_performance.global_reduction_seconds += secondsSince(start);
        return ready;
    }

    void sumGlobal(const double* local, double* global, int count) const {
        const Clock::time_point start = Clock::now();
        parallel.sum(local, global, count);
        ++current_performance.global_reductions;
        current_performance.global_reduction_seconds += secondsSince(start);
    }

    int coarseIndex(Index global_cell) const {
        if (global_cell < 0 || static_cast<std::size_t>(global_cell) >= amg_cell_to_coarse.size()) {
            throw std::out_of_range("AMG aggregate cell id is invalid");
        }
        return amg_cell_to_coarse[static_cast<std::size_t>(global_cell)];
    }

    void buildCoarseMapping() {
        const Index global_cells = mesh.globalCellCount();
        std::vector<int> adjacency(static_cast<std::size_t>(global_cells) * 6U, invalid_index);
        // 每个 owned cell 恰由一个 rank 发布其图邻接。MAX 归约在 -1 哨兵与合法
        // 非负 global ID 之间得到唯一结果，随后每个 rank 都拥有相同的粗化图。
        for (Index cell : detail::meshData(mesh).owned_cells) {
            const Index global_cell = detail::globalCellId(mesh, cell);
            for (Index slot = 0; slot < 6; ++slot) {
                const Index neighbour = detail::meshData(mesh).cell_neighbours[static_cast<std::size_t>(cell)]
                    [static_cast<std::size_t>(slot)];
                adjacency[6U * static_cast<std::size_t>(global_cell) + static_cast<std::size_t>(slot)] =
                    neighbour == invalid_index ? invalid_index : detail::globalCellId(mesh, neighbour);
            }
        }
        std::vector<int> global_adjacency(adjacency.size(), invalid_index);
        const Clock::time_point adjacency_start = Clock::now();
        detail::checkMpi(MPI_Allreduce(
            adjacency.data(), global_adjacency.data(),
            detail::mpiCount(adjacency.size(), "AMG graph adjacency"), MPI_INT, MPI_MAX,
            parallel.communicator), "MPI_Allreduce(AMG graph adjacency)");
        ++pending_performance.global_reductions;
        pending_performance.global_reduction_seconds += secondsSince(adjacency_start);

        amg_cell_to_coarse.resize(static_cast<std::size_t>(global_cells));
        for (Index cell = 0; cell < global_cells; ++cell) {
            amg_cell_to_coarse[static_cast<std::size_t>(cell)] = cell;
        }
        int groups = global_cells;
        int levels = 1;
        while (groups > config.amg_coarse_size && levels < config.amg_max_levels) {
            std::vector<std::vector<int>> graph(static_cast<std::size_t>(groups));
            for (Index cell = 0; cell < global_cells; ++cell) {
                const int source = amg_cell_to_coarse[static_cast<std::size_t>(cell)];
                for (Index slot = 0; slot < 6; ++slot) {
                    const int neighbour = global_adjacency[
                        6U * static_cast<std::size_t>(cell) + static_cast<std::size_t>(slot)];
                    if (neighbour == invalid_index) continue;
                    const int target = amg_cell_to_coarse[static_cast<std::size_t>(neighbour)];
                    if (target != source) graph[static_cast<std::size_t>(source)].push_back(target);
                }
            }
            for (auto& links : graph) {
                std::sort(links.begin(), links.end());
                links.erase(std::unique(links.begin(), links.end()), links.end());
            }
            std::vector<int> next(static_cast<std::size_t>(groups), invalid_index);
            int next_groups = 0;
            for (int group = 0; group < groups; ++group) {
                if (next[static_cast<std::size_t>(group)] != invalid_index) continue;
                next[static_cast<std::size_t>(group)] = next_groups;
                for (int neighbour : graph[static_cast<std::size_t>(group)]) {
                    if (next[static_cast<std::size_t>(neighbour)] == invalid_index) {
                        next[static_cast<std::size_t>(neighbour)] = next_groups;
                        break;
                    }
                }
                ++next_groups;
            }
            if (next_groups >= groups) break;
            for (int& aggregate : amg_cell_to_coarse) {
                aggregate = next[static_cast<std::size_t>(aggregate)];
            }
            groups = next_groups;
            ++levels;
        }
        const std::int64_t coarse_count = groups;
        // 当前全局粗矩阵在每个 rank 复制，以避免 root 串行通信。限制实际行数，
        // 防止错误配置把 coarse_count² 的准备缓冲膨胀为不可控内存。
        constexpr std::int64_t maximum_replicated_coarse_rows = 2048;
        if (coarse_count <= 0 || coarse_count > maximum_replicated_coarse_rows) {
            throw std::runtime_error("distributed AMG coarse space is invalid");
        }
        amg_aggregate.resize(static_cast<std::size_t>(matrix.rows()));
        const auto& owned = detail::meshData(mesh).owned_cells;
        for (std::size_t row = 0; row < owned.size(); ++row) {
            amg_aggregate[row] = coarseIndex(detail::globalCellId(mesh, owned[row]));
        }
        amg_coarse_rhs.resize(coarse_count);
        amg_global_rhs.resize(coarse_count);
        amg_coarse_correction.resize(coarse_count);
        amg_local_coarse.assign(
            static_cast<std::size_t>(coarse_count * coarse_count), 0.0);
        amg_global_coarse.resize(amg_local_coarse.size());
    }

    bool updateDistributedAmg(bool rebuild_pattern) {
        if (rebuild_pattern) buildCoarseMapping();
        amg_inverse_diagonal.resize(matrix.rows());
        bool diagonal_ok = true;
        for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
            const double diagonal = matrix.coeff(row, row);
            if (!std::isfinite(diagonal) || std::abs(diagonal) <= breakdown_tolerance) {
                diagonal_ok = false;
                amg_inverse_diagonal[row] = 0.0;
            } else {
                amg_inverse_diagonal[row] = 1.0 / diagonal;
            }
        }
        const Clock::time_point diagonal_reduction_start = Clock::now();
        const bool global_diagonal_ok =
            parallel.maximum(diagonal_ok ? 0 : 1) == 0;
        ++pending_performance.global_reductions;
        pending_performance.global_reduction_seconds +=
            secondsSince(diagonal_reduction_start);
        if (!global_diagonal_ok) {
            amg_ready = false;
            return false;
        }
        std::fill(amg_local_coarse.begin(), amg_local_coarse.end(), 0.0);
        const std::size_t coarse_count = static_cast<std::size_t>(amg_coarse_rhs.size());
        const auto add = [&](int row, int column, double value) {
            amg_local_coarse[static_cast<std::size_t>(row) * coarse_count +
                             static_cast<std::size_t>(column)] += value;
        };
        for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
            for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix, column);
                 entry; ++entry) {
                add(amg_aggregate[static_cast<std::size_t>(entry.row())],
                    amg_aggregate[static_cast<std::size_t>(entry.col())], entry.value());
            }
        }
        for (const RemoteCoupling& coupling : remote) {
            add(amg_aggregate[static_cast<std::size_t>(coupling.row)],
                coarseIndex(detail::globalCellId(mesh, coupling.ghost_cell)),
                coupling.coefficient);
        }
        const Clock::time_point reduction_start = Clock::now();
        parallel.sum(amg_local_coarse.data(), amg_global_coarse.data(),
                     detail::mpiCount(amg_local_coarse.size(), "AMG coarse matrix"));
        ++pending_performance.global_reductions;
        pending_performance.global_reduction_seconds += secondsSince(reduction_start);

        std::vector<Eigen::Triplet<double>> entries;
        entries.reserve(coarse_count * 7U);
        for (std::size_t row = 0; row < coarse_count; ++row) {
            for (std::size_t column = 0; column < coarse_count; ++column) {
                entries.emplace_back(
                    row, column, amg_global_coarse[row * coarse_count + column]);
            }
        }
        amg_coarse_matrix.resize(coarse_count, coarse_count);
        amg_coarse_matrix.setFromTriplets(entries.begin(), entries.end());
        amg_coarse_matrix.makeCompressed();
        if (rebuild_pattern) amg_coarse_solver.analyzePattern(amg_coarse_matrix);
        amg_coarse_solver.factorize(amg_coarse_matrix);
        amg_residual.resize(matrix.rows());
        amg_product.resize(matrix.rows());
        amg_ready = amg_coarse_solver.info() == Eigen::Success;
        return amg_ready;
    }

    bool applyDistributedAmg(
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output)
    {
        if (!amg_ready || &input == &output) return false;
        constexpr double weight = 2.0 / 3.0;
        output.noalias() = weight * amg_inverse_diagonal.cwiseProduct(input);
        for (int sweep = 1; sweep < config.amg_smoothing_steps; ++sweep) {
            apply(output, amg_product);
            amg_residual = input - amg_product;
            output.noalias() += weight *
                amg_inverse_diagonal.cwiseProduct(amg_residual);
        }
        apply(output, amg_product);
        amg_residual = input - amg_product;
        amg_coarse_rhs.setZero();
        for (Eigen::Index row = 0; row < amg_residual.size(); ++row) {
            amg_coarse_rhs[amg_aggregate[static_cast<std::size_t>(row)]] +=
                amg_residual[row];
        }
        sumGlobal(amg_coarse_rhs.data(), amg_global_rhs.data(),
                  static_cast<int>(amg_coarse_rhs.size()));
        amg_coarse_correction = amg_coarse_solver.solve(amg_global_rhs);
        if (amg_coarse_solver.info() != Eigen::Success ||
            !amg_coarse_correction.allFinite()) return false;
        for (Eigen::Index row = 0; row < output.size(); ++row) {
            output[row] += amg_coarse_correction[
                amg_aggregate[static_cast<std::size_t>(row)]];
        }
        for (int sweep = 0; sweep < config.amg_smoothing_steps; ++sweep) {
            apply(output, amg_product);
            amg_residual = input - amg_product;
            output.noalias() += weight *
                amg_inverse_diagonal.cwiseProduct(amg_residual);
        }
        return output.allFinite();
    }

    double dotGlobal(
        const Eigen::VectorXd& left,
        const Eigen::VectorXd& right) const
    {
        const double local = left.dot(right);
        double global = 0.0;
        sumGlobal(&local, &global, 1);
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
        sumGlobal(local, global, 2);
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
        sumGlobal(local, global, 3);
        global_first = global[0];
        global_second = global[1];
        global_success = global[2] == 0.0;
    }

    void apply(const Eigen::VectorXd& input, Eigen::VectorXd& output) {
        const Clock::time_point start = Clock::now();
        ++current_performance.sparse_matvecs;
        const Clock::time_point halo_start = Clock::now();
        krylov_halo.begin(input);
        current_performance.halo_seconds += secondsSince(halo_start);
        output.noalias() = interior_matrix * input;
        const Clock::time_point halo_wait_start = Clock::now();
        krylov_halo.finish();
        ++current_performance.halo_exchanges;
        current_performance.halo_seconds += secondsSince(halo_wait_start);
        output.noalias() += boundary_matrix * input;
        for (const RemoteCoupling& coupling : remote) {
            output[coupling.row] += coupling.coefficient *
                krylov_halo.value(coupling.ghost_cell);
        }
        current_performance.sparse_matvec_seconds += secondsSince(start);
    }

    bool precondition(
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output)
    {
        const Clock::time_point start = Clock::now();
        bool local_success = false;
        if (!hasPreconditioner(config)) {
            output = input;
            local_success = true;
        } else if (usesAmg(config)) {
            local_success = applyDistributedAmg(input, output);
        } else if (config.solver == LinearSolverType::ConjugateGradient) {
            output = incomplete_cholesky.solve(input);
            local_success = incomplete_cholesky.info() == Eigen::Success;
        } else {
            output = ilut.solve(input);
            local_success = ilut.info() == Eigen::Success;
        }
        if (hasPreconditioner(config)) {
            ++current_performance.preconditioner_applications;
            current_performance.preconditioner_apply_seconds += secondsSince(start);
        }
        if (!local_success || !output.allFinite()) {
            // 失败 rank 仍需参加下一次全局归约；零向量避免把 NaN 传播给其他 rank。
            output.setZero();
            return false;
        }
        return true;
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
        SolveResult result{
            status,
            iterations,
            initial_residual,
            final_residual,
            final_residual / scale,
        };
        result.performance = current_performance;
        result.performance.linear_solves = 1;
        result.performance.krylov_iterations = static_cast<std::uint64_t>(iterations);
        return result;
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
        double rho = dotGlobal(shadow, residual);
        SolveStatus status = SolveStatus::MaxIterations;
        int iterations = 0;
        for (int iteration = 1; iteration <= iteration_limit; ++iteration) {
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
            // 将本轮真实残差范数和下一轮 rho 合并到同一次 Allreduce。
            // 下一轮不再单独归约 rho，正常 BiCGSTAB 每轮少一个全局同步点。
            const double local_residual_products[2] = {
                residual.squaredNorm(), shadow.dot(residual)};
            double global_residual_products[2]{};
            sumGlobal(local_residual_products, global_residual_products, 2);
            const double residual_norm =
                std::sqrt(std::max(global_residual_products[0], 0.0));
            if (invalid(residual_norm)) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (residual_norm <= target) {
                status = SolveStatus::Converged;
                break;
            }
            previous_rho = rho;
            rho = global_residual_products[1];
        }
        return finish(status, iterations, initial_residual, scale, b, x);
    }

    const Mesh& mesh;
    ParallelContext parallel;
    LinearSolverConfig config;
    KrylovHalo krylov_halo;
    // 系数和局部矩阵均由求解器拥有快照，Equation/Assembly 可安全地在调用后销毁。
    Eigen::SparseMatrix<double> matrix;
    Eigen::SparseMatrix<double> interior_matrix;
    Eigen::SparseMatrix<double> boundary_matrix;
    std::vector<char> boundary_rows;
    std::vector<RemoteCoupling> remote;
    Eigen::IncompleteCholesky<double> incomplete_cholesky;
    Eigen::IncompleteLUT<double> ilut;
    std::vector<int> amg_cell_to_coarse;
    std::vector<int> amg_aggregate;
    std::vector<double> amg_local_coarse;
    std::vector<double> amg_global_coarse;
    Eigen::SparseMatrix<double> amg_coarse_matrix;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> amg_coarse_solver;
    Eigen::VectorXd amg_inverse_diagonal;
    Eigen::VectorXd amg_residual;
    Eigen::VectorXd amg_product;
    Eigen::VectorXd amg_coarse_rhs;
    Eigen::VectorXd amg_global_rhs;
    Eigen::VectorXd amg_coarse_correction;
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
    bool equation_ready = false;
    bool spmv_pattern_ready = false;
    int amg_updates_since_factorization = 0;
    bool amg_ready = false;
    PerformanceCounters pending_performance;
    mutable PerformanceCounters current_performance;
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
    state.current_performance = state.pending_performance;
    state.pending_performance = {};
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
    state.sumGlobal(local_norms, global_norms, 2);
    const double initial_residual = std::sqrt(std::max(global_norms[0], 0.0));
    const double rhs_norm = std::sqrt(std::max(global_norms[1], 0.0));
    const double scale = std::max({initial_residual, rhs_norm, 1e-30});
    const double target = std::max(
        state.config.absolute_tolerance,
        state.config.relative_tolerance * scale);
    if (initial_residual <= target) {
        SolveResult result{
            SolveStatus::Converged, 0, initial_residual,
            initial_residual, initial_residual / scale,
        };
        result.performance = state.current_performance;
        result.performance.linear_solves = 1;
        return result;
    }
    if (!state.factorization_succeeded) {
        SolveResult result{
            SolveStatus::NumericalFailure, 0, initial_residual,
            initial_residual, initial_residual / scale,
        };
        result.performance = state.current_performance;
        result.performance.linear_solves = 1;
        return result;
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
    result.performance = state.current_performance;
    result.performance.linear_solves = 1;
    result.performance.krylov_iterations = static_cast<std::uint64_t>(completed);
    return result;
}

}  // babelsim 命名空间
