#include "internal/mesh_access.h"
#include "babelsim/distributed_solver.h"

#include "babelsim/mpi_support.h"
#include "algebra/inplace_preconditioner.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace babelsim {
namespace {

#ifndef BABELSIM_ASYNC_KRYLOV_HALO
#define BABELSIM_ASYNC_KRYLOV_HALO 1
#endif

#ifndef BABELSIM_CSR_SPMV
#define BABELSIM_CSR_SPMV 1
#endif

#ifndef BABELSIM_INPLACE_PRECONDITIONER
#define BABELSIM_INPLACE_PRECONDITIONER 1
#endif

constexpr double breakdown_tolerance = 0.0;
using Clock = std::chrono::steady_clock;
using SparseMatrix = Eigen::SparseMatrix<double>;

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

// 分布式 Krylov 热路径只需要对局部矩阵做 y=A*x。Eigen 的列压缩矩阵适合
// 预条件器，但对按行的有限体积 stencil 会产生两次稀疏遍历（interior 和
// boundary），且每行访问跨越整列存储。这个私有 CSR 视图复用 Eigen 的系数
// pattern，仅缓存 row/column/source-position；每次方程更新只覆盖数值，不重新
// 分配索引。它不改变公共 Equation/Field API，也不参与预条件器因子化。
struct CsrSpmvBlock {
    struct Entry {
        int column = 0;
        int source_position = 0;
    };

    void build(
        const Eigen::SparseMatrix<double>& source,
        const std::vector<char>& boundary_rows,
        bool boundary)
    {
        if (source.rows() != static_cast<Eigen::Index>(boundary_rows.size())) {
            throw std::invalid_argument("CSR SpMV row mask has the wrong size");
        }
        row_offsets.clear();
        columns.clear();
        source_positions.clear();
        values.clear();
        active_rows.clear();
        inactive_rows.clear();
        std::vector<std::vector<Entry>> rows(static_cast<std::size_t>(source.rows()));
        for (Eigen::Index column = 0; column < source.outerSize(); ++column) {
            Eigen::Index source_position = source.outerIndexPtr()[column];
            for (Eigen::SparseMatrix<double>::InnerIterator entry(source, column);
                 entry; ++entry) {
                if ((boundary_rows[static_cast<std::size_t>(entry.row())] != 0) != boundary)
                    { ++source_position; continue; }
                rows[static_cast<std::size_t>(entry.row())].push_back(
                    {static_cast<int>(entry.col()), static_cast<int>(source_position)});
                ++source_position;
            }
        }
        row_offsets.assign(static_cast<std::size_t>(source.rows()) + 1U, 0);
        for (Eigen::Index row = 0; row < source.rows(); ++row) {
            auto& entries = rows[static_cast<std::size_t>(row)];
            std::sort(entries.begin(), entries.end(),
                [](const Entry& left, const Entry& right) {
                    return left.column < right.column;
                });
            row_offsets[static_cast<std::size_t>(row + 1)] =
                row_offsets[static_cast<std::size_t>(row)] +
                static_cast<Eigen::Index>(entries.size());
            if (!entries.empty()) {
                active_rows.push_back(static_cast<int>(row));
            } else {
                inactive_rows.push_back(static_cast<int>(row));
            }
        }
        columns.resize(static_cast<std::size_t>(row_offsets.back()));
        source_positions.resize(columns.size());
        values.resize(columns.size());
        std::size_t position = 0;
        const double* source_values = source.valuePtr();
        for (const auto& entries : rows) {
            for (const Entry& entry : entries) {
                columns[position] = entry.column;
                source_positions[position] = entry.source_position;
                values[position] = source_values[static_cast<std::size_t>(entry.source_position)];
                ++position;
            }
        }
    }

    void update(const Eigen::SparseMatrix<double>& source) {
        const double* source_values = source.valuePtr();
        for (std::size_t position = 0; position < values.size(); ++position) {
            values[position] = source_values[
                static_cast<std::size_t>(source_positions[position])];
        }
    }

    void multiply(
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output,
        bool add) const
    {
        const double* x = input.data();
        double* y = output.data();
        const int* offsets = row_offsets.data();
        const int* columns_data = columns.data();
        const double* values_data = values.data();
        // The interior block fully overwrites active rows.  Only rows without
        // an interior contribution need clearing before the boundary block is
        // accumulated; this avoids streaming over the whole vector each SpMV.
        if (!add) {
            for (const int row_value : inactive_rows) {
                y[static_cast<std::size_t>(row_value)] = 0.0;
            }
        }
        for (const int row_value : active_rows) {
            const std::size_t row = static_cast<std::size_t>(row_value);
            double sum = 0.0;
            const int begin = offsets[row];
            const int end = offsets[row + 1U];
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#endif
            for (int position = begin; position < end; ++position) {
                sum += values_data[position] * x[columns_data[position]];
            }
            if (add) y[row] += sum;
            else y[row] = sum;
        }
    }

private:
    std::vector<int> row_offsets;
    std::vector<int> columns;
    std::vector<int> source_positions;
    std::vector<int> active_rows;
    std::vector<int> inactive_rows;
    std::vector<double> values;
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
        // 每个 owned 值都会被覆盖；第一层 ghost 也由下面的交换完整写入。
        // 不再为每次 SpMV 清零整个局部 cell 布局，避免与 Krylov 向量长度无关的
        // O(ghost_cells) 内存流量。partitionMesh 保证跨分区耦合只连接第一层 ghost。
        for (Index cell : detail::meshData(*m_mesh).owned_cells) {
            m_values[static_cast<std::size_t>(cell)] = owned_values[detail::ownedIndex(*m_mesh, cell)];
        }
#if BABELSIM_ASYNC_KRYLOV_HALO
        m_exchange.beginFirstLayer(m_values);
#else
        // A/B control: retain the pre-overlap blocking path while keeping the
        // same cached communication layout and numerical operations.
        m_exchange.exchangeFirstLayer(m_values);
#endif
        m_active = true;
    }

    void finish() {
        if (!m_active) throw std::logic_error("Krylov halo exchange is not active");
#if BABELSIM_ASYNC_KRYLOV_HALO
        m_exchange.finishFirstLayer(m_values);
#endif
        m_active = false;
    }

    double value(Index ghost_cell) const {
        return m_values.at(static_cast<std::size_t>(ghost_cell));
    }

    std::size_t bytes() const { return m_exchange.plannedBytes(1, true); }

private:
    const Mesh* m_mesh;
    HaloExchange m_exchange;
    std::vector<double> m_values;
    bool m_active = false;
};

// A coarse vector is distributed by aggregate owner.  Unlike the old AMG
// implementation, a rank stores only the coarse rows it owns and exchanges
// the off-rank columns requested by those rows.  The request graph is built
// once for a fixed matrix pattern; values use the same packed layout on every
// coarse Jacobi sweep.
class DistributedCoarseHalo {
public:
    DistributedCoarseHalo() = default;

    DistributedCoarseHalo(ParallelContext parallel, std::vector<int> owners)
        : m_parallel(parallel), m_owners(std::move(owners)) {}

    void reset(ParallelContext parallel, std::vector<int> owners) {
        m_parallel = parallel;
        m_owners = std::move(owners);
        m_local_index.assign(m_owners.size(), -1);
        m_remote_slot.assign(m_owners.size(), -1);
        m_local_global_ids.clear();
        m_requested_global_ids.clear();
        m_send_indices.clear();
        m_send_buffer.clear();
        m_receive_buffer.clear();
        m_active = false;
    }

    void setLocalIds(const std::vector<int>& ids) {
        m_local_global_ids = ids;
        m_local_index.assign(m_owners.size(), -1);
        for (std::size_t i = 0; i < ids.size(); ++i) {
            if (ids[i] < 0 || static_cast<std::size_t>(ids[i]) >= m_local_index.size()) {
                throw std::invalid_argument("AMG local coarse id is invalid");
            }
            m_local_index[static_cast<std::size_t>(ids[i])] = static_cast<int>(i);
        }
    }

    void build(const std::vector<std::vector<int>>& columns) {
        m_remote_slot.assign(m_owners.size(), -1);
        std::vector<std::vector<int>> requested(static_cast<std::size_t>(m_parallel.size));
        for (const auto& row : columns) {
            for (const int column : row) {
                if (column < 0 || static_cast<std::size_t>(column) >= m_owners.size()) {
                    throw std::invalid_argument("AMG coarse column id is invalid");
                }
                if (m_owners[static_cast<std::size_t>(column)] == m_parallel.rank) continue;
                requested[static_cast<std::size_t>(m_owners[static_cast<std::size_t>(column)])].push_back(column);
            }
        }
        for (auto& ids : requested) {
            std::sort(ids.begin(), ids.end());
            ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        }
        m_receive_counts.assign(static_cast<std::size_t>(m_parallel.size), 0);
        m_receive_offsets.assign(static_cast<std::size_t>(m_parallel.size), 0);
        std::vector<int> outgoing_ids;
        for (int peer = 0; peer < m_parallel.size; ++peer) {
            const std::size_t p = static_cast<std::size_t>(peer);
            m_receive_offsets[p] = detail::mpiCount(outgoing_ids.size(), "coarse halo request offset");
            m_receive_counts[p] = detail::mpiCount(requested[p].size(), "coarse halo request count");
            outgoing_ids.insert(outgoing_ids.end(), requested[p].begin(), requested[p].end());
        }
        m_send_counts.assign(static_cast<std::size_t>(m_parallel.size), 0);
        detail::checkMpi(MPI_Alltoall(
            m_receive_counts.data(), 1, MPI_INT, m_send_counts.data(), 1, MPI_INT,
            m_parallel.communicator), "MPI_Alltoall(AMG coarse halo counts)");
        m_send_offsets.assign(static_cast<std::size_t>(m_parallel.size), 0);
        std::size_t incoming_size = 0;
        for (int peer = 0; peer < m_parallel.size; ++peer) {
            const std::size_t p = static_cast<std::size_t>(peer);
            m_send_offsets[p] = detail::mpiCount(incoming_size, "coarse halo send offset");
            incoming_size += static_cast<std::size_t>(m_send_counts[p]);
        }
        std::vector<int> incoming_ids(incoming_size);
        int dummy = 0;
        detail::checkMpi(MPI_Alltoallv(
            outgoing_ids.empty() ? &dummy : outgoing_ids.data(), m_receive_counts.data(),
            m_receive_offsets.data(), MPI_INT,
            incoming_ids.empty() ? &dummy : incoming_ids.data(), m_send_counts.data(),
            m_send_offsets.data(), MPI_INT, m_parallel.communicator),
            "MPI_Alltoallv(AMG coarse halo requests)");
        m_send_indices.resize(incoming_ids.size());
        for (std::size_t i = 0; i < incoming_ids.size(); ++i) {
            const int id = incoming_ids[i];
            if (id < 0 || static_cast<std::size_t>(id) >= m_local_index.size() ||
                m_local_index[static_cast<std::size_t>(id)] < 0) {
                throw std::runtime_error("AMG coarse halo request is not locally owned");
            }
            m_send_indices[i] = m_local_index[static_cast<std::size_t>(id)];
        }
        m_requested_global_ids = std::move(outgoing_ids);
        for (std::size_t i = 0; i < m_requested_global_ids.size(); ++i) {
            const int id = m_requested_global_ids[i];
            if (id >= 0 && static_cast<std::size_t>(id) < m_remote_slot.size()) {
                m_remote_slot[static_cast<std::size_t>(id)] = static_cast<int>(i);
            }
        }
        m_send_buffer.resize(m_send_indices.size());
        m_receive_buffer.resize(m_requested_global_ids.size());
    }

    void begin(const Eigen::VectorXd& local_values) {
        if (m_parallel.size == 1) return;
        if (m_active || local_values.size() != static_cast<Eigen::Index>(m_local_global_ids.size())) {
            throw std::logic_error("AMG coarse halo begin is invalid");
        }
        for (std::size_t i = 0; i < m_send_indices.size(); ++i) {
            m_send_buffer[i] = local_values[m_send_indices[i]];
        }
        detail::checkMpi(MPI_Ialltoallv(
            m_send_buffer.empty() ? &m_dummy : m_send_buffer.data(),
            sendCountsValues().data(), sendOffsetsValues().data(), MPI_DOUBLE,
            m_receive_buffer.empty() ? &m_dummy : m_receive_buffer.data(),
            receiveCountsValues().data(), receiveOffsetsValues().data(), MPI_DOUBLE,
            m_parallel.communicator, &m_request), "MPI_Ialltoallv(AMG coarse values)");
        m_active = true;
    }

    void finish() {
        if (m_parallel.size == 1) return;
        if (!m_active) throw std::logic_error("AMG coarse halo finish without begin");
        detail::checkMpi(MPI_Wait(&m_request, MPI_STATUS_IGNORE), "MPI_Wait(AMG coarse values)");
        m_active = false;
    }

    int localSlot(int global_id) const {
        if (global_id < 0 || static_cast<std::size_t>(global_id) >= m_local_index.size()) return -1;
        return m_local_index[static_cast<std::size_t>(global_id)];
    }

    int remoteSlot(int global_id) const {
        if (global_id < 0 || static_cast<std::size_t>(global_id) >= m_remote_slot.size()) return -1;
        return m_remote_slot[static_cast<std::size_t>(global_id)];
    }

    double remoteValue(int slot) const { return m_receive_buffer.at(static_cast<std::size_t>(slot)); }
    const std::vector<int>& localGlobalIds() const { return m_local_global_ids; }

private:
    const std::vector<int>& sendCountsValues() const {
        return m_send_counts;
    }
    const std::vector<int>& sendOffsetsValues() const {
        return m_send_offsets;
    }
    const std::vector<int>& receiveCountsValues() const {
        return m_receive_counts;
    }
    const std::vector<int>& receiveOffsetsValues() const {
        return m_receive_offsets;
    }

    ParallelContext m_parallel;
    std::vector<int> m_owners;
    std::vector<int> m_local_global_ids;
    std::vector<int> m_local_index;
    std::vector<int> m_remote_slot;
    std::vector<int> m_requested_global_ids;
    std::vector<int> m_send_indices;
    std::vector<int> m_send_counts;
    std::vector<int> m_send_offsets;
    std::vector<int> m_receive_counts;
    std::vector<int> m_receive_offsets;
    std::vector<double> m_send_buffer;
    std::vector<double> m_receive_buffer;
    MPI_Request m_request = MPI_REQUEST_NULL;
    double m_dummy = 0.0;
    bool m_active = false;
};

struct CoarseContribution {
    int row = 0;
    int column = 0;
    double value = 0.0;
};

struct CoarseRowEntry {
    int column = 0;
    double value = 0.0;
    int local_column = -1;
    int remote_slot = -1;
};

static_assert(std::is_trivially_copyable<CoarseContribution>::value,
              "AMG contribution must be MPI-byte-copyable");

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
        std::sort(remote.begin(), remote.end(), [](const RemoteCoupling& left,
                                                   const RemoteCoupling& right) {
            if (left.row != right.row) return left.row < right.row;
            return left.ghost_cell < right.ghost_cell;
        });
        remote_rows.reserve(remote.size());
        remote_ghost_cells.reserve(remote.size());
        remote_coefficients.assign(remote.size(), 0.0);
        for (const RemoteCoupling& coupling : remote) {
            remote_rows.push_back(static_cast<int>(coupling.row));
            remote_ghost_cells.push_back(coupling.ghost_cell);
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
        for (std::size_t index = 0; index < remote.size(); ++index) {
            remote_coefficients[index] = remote[index].coefficient;
        }
        equation_ready = true;
    }

    void setMatrix(const Eigen::SparseMatrix<double>& value) {
        if (value.rows() != detail::ownedCellCount(mesh) ||
            value.cols() != detail::ownedCellCount(mesh)) {
            throw std::invalid_argument("distributed local matrix size is invalid");
        }
        if (!spmv_pattern_ready) {
            // The first matrix establishes the immutable sparse pattern.  All
            // later SIMPLE iterations reuse it and only replace contiguous
            // coefficient values below; this avoids an Eigen sparse
            // assignment and allocator traffic on every equation update.
            matrix = value;
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
#if BABELSIM_CSR_SPMV
            interior_csr.build(matrix, boundary_rows, false);
            boundary_csr.build(matrix, boundary_rows, true);
#endif
            spmv_pattern_ready = true;
        } else {
            // SparseAssembly owns a fixed stencil.  Refuse a changed pattern
            // instead of silently pairing a coefficient with the wrong entry.
            if (value.nonZeros() != matrix.nonZeros() ||
                !std::equal(
                    value.outerIndexPtr(),
                    value.outerIndexPtr() + value.outerSize() + 1,
                    matrix.outerIndexPtr()) ||
                !std::equal(
                    value.innerIndexPtr(),
                    value.innerIndexPtr() + value.nonZeros(),
                    matrix.innerIndexPtr())) {
                throw std::logic_error("distributed sparse matrix pattern changed");
            }
            std::copy_n(value.valuePtr(), value.nonZeros(), matrix.valuePtr());
#if !BABELSIM_CSR_SPMV
            // The split Eigen matrices are only needed by the fallback SpMV
            // path.  The default CSR path reads the shared coefficient array.
            for (Eigen::Index column = 0; column < value.outerSize(); ++column) {
                for (Eigen::SparseMatrix<double>::InnerIterator entry(value, column);
                     entry; ++entry) {
                    Eigen::SparseMatrix<double>& target =
                        boundary_rows[static_cast<std::size_t>(entry.row())]
                        ? boundary_matrix : interior_matrix;
                    target.coeffRef(entry.row(), entry.col()) = entry.value();
                }
            }
#endif
#if BABELSIM_CSR_SPMV
            interior_csr.update(matrix);
            boundary_csr.update(matrix);
#endif
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
            ilut.setDroptol(config.ilut_drop_tolerance);
            ilut.setFillfactor(config.ilut_fill_factor);
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
        // Publish variable-length, de-duplicated graph rows.  The old fixed
        // six-neighbour padding made a polyhedral cell silently lose edges.
        // Rows are gathered once during AMG setup; Krylov iterations never see
        // this representation.
        std::vector<int> packed;
        for (Index cell : detail::meshData(mesh).owned_cells) {
            const Index global_cell = detail::globalCellId(mesh, cell);
            std::vector<int> row;
            for (Index neighbour : mesh.cellNeighbours(cell)) {
                if (neighbour != invalid_index) row.push_back(detail::globalCellId(mesh, neighbour));
            }
            std::sort(row.begin(), row.end());
            row.erase(std::unique(row.begin(), row.end()), row.end());
            packed.push_back(global_cell);
            packed.push_back(static_cast<int>(row.size()));
            packed.insert(packed.end(), row.begin(), row.end());
        }
        int local_count = static_cast<int>(packed.size());
        std::vector<int> counts(static_cast<std::size_t>(parallel.size));
        detail::checkMpi(MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                                       parallel.communicator), "MPI_Allgather(AMG graph row sizes)");
        std::vector<int> displacements(counts.size(), 0);
        for (std::size_t rank = 1; rank < counts.size(); ++rank) {
            displacements[rank] = displacements[rank - 1] + counts[rank - 1];
        }
        const int total_count = displacements.empty() ? 0 :
            displacements.back() + counts.back();
        std::vector<int> gathered(static_cast<std::size_t>(total_count));
        const Clock::time_point adjacency_start = Clock::now();
        detail::checkMpi(MPI_Allgatherv(
            packed.data(), local_count, MPI_INT, gathered.data(), counts.data(),
            displacements.data(), MPI_INT, parallel.communicator),
            "MPI_Allgatherv(AMG graph rows)");
        ++pending_performance.global_reductions;
        pending_performance.global_reduction_seconds += secondsSince(adjacency_start);

        std::vector<std::vector<int>> global_adjacency(static_cast<std::size_t>(global_cells));
        for (std::size_t cursor = 0; cursor < gathered.size();) {
            if (cursor + 2 > gathered.size()) throw std::runtime_error("truncated AMG graph row");
            const int row = gathered[cursor++];
            const int degree = gathered[cursor++];
            if (row < 0 || row >= global_cells || degree < 0 ||
                cursor + static_cast<std::size_t>(degree) > gathered.size()) {
                throw std::runtime_error("invalid AMG graph row");
            }
            auto& links = global_adjacency[static_cast<std::size_t>(row)];
            links.insert(links.end(), gathered.begin() + static_cast<std::ptrdiff_t>(cursor),
                         gathered.begin() + static_cast<std::ptrdiff_t>(cursor + degree));
            cursor += static_cast<std::size_t>(degree);
        }
        for (auto& links : global_adjacency) {
            std::sort(links.begin(), links.end());
            links.erase(std::unique(links.begin(), links.end()), links.end());
        }

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
                for (const int neighbour : global_adjacency[static_cast<std::size_t>(cell)]) {
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
        // 只复制 aggregate 的 owner 元数据；粗矩阵和粗向量随后按 owner 分布。
        // 这个上限保护的是全局映射本身，而不是 coarse_count² 的密集矩阵。
        constexpr std::int64_t maximum_coarse_metadata = 10'000'000;
        if (coarse_count <= 0 || coarse_count > maximum_coarse_metadata) {
            throw std::runtime_error("distributed AMG coarse space is invalid");
        }
        // Only a genuinely tiny coarsest problem uses the existing local
        // ILUT fallback.  Normal cases (including the configured coarse size
        // 48) stay on the owner-distributed coarse operator below.
        amg_tiny_local_fallback = coarse_count <= 64 && matrix.rows() <= 64;
        amg_aggregate.resize(static_cast<std::size_t>(matrix.rows()));
        const auto& owned = detail::meshData(mesh).owned_cells;
        for (std::size_t row = 0; row < owned.size(); ++row) {
            amg_aggregate[row] = coarseIndex(detail::globalCellId(mesh, owned[row]));
        }
        amg_coarse_owner.assign(static_cast<std::size_t>(coarse_count), parallel.size);
        for (const int aggregate : amg_aggregate) {
            amg_coarse_owner[static_cast<std::size_t>(aggregate)] =
                std::min(amg_coarse_owner[static_cast<std::size_t>(aggregate)], parallel.rank);
        }
        const Clock::time_point owner_start = Clock::now();
        detail::checkMpi(MPI_Allreduce(
            MPI_IN_PLACE, amg_coarse_owner.data(), static_cast<int>(coarse_count), MPI_INT,
            MPI_MIN, parallel.communicator), "MPI_Allreduce(AMG aggregate owners)");
        ++pending_performance.global_reductions;
        pending_performance.global_reduction_seconds += secondsSince(owner_start);
        amg_coarse_local_ids.clear();
        for (int aggregate = 0; aggregate < static_cast<int>(coarse_count); ++aggregate) {
            if (amg_coarse_owner[static_cast<std::size_t>(aggregate)] == parallel.rank) {
                amg_coarse_local_ids.push_back(aggregate);
            }
        }
        amg_coarse_local_index.assign(static_cast<std::size_t>(coarse_count), -1);
        for (std::size_t local = 0; local < amg_coarse_local_ids.size(); ++local) {
            amg_coarse_local_index[static_cast<std::size_t>(amg_coarse_local_ids[local])] =
                static_cast<int>(local);
        }
        amg_coarse_halo.reset(parallel, amg_coarse_owner);
        amg_coarse_halo.setLocalIds(amg_coarse_local_ids);
        amg_coarse_rhs.resize(static_cast<Eigen::Index>(amg_coarse_local_ids.size()));
        amg_coarse_correction.resize(static_cast<Eigen::Index>(amg_coarse_local_ids.size()));
        amg_coarse_product.resize(static_cast<Eigen::Index>(amg_coarse_local_ids.size()));
    }

    std::vector<CoarseContribution> exchangeCoarseContributions(
        const std::vector<std::vector<CoarseContribution>>& outgoing) const
    {
        if (outgoing.size() != static_cast<std::size_t>(parallel.size)) {
            throw std::invalid_argument("AMG contribution peer count is invalid");
        }
        std::vector<int> send_counts(static_cast<std::size_t>(parallel.size), 0);
        std::vector<int> send_offsets(static_cast<std::size_t>(parallel.size), 0);
        std::vector<CoarseContribution> packed;
        for (int peer = 0; peer < parallel.size; ++peer) {
            const std::size_t p = static_cast<std::size_t>(peer);
            send_offsets[p] = detail::mpiCount(
                packed.size() * sizeof(CoarseContribution), "AMG contribution send offset");
            send_counts[p] = detail::mpiCount(
                outgoing[p].size() * sizeof(CoarseContribution), "AMG contribution send count");
            packed.insert(packed.end(), outgoing[p].begin(), outgoing[p].end());
        }
        std::vector<int> receive_counts(static_cast<std::size_t>(parallel.size), 0);
        detail::checkMpi(MPI_Alltoall(
            send_counts.data(), 1, MPI_INT, receive_counts.data(), 1, MPI_INT,
            parallel.communicator), "MPI_Alltoall(AMG contribution counts)");
        std::vector<int> receive_offsets(static_cast<std::size_t>(parallel.size), 0);
        std::size_t receive_bytes = 0;
        for (int peer = 0; peer < parallel.size; ++peer) {
            const std::size_t p = static_cast<std::size_t>(peer);
            receive_offsets[p] = detail::mpiCount(receive_bytes, "AMG contribution receive offset");
            receive_bytes += static_cast<std::size_t>(receive_counts[p]);
        }
        if (receive_bytes % sizeof(CoarseContribution) != 0) {
            throw std::runtime_error("AMG contribution byte count is not aligned");
        }
        std::vector<CoarseContribution> received(receive_bytes / sizeof(CoarseContribution));
        CoarseContribution dummy{};
        detail::checkMpi(MPI_Alltoallv(
            packed.empty() ? &dummy : packed.data(), send_counts.data(), send_offsets.data(), MPI_BYTE,
            received.empty() ? &dummy : received.data(), receive_counts.data(), receive_offsets.data(),
            MPI_BYTE, parallel.communicator), "MPI_Alltoallv(AMG contributions)");
        return received;
    }

    bool updateDistributedAmg(bool rebuild_pattern) {
        if (rebuild_pattern) buildCoarseMapping();
        amg_inverse_diagonal.resize(matrix.rows());
        bool diagonal_ok = true;
        amg_diagonal_only = true;
        for (const RemoteCoupling& coupling : remote) {
            if (coupling.coefficient != 0.0) amg_diagonal_only = false;
        }
        for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
            const double diagonal = matrix.coeff(row, row);
            if (!std::isfinite(diagonal) || std::abs(diagonal) <= breakdown_tolerance) {
                diagonal_ok = false;
                amg_inverse_diagonal[row] = 0.0;
            } else {
                amg_inverse_diagonal[row] = 1.0 / diagonal;
            }
        }
        for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
            for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix, column);
                 entry; ++entry) {
                if (entry.row() != entry.col() && entry.value() != 0.0) {
                    amg_diagonal_only = false;
                }
            }
        }
        const Clock::time_point diagonal_reduction_start = Clock::now();
        const bool global_diagonal_ok = parallel.maximum(diagonal_ok ? 0 : 1) == 0;
        ++pending_performance.global_reductions;
        pending_performance.global_reduction_seconds +=
            secondsSince(diagonal_reduction_start);
        if (!global_diagonal_ok) {
            amg_ready = false;
            return false;
        }

        std::vector<std::vector<CoarseContribution>> outgoing(
            static_cast<std::size_t>(parallel.size));
        const auto add = [&](int row, int column, double value) {
            if (!std::isfinite(value)) return;
            const int owner = amg_coarse_owner.at(static_cast<std::size_t>(row));
            outgoing[static_cast<std::size_t>(owner)].push_back({row, column, value});
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
        const std::vector<CoarseContribution> received =
            exchangeCoarseContributions(outgoing);
        std::vector<std::map<int, double>> rows(amg_coarse_local_ids.size());
        for (const CoarseContribution& contribution : received) {
            if (contribution.row < 0 ||
                static_cast<std::size_t>(contribution.row) >= amg_coarse_owner.size() ||
                amg_coarse_owner[static_cast<std::size_t>(contribution.row)] != parallel.rank) {
                throw std::runtime_error("AMG contribution arrived at the wrong owner");
            }
            const int local_row = amg_coarse_local_index.at(
                static_cast<std::size_t>(contribution.row));
            rows[static_cast<std::size_t>(local_row)][contribution.column] += contribution.value;
        }
        if (rebuild_pattern) {
            amg_coarse_rows.resize(rows.size());
            std::vector<std::vector<int>> columns(rows.size());
            for (std::size_t row = 0; row < rows.size(); ++row) {
                auto& target = amg_coarse_rows[row];
                target.clear();
                target.reserve(rows[row].size());
                columns[row].reserve(rows[row].size());
                for (const auto& entry : rows[row]) {
                    target.push_back({entry.first, entry.second, -1, -1});
                    columns[row].push_back(entry.first);
                }
            }
            amg_coarse_halo.build(columns);
            amg_coarse_prolongation_halo.reset(parallel, amg_coarse_owner);
            amg_coarse_prolongation_halo.setLocalIds(amg_coarse_local_ids);
            std::vector<std::vector<int>> prolongation_requests(1);
            prolongation_requests.front().reserve(amg_aggregate.size());
            for (const int aggregate : amg_aggregate) {
                prolongation_requests.front().push_back(aggregate);
            }
            amg_coarse_prolongation_halo.build(prolongation_requests);
            for (auto& target : amg_coarse_rows) {
                for (auto& entry : target) {
                    entry.local_column = amg_coarse_halo.localSlot(entry.column);
                    entry.remote_slot = amg_coarse_halo.remoteSlot(entry.column);
                    if (entry.local_column < 0 && entry.remote_slot < 0) {
                        throw std::runtime_error("AMG coarse matrix column has no owner");
                    }
                }
            }
        } else {
            if (rows.size() != amg_coarse_rows.size()) {
                throw std::logic_error("distributed AMG coarse row pattern changed");
            }
            for (std::size_t row = 0; row < rows.size(); ++row) {
                if (rows[row].size() != amg_coarse_rows[row].size()) {
                    throw std::logic_error("distributed AMG coarse column pattern changed");
                }
                for (auto& entry : amg_coarse_rows[row]) {
                    const auto value = rows[row].find(entry.column);
                    if (value == rows[row].end()) {
                        throw std::logic_error("distributed AMG coarse column pattern changed");
                    }
                    entry.value = value->second;
                }
            }
        }
        amg_coarse_inverse_diagonal.assign(amg_coarse_rows.size(), 0.0);
        bool coarse_diagonal_ok = true;
        for (std::size_t row = 0; row < amg_coarse_rows.size(); ++row) {
            const int global_row = amg_coarse_local_ids[row];
            for (const CoarseRowEntry& entry : amg_coarse_rows[row]) {
                if (entry.column == global_row) {
                    if (!std::isfinite(entry.value) ||
                        std::abs(entry.value) <= breakdown_tolerance) {
                        coarse_diagonal_ok = false;
                    } else {
                        amg_coarse_inverse_diagonal[row] = 1.0 / entry.value;
                    }
                    break;
                }
            }
            if (amg_coarse_inverse_diagonal[row] == 0.0) coarse_diagonal_ok = false;
        }
        if (amg_tiny_local_fallback) {
            amg_tiny_ilut.setDroptol(config.ilut_drop_tolerance);
            amg_tiny_ilut.setFillfactor(config.ilut_fill_factor);
            if (rebuild_pattern) amg_tiny_ilut.compute(matrix);
            else amg_tiny_ilut.factorize(matrix);
        }
        amg_residual.resize(matrix.rows());
        amg_product.resize(matrix.rows());
        amg_ready = coarse_diagonal_ok &&
            (!amg_tiny_local_fallback || amg_tiny_ilut.info() == Eigen::Success);
        return amg_ready;
    }

    void multiplyCoarse(
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output,
        DistributedCoarseHalo& halo) {
        halo.begin(input);
        output.setZero();
        // Rows with only local aggregate columns can be evaluated while the
        // remote column values are in flight.
        for (std::size_t row = 0; row < amg_coarse_rows.size(); ++row) {
            double sum = 0.0;
            for (const CoarseRowEntry& entry : amg_coarse_rows[row]) {
                if (entry.local_column >= 0) sum += entry.value * input[entry.local_column];
            }
            output[static_cast<Eigen::Index>(row)] = sum;
        }
        halo.finish();
        for (std::size_t row = 0; row < amg_coarse_rows.size(); ++row) {
            for (const CoarseRowEntry& entry : amg_coarse_rows[row]) {
                if (entry.remote_slot >= 0) {
                    output[static_cast<Eigen::Index>(row)] +=
                        entry.value * halo.remoteValue(entry.remote_slot);
                }
            }
        }
    }

    bool applyDistributedAmg(
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output)
    {
        if (!amg_ready || &input == &output) return false;
        if (amg_tiny_local_fallback) {
#if BABELSIM_INPLACE_PRECONDITIONER
            amg_tiny_ilut.solveInPlace(input, output);
#else
            output = amg_tiny_ilut.solve(input);
#endif
            return amg_tiny_ilut.info() == Eigen::Success && output.allFinite();
        }
        if (amg_diagonal_only) {
            output.array() = amg_inverse_diagonal.array() * input.array();
            return output.allFinite();
        }
        constexpr double weight = 2.0 / 3.0;
        output.array() = weight * amg_inverse_diagonal.array() * input.array();
        for (int sweep = 1; sweep < config.amg_smoothing_steps; ++sweep) {
            apply(output, amg_product);
            amg_residual.noalias() = input;
            amg_residual.noalias() -= amg_product;
            output.array() += weight * amg_inverse_diagonal.array() * amg_residual.array();
        }
        apply(output, amg_product);
        amg_residual.noalias() = input;
        amg_residual.noalias() -= amg_product;

        std::vector<std::vector<CoarseContribution>> outgoing(
            static_cast<std::size_t>(parallel.size));
        amg_coarse_rhs.setZero();
        for (Eigen::Index row = 0; row < amg_residual.size(); ++row) {
            const int aggregate = amg_aggregate[static_cast<std::size_t>(row)];
            const double value = amg_residual[row];
            const int owner = amg_coarse_owner[static_cast<std::size_t>(aggregate)];
            if (owner == parallel.rank) {
                amg_coarse_rhs[amg_coarse_local_index[static_cast<std::size_t>(aggregate)]] += value;
            } else {
                outgoing[static_cast<std::size_t>(owner)].push_back({aggregate, 0, value});
            }
        }
        const std::vector<CoarseContribution> received =
            exchangeCoarseContributions(outgoing);
        for (const CoarseContribution& contribution : received) {
            if (contribution.row < 0 ||
                static_cast<std::size_t>(contribution.row) >= amg_coarse_owner.size() ||
                amg_coarse_owner[static_cast<std::size_t>(contribution.row)] != parallel.rank) {
                throw std::runtime_error("AMG coarse rhs arrived at the wrong owner");
            }
            amg_coarse_rhs[amg_coarse_local_index[static_cast<std::size_t>(contribution.row)]] +=
                contribution.value;
        }

        amg_coarse_correction.setZero();
        amg_coarse_product.resize(amg_coarse_rhs.size());
        // Fixed-point coarse iterations avoid a second global Krylov solve and
        // never replicate the coarse matrix/vector.  The matrix is deliberately
        // kept sparse and owner distributed; halo traffic is only for columns
        // touched by the local coarse rows.
        constexpr int coarse_iterations = 24;
        constexpr double coarse_weight = 0.72;
        for (int sweep = 0; sweep < coarse_iterations; ++sweep) {
            multiplyCoarse(amg_coarse_correction, amg_coarse_product, amg_coarse_halo);
            for (Eigen::Index row = 0; row < amg_coarse_rhs.size(); ++row) {
                amg_coarse_correction[row] += coarse_weight *
                    amg_coarse_inverse_diagonal[static_cast<std::size_t>(row)] *
                    (amg_coarse_rhs[row] - amg_coarse_product[row]);
            }
        }
        if (!amg_coarse_correction.allFinite()) return false;
        amg_coarse_prolongation_halo.begin(amg_coarse_correction);
        amg_coarse_prolongation_halo.finish();
        for (Eigen::Index row = 0; row < output.size(); ++row) {
            const int aggregate = amg_aggregate[static_cast<std::size_t>(row)];
            const int local = amg_coarse_halo.localSlot(aggregate);
            output[row] += local >= 0
                ? amg_coarse_correction[local]
                : amg_coarse_prolongation_halo.remoteValue(
                    amg_coarse_prolongation_halo.remoteSlot(aggregate));
        }
        for (int sweep = 0; sweep < config.amg_smoothing_steps; ++sweep) {
            apply(output, amg_product);
            amg_residual.noalias() = input;
            amg_residual.noalias() -= amg_product;
            output.array() += weight * amg_inverse_diagonal.array() * amg_residual.array();
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
        double local_third,
        bool local_success,
        double& global_first,
        double& global_second,
        double& global_third,
        bool& global_success) const
    {
        const double local[4] = {
            local_success ? local_first : 0.0,
            local_success ? local_second : 0.0,
            local_success ? local_third : 0.0,
            local_success ? 0.0 : 1.0,
        };
        double global[4]{};
        sumGlobal(local, global, 4);
        global_first = global[0];
        global_second = global[1];
        global_third = global[2];
        global_success = global[3] == 0.0;
    }

    void apply(const Eigen::VectorXd& input, Eigen::VectorXd& output) {
        const Clock::time_point start = Clock::now();
        ++current_performance.sparse_matvecs;
        const Clock::time_point halo_start = Clock::now();
        krylov_halo.begin(input);
        current_performance.halo_seconds += secondsSince(halo_start);
        current_performance.halo_bytes += krylov_halo.bytes();
#if BABELSIM_CSR_SPMV
        interior_csr.multiply(input, output, false);
#else
        output.noalias() = interior_matrix * input;
#endif
        const Clock::time_point halo_wait_start = Clock::now();
        krylov_halo.finish();
        ++current_performance.halo_exchanges;
        current_performance.halo_seconds += secondsSince(halo_wait_start);
#if BABELSIM_CSR_SPMV
        boundary_csr.multiply(input, output, true);
#else
        output.noalias() += boundary_matrix * input;
#endif
        for (std::size_t index = 0; index < remote_rows.size(); ++index) {
            output[remote_rows[index]] += remote_coefficients[index] *
                krylov_halo.value(remote_ghost_cells[index]);
        }
        current_performance.sparse_matvec_seconds += secondsSince(start);
    }

    bool precondition(
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output)
    {
        const Clock::time_point start = Clock::now();
        bool local_success = false;
        // A BiCGSTAB exact/near-exact first step can produce a zero s vector.
        // Treat that zero input as an identity preconditioner so the fused
        // products reduction can perform the exact early exit without asking
        // ILUT/AMG to factor a zero RHS.
        if (input.squaredNorm() == 0.0) {
            output.setZero();
            local_success = true;
        } else if (!hasPreconditioner(config)) {
            output = input;
            local_success = true;
        } else if (usesAmg(config)) {
            local_success = applyDistributedAmg(input, output);
        } else if (config.solver == LinearSolverType::ConjugateGradient) {
#if BABELSIM_INPLACE_PRECONDITIONER
            incomplete_cholesky.solveInPlace(input, output);
            local_success = incomplete_cholesky.info() == Eigen::Success;
#else
            output = incomplete_cholesky.solve(input);
            local_success = incomplete_cholesky.info() == Eigen::Success;
#endif
        } else {
#if BABELSIM_INPLACE_PRECONDITIONER
            ilut.solveInPlace(input, output);
            local_success = ilut.info() == Eigen::Success;
#else
            output = ilut.solve(input);
            local_success = ilut.info() == Eigen::Success;
#endif
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
        const double target = residualTarget(config, scale);
        if (!std::isfinite(final_residual)) {
            status = SolveStatus::NumericalFailure;
        } else if (residualConverged(final_residual, target)) {
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
            intermediate.noalias() = residual - alpha * direction_product;
            iterations = iteration;
            // The ||s|| check is carried into the products reduction below,
            // together with t·s and t·t.  This keeps the exact early-exit
            // criterion without adding a fourth global synchronization.
            const bool local_intermediate_precondition_success =
                precondition(intermediate, preconditioned_intermediate);
            apply(preconditioned_intermediate, intermediate_product);
            const double local_products[2] = {
                intermediate_product.dot(intermediate),
                intermediate_product.squaredNorm(),
            };
            const double local_intermediate_squared_norm = intermediate.squaredNorm();
            double global_products[3]{};
            bool global_intermediate_precondition_success = false;
            productsGlobalWithStatus(
                local_products[0], local_products[1],
                local_intermediate_squared_norm,
                local_intermediate_precondition_success,
                global_products[0], global_products[1],
                global_products[2],
                global_intermediate_precondition_success);
            if (!global_intermediate_precondition_success) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            if (invalid(global_products[0]) || invalid(global_products[1])) {
                status = SolveStatus::NumericalFailure;
                break;
            }
            const double intermediate_norm =
                std::sqrt(std::max(global_products[2], 0.0));
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
            residual.noalias() = intermediate - omega * intermediate_product;
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
#if BABELSIM_CSR_SPMV
    CsrSpmvBlock interior_csr;
    CsrSpmvBlock boundary_csr;
#endif
    std::vector<char> boundary_rows;
    std::vector<RemoteCoupling> remote;
    std::vector<int> remote_rows;
    std::vector<Index> remote_ghost_cells;
    std::vector<double> remote_coefficients;
#if BABELSIM_INPLACE_PRECONDITIONER
    detail::InPlaceIncompleteCholesky incomplete_cholesky;
    detail::InPlaceIncompleteLut ilut;
    detail::InPlaceIncompleteLut amg_tiny_ilut;
#else
    Eigen::IncompleteCholesky<double> incomplete_cholesky;
    Eigen::IncompleteLUT<double> ilut;
    Eigen::IncompleteLUT<double> amg_tiny_ilut;
#endif
    std::vector<int> amg_cell_to_coarse;
    std::vector<int> amg_aggregate;
    std::vector<int> amg_coarse_owner;
    std::vector<int> amg_coarse_local_ids;
    std::vector<int> amg_coarse_local_index;
    std::vector<std::vector<CoarseRowEntry>> amg_coarse_rows;
    std::vector<double> amg_coarse_inverse_diagonal;
    DistributedCoarseHalo amg_coarse_halo;
    DistributedCoarseHalo amg_coarse_prolongation_halo;
    Eigen::VectorXd amg_inverse_diagonal;
    Eigen::VectorXd amg_residual;
    Eigen::VectorXd amg_product;
    Eigen::VectorXd amg_coarse_rhs;
    Eigen::VectorXd amg_coarse_correction;
    Eigen::VectorXd amg_coarse_product;
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
    bool amg_diagonal_only = false;
    bool amg_tiny_local_fallback = false;
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
    const double scale = residualScale(initial_residual, rhs_norm);
    const double target = residualTarget(state.config, scale);
    if (residualConverged(initial_residual, target)) {
        SolveResult result{
            SolveStatus::Converged, 0, initial_residual,
            initial_residual, initial_residual / scale,
        };
        result.performance = state.current_performance;
        result.performance.linear_solves = 1;
        return result;
    }
    if (!std::isfinite(initial_residual) || !state.factorization_succeeded) {
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
