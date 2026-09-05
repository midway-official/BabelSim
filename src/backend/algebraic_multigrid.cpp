#include "backend/algebraic_multigrid.h"

#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace babelsim::detail {
namespace {

using SparseMatrix = Eigen::SparseMatrix<double>;
constexpr double diagonal_tolerance = 1e-30;
constexpr double smoothing_weight = 2.0 / 3.0;

bool finite(const Eigen::VectorXd& values) {
    return values.allFinite();
}

struct Level {
    SparseMatrix matrix;
    SparseMatrix prolongation;
    Eigen::VectorXd inverse_diagonal;
    Eigen::VectorXd residual;
    Eigen::VectorXd right_hand_side;
    Eigen::VectorXd correction;
    std::unique_ptr<Eigen::SparseLU<SparseMatrix>> direct_solver;
};

std::vector<int> aggregates(const SparseMatrix& matrix, int& count) {
    const Eigen::Index rows = matrix.rows();
    std::vector<int> strongest(static_cast<std::size_t>(rows), -1);
    std::vector<double> strength(static_cast<std::size_t>(rows), -1.0);
    for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(matrix, column); entry; ++entry) {
            const Eigen::Index row = entry.row();
            if (row == column) continue;
            const double value = std::abs(entry.value());
            if (value > strength[static_cast<std::size_t>(row)]) {
                strength[static_cast<std::size_t>(row)] = value;
                strongest[static_cast<std::size_t>(row)] = static_cast<int>(column);
            }
            if (value > strength[static_cast<std::size_t>(column)]) {
                strength[static_cast<std::size_t>(column)] = value;
                strongest[static_cast<std::size_t>(column)] = static_cast<int>(row);
            }
        }
    }

    std::vector<int> result(static_cast<std::size_t>(rows), -1);
    count = 0;
    for (Eigen::Index row = 0; row < rows; ++row) {
        if (result[static_cast<std::size_t>(row)] != -1) continue;
        result[static_cast<std::size_t>(row)] = count;
        const int neighbour = strongest[static_cast<std::size_t>(row)];
        if (neighbour >= 0 && neighbour < rows &&
            result[static_cast<std::size_t>(neighbour)] == -1) {
            result[static_cast<std::size_t>(neighbour)] = count;
        }
        ++count;
    }
    return result;
}

SparseMatrix prolongation(const std::vector<int>& aggregate, int coarse_rows) {
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(aggregate.size());
    for (std::size_t row = 0; row < aggregate.size(); ++row) {
        entries.emplace_back(static_cast<Eigen::Index>(row), aggregate[row], 1.0);
    }
    SparseMatrix result(static_cast<Eigen::Index>(aggregate.size()), coarse_rows);
    result.setFromTriplets(entries.begin(), entries.end());
    result.makeCompressed();
    return result;
}

SparseMatrix galerkin(const SparseMatrix& matrix, const SparseMatrix& interpolation) {
    SparseMatrix coarse = interpolation.transpose() * matrix * interpolation;
    coarse.prune(0.0);
    coarse.makeCompressed();
    return coarse;
}

bool initializeDiagonal(Level& level) {
    const Eigen::Index rows = level.matrix.rows();
    level.inverse_diagonal.resize(rows);
    level.residual.resize(rows);
    level.right_hand_side.resize(rows);
    level.correction.resize(rows);
    for (Eigen::Index row = 0; row < rows; ++row) {
        const double diagonal = level.matrix.coeff(row, row);
        if (!std::isfinite(diagonal) || std::abs(diagonal) <= diagonal_tolerance) {
            return false;
        }
        level.inverse_diagonal[row] = 1.0 / diagonal;
    }
    return true;
}

}  // 匿名命名空间

struct AlgebraicMultigrid::Implementation {
    explicit Implementation(LinearSolverConfig value)
        : config(std::move(value))
    {}

    void smooth(Level& level, const Eigen::VectorXd& right_hand_side, Eigen::VectorXd& solution) {
        for (int sweep = 0; sweep < config.amg_smoothing_steps; ++sweep) {
            level.residual.noalias() = right_hand_side - level.matrix * solution;
            solution.noalias() += smoothing_weight *
                level.inverse_diagonal.cwiseProduct(level.residual);
        }
    }

    bool vCycle(std::size_t index, const Eigen::VectorXd& right_hand_side, Eigen::VectorXd& solution) {
        Level& level = levels[index];
        if (index + 1U == levels.size()) {
            solution = level.direct_solver->solve(right_hand_side);
            return level.direct_solver->info() == Eigen::Success && finite(solution);
        }

        smooth(level, right_hand_side, solution);
        level.residual.noalias() = right_hand_side - level.matrix * solution;
        Level& coarse = levels[index + 1U];
        coarse.right_hand_side.noalias() = level.prolongation.transpose() * level.residual;
        coarse.correction.setZero();
        if (!vCycle(index + 1U, coarse.right_hand_side, coarse.correction)) return false;
        solution.noalias() += level.prolongation * coarse.correction;
        smooth(level, right_hand_side, solution);
        return finite(solution);
    }

    bool factorizeCoarsest(bool analyze_pattern) {
        Level& coarse = levels.back();
        if (!coarse.direct_solver) {
            coarse.direct_solver = std::make_unique<Eigen::SparseLU<SparseMatrix>>();
            analyze_pattern = true;
        }
        if (analyze_pattern) coarse.direct_solver->analyzePattern(coarse.matrix);
        coarse.direct_solver->factorize(coarse.matrix);
        return coarse.direct_solver->info() == Eigen::Success;
    }

    bool updateHierarchy(const SparseMatrix& matrix, bool rebuild) {
        if (matrix.rows() != matrix.cols() || matrix.rows() == 0) {
            throw std::invalid_argument("AMG matrix must be non-empty and square");
        }
        if (rebuild) {
            levels.clear();
            levels.push_back({});
            levels.front().matrix = matrix;
            for (int depth = 1; depth < config.amg_max_levels; ++depth) {
                Level& fine = levels.back();
                if (fine.matrix.rows() <= config.amg_coarse_size) break;
                int coarse_rows = 0;
                const std::vector<int> aggregate = aggregates(fine.matrix, coarse_rows);
                if (coarse_rows <= 0 || coarse_rows >= fine.matrix.rows()) break;
                fine.prolongation = prolongation(aggregate, coarse_rows);
                // vector 扩容会使 fine 引用失效，因此先完成 Galerkin 计算再插入层级。
                SparseMatrix coarse = galerkin(fine.matrix, fine.prolongation);
                levels.push_back({});
                levels.back().matrix = std::move(coarse);
            }
        } else {
            levels.front().matrix = matrix;
            for (std::size_t index = 0; index + 1U < levels.size(); ++index) {
                levels[index + 1U].matrix =
                    galerkin(levels[index].matrix, levels[index].prolongation);
            }
        }

        for (Level& level : levels) {
            if (!initializeDiagonal(level)) return false;
        }
        return factorizeCoarsest(rebuild);
    }

    LinearSolverConfig config;
    std::vector<Level> levels;
    bool prepared = false;
};

AlgebraicMultigrid::AlgebraicMultigrid(LinearSolverConfig config)
    : m_implementation(std::make_unique<Implementation>(std::move(config)))
{}

AlgebraicMultigrid::~AlgebraicMultigrid() = default;
AlgebraicMultigrid::AlgebraicMultigrid(AlgebraicMultigrid&&) noexcept = default;
AlgebraicMultigrid& AlgebraicMultigrid::operator=(AlgebraicMultigrid&&) noexcept = default;

void AlgebraicMultigrid::compute(const SparseMatrix& matrix) {
    if (!m_implementation) throw std::logic_error("AMG is moved-from");
    m_implementation->prepared = m_implementation->updateHierarchy(matrix, true);
}

void AlgebraicMultigrid::factorize(const SparseMatrix& matrix) {
    if (!m_implementation || !m_implementation->prepared) {
        throw std::logic_error("AMG pattern must be computed before factorization");
    }
    m_implementation->prepared = m_implementation->updateHierarchy(matrix, false);
}

bool AlgebraicMultigrid::ready() const {
    return m_implementation && m_implementation->prepared;
}

bool AlgebraicMultigrid::apply(const Eigen::VectorXd& input, Eigen::VectorXd& output) {
    if (!ready() || input.size() != m_implementation->levels.front().matrix.rows() ||
        &input == &output) {
        return false;
    }
    output.setZero(input.size());
    return m_implementation->vCycle(0, input, output);
}

SolveResult AlgebraicMultigrid::solve(
    const Eigen::VectorXd& right_hand_side,
    Eigen::VectorXd& solution)
{
    if (!ready() || right_hand_side.size() != m_implementation->levels.front().matrix.rows()) {
        throw std::invalid_argument("AMG linear system is not prepared");
    }
    Implementation& state = *m_implementation;
    if (!state.config.warm_start || solution.size() != right_hand_side.size()) {
        solution.setZero(right_hand_side.size());
    }
    Level& finest = state.levels.front();
    finest.residual.noalias() = right_hand_side - finest.matrix * solution;
    const double initial_residual = finest.residual.norm();
    const double scale = std::max({initial_residual, right_hand_side.norm(), 1e-30});
    const double target = std::max(
        state.config.absolute_tolerance, state.config.relative_tolerance * scale);
    if (!std::isfinite(initial_residual) || !finite(solution)) {
        return {SolveStatus::NumericalFailure, 0, initial_residual, initial_residual,
                initial_residual / scale};
    }
    for (int iteration = 0; iteration <= state.config.max_iterations; ++iteration) {
        const double residual = finest.residual.norm();
        if (!std::isfinite(residual)) {
            return {SolveStatus::NumericalFailure, iteration, initial_residual, residual, residual / scale};
        }
        if (residual <= target) {
            return {SolveStatus::Converged, iteration, initial_residual, residual, residual / scale};
        }
        if (iteration == state.config.max_iterations) break;
        finest.correction.setZero();
        // residual 是 smooth 的输出缓冲，不能同时作为 V-cycle 的右端项。
        finest.right_hand_side = finest.residual;
        if (!state.vCycle(0, finest.right_hand_side, finest.correction)) {
            return {SolveStatus::NumericalFailure, iteration, initial_residual, residual, residual / scale};
        }
        solution += finest.correction;
        finest.residual.noalias() = right_hand_side - finest.matrix * solution;
    }
    const double residual = finest.residual.norm();
    return {SolveStatus::MaxIterations, state.config.max_iterations, initial_residual,
            residual, residual / scale};
}

}  // babelsim::detail 命名空间
