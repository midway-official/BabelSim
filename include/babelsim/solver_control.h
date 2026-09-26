#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace babelsim {

// 线性系统的配置和结果属于数值控制，不暴露任何稀疏矩阵或后端实现。
enum class LinearSolverType {
    ConjugateGradient,
    BiCGSTAB,
    GMRES,
    FGMRES,
};

enum class PreconditionerType {
    None,
    IncompleteCholesky,
    Hypre,
    GAMG,
    BlockJacobi,
    ASM,
    Jacobi,
};

enum class SolveStatus {
    Converged,
    MaxIterations,
    NumericalFailure,
};

// Rank-local account of each persistent PETSc system. The values are
// intentionally backend-neutral text/counters so no PETSc handle escapes.
struct EquationPerformanceCounters {
    std::string identity;
    std::string ksp_type;
    std::string pc_type;
    std::string last_status;
    std::uint64_t matrix_pattern_builds = 0;
    std::uint64_t matrix_value_updates = 0;
    std::uint64_t rhs_only_solves = 0;
    std::uint64_t linear_solves = 0;
    std::uint64_t krylov_iterations = 0;
    std::uint64_t preconditioner_setups = 0;
    std::uint64_t true_residual_checks = 0;
    double last_initial_residual = 0.0;
    double last_final_residual = 0.0;
    double last_relative_residual = 0.0;
    double pattern_build_seconds = 0.0;
    double matrix_update_seconds = 0.0;
    double preconditioner_seconds = 0.0;
    double rhs_seconds = 0.0;
    double solve_seconds = 0.0;
    double residual_check_seconds = 0.0;
};

// 计算后端的只读性能快照。它只描述工作量与时间，不暴露矩阵、MPI 或存储实现；
// Runtime 可用它同时报告总求解时间和单位时间步/外迭代成本。
struct PerformanceCounters {
    std::uint64_t linear_solves = 0;
    std::uint64_t krylov_iterations = 0;
    std::uint64_t sparse_matvecs = 0;
    std::uint64_t halo_exchanges = 0;
    std::uint64_t halo_bytes = 0;
    std::uint64_t global_reductions = 0;
    std::uint64_t equation_assemblies = 0;
    std::uint64_t matrix_pattern_builds = 0;
    std::uint64_t matrix_value_updates = 0;
    std::uint64_t rhs_only_solves = 0;
    std::uint64_t true_residual_checks = 0;
    std::uint64_t preconditioner_setups = 0;
    std::uint64_t preconditioner_applications = 0;
    // Result I/O is recorded separately from numerical work so one-time and
    // per-write costs are not mistaken for Krylov or assembly hotspots.
    std::uint64_t output_writes = 0;
    double elapsed_seconds = 0.0;
    double assembly_seconds = 0.0;
    double pattern_build_seconds = 0.0;
    double matrix_update_seconds = 0.0;
    double rhs_seconds = 0.0;
    double residual_check_seconds = 0.0;
    double preconditioner_seconds = 0.0;
    double preconditioner_apply_seconds = 0.0;
    double linear_solve_seconds = 0.0;
    double sparse_matvec_seconds = 0.0;
    double halo_seconds = 0.0;
    double global_reduction_seconds = 0.0;
    double output_seconds = 0.0;
    std::vector<EquationPerformanceCounters> equation_systems;

    PerformanceCounters& operator+=(const PerformanceCounters& other) {
        linear_solves += other.linear_solves;
        krylov_iterations += other.krylov_iterations;
        sparse_matvecs += other.sparse_matvecs;
        halo_exchanges += other.halo_exchanges;
        halo_bytes += other.halo_bytes;
        global_reductions += other.global_reductions;
        equation_assemblies += other.equation_assemblies;
        matrix_pattern_builds += other.matrix_pattern_builds;
        matrix_value_updates += other.matrix_value_updates;
        rhs_only_solves += other.rhs_only_solves;
        true_residual_checks += other.true_residual_checks;
        preconditioner_setups += other.preconditioner_setups;
        preconditioner_applications += other.preconditioner_applications;
        output_writes += other.output_writes;
        elapsed_seconds += other.elapsed_seconds;
        assembly_seconds += other.assembly_seconds;
        pattern_build_seconds += other.pattern_build_seconds;
        matrix_update_seconds += other.matrix_update_seconds;
        rhs_seconds += other.rhs_seconds;
        residual_check_seconds += other.residual_check_seconds;
        preconditioner_seconds += other.preconditioner_seconds;
        preconditioner_apply_seconds += other.preconditioner_apply_seconds;
        linear_solve_seconds += other.linear_solve_seconds;
        sparse_matvec_seconds += other.sparse_matvec_seconds;
        halo_seconds += other.halo_seconds;
        global_reduction_seconds += other.global_reduction_seconds;
        output_seconds += other.output_seconds;
        for (const auto& source : other.equation_systems) {
            auto target = std::find_if(equation_systems.begin(), equation_systems.end(),
                [&](const EquationPerformanceCounters& value) { return value.identity == source.identity; });
            if (target == equation_systems.end()) equation_systems.push_back(source);
            else {
                target->ksp_type = source.ksp_type;
                target->pc_type = source.pc_type;
                target->last_status = source.last_status;
                target->matrix_pattern_builds += source.matrix_pattern_builds;
                target->matrix_value_updates += source.matrix_value_updates;
                target->rhs_only_solves += source.rhs_only_solves;
                target->linear_solves += source.linear_solves;
                target->krylov_iterations += source.krylov_iterations;
                target->preconditioner_setups += source.preconditioner_setups;
                target->true_residual_checks += source.true_residual_checks;
                target->last_initial_residual = source.last_initial_residual;
                target->last_final_residual = source.last_final_residual;
                target->last_relative_residual = source.last_relative_residual;
                target->pattern_build_seconds += source.pattern_build_seconds;
                target->matrix_update_seconds += source.matrix_update_seconds;
                target->preconditioner_seconds += source.preconditioner_seconds;
                target->rhs_seconds += source.rhs_seconds;
                target->solve_seconds += source.solve_seconds;
                target->residual_check_seconds += source.residual_check_seconds;
            }
        }
        return *this;
    }
};

struct LinearSolverConfig {
    LinearSolverType solver = LinearSolverType::BiCGSTAB;
    PreconditionerType preconditioner = PreconditionerType::BlockJacobi;
    double absolute_tolerance = 1e-12;
    double relative_tolerance = 1e-8;
    int max_iterations = 1000;
    bool warm_start = false;
    // GAMG hierarchy controls; Hypre uses its native defaults unless PETSc
    // options are explicitly added to the case-level typed configuration.
    int amg_max_levels = 12;
    int amg_coarse_size = 48;
    int amg_smoothing_steps = 2;

    void validate() const;
};

inline bool operator==(const LinearSolverConfig& a, const LinearSolverConfig& b) {
    return a.solver==b.solver && a.preconditioner==b.preconditioner &&
        a.absolute_tolerance==b.absolute_tolerance && a.relative_tolerance==b.relative_tolerance &&
        a.max_iterations==b.max_iterations && a.warm_start==b.warm_start &&
        a.amg_max_levels==b.amg_max_levels && a.amg_coarse_size==b.amg_coarse_size &&
        a.amg_smoothing_steps==b.amg_smoothing_steps;
}

// 所有后端使用同一合同：原始（未预条件）真残差的全局 L2 范数，
// ||b-Ax|| <= max(atol, rtol * max(||r0||, ||b||))。
// 不引入有量纲的隐藏下限，也不把舍入/停滞当作收敛。
inline double residualScale(double initial, double rhs) {
    const double value = std::max(initial, rhs);
    return value > 0.0 ? value : 1.0;
}

inline double residualTarget(const LinearSolverConfig& config, double scale) {
    return std::max(config.absolute_tolerance, config.relative_tolerance * scale);
}

inline bool residualConverged(double residual, double target) {
    return std::isfinite(residual) && std::isfinite(target) && residual <= target;
}

struct EquationResidual {
    // 原方程（无欠松弛、无参考点惩罚）的全局 L2 残差。
    double norm = 0.0;
    double scale = 0.0; // ||A*x|| + ||b||，对当前非线性系数重新装配
    double relative() const {
        if (!std::isfinite(norm) || !std::isfinite(scale))
            return std::numeric_limits<double>::infinity();
        return scale > 0.0 ? norm / scale : norm;
    }
    bool converged(double absolute_tolerance, double relative_tolerance) const {
        return std::isfinite(scale) && residualConverged(norm,
            std::max(absolute_tolerance, relative_tolerance * scale));
    }
};

struct SolveResult {
    // 矢量方程的公开结果：各分量必须全部收敛；迭代数求和，绝对残差合成 L2，
    // 相对残差取最差分量。SIMPLE 内部仍可取得各分量结果用于详细诊断。
    SolveStatus status = SolveStatus::NumericalFailure;
    int iterations = 0;
    double initial_residual = 0.0;
    double final_residual = 0.0;
    double relative_residual = 0.0;
    PerformanceCounters performance;

    SolveResult() = default;
    SolveResult(
        SolveStatus solve_status,
        int iteration_count,
        double initial,
        double final,
        double relative)
        : status(solve_status),
          iterations(iteration_count),
          initial_residual(initial),
          final_residual(final),
          relative_residual(relative)
    {}

    bool converged() const { return status == SolveStatus::Converged; }
    bool healthy() const {
        return status != SolveStatus::NumericalFailure &&
            std::isfinite(initial_residual) &&
            std::isfinite(final_residual) &&
            std::isfinite(relative_residual);
    }
};


}  // babelsim 命名空间
