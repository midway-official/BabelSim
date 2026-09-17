#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace babelsim {

// 线性系统的配置和结果属于数值控制，不暴露任何稀疏矩阵或后端实现。
enum class LinearSolverType {
    ConjugateGradient,
    BiCGSTAB,
};

enum class PreconditionerType {
    None,
    IncompleteCholesky,
    ILUT,
    AlgebraicMultigrid,
};

enum class SolveStatus {
    Converged,
    MaxIterations,
    NumericalFailure,
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
    std::uint64_t preconditioner_setups = 0;
    std::uint64_t preconditioner_applications = 0;
    // Result I/O is recorded separately from numerical work so one-time and
    // per-write costs are not mistaken for Krylov or assembly hotspots.
    std::uint64_t output_writes = 0;
    double elapsed_seconds = 0.0;
    double assembly_seconds = 0.0;
    double preconditioner_seconds = 0.0;
    double preconditioner_apply_seconds = 0.0;
    double linear_solve_seconds = 0.0;
    double sparse_matvec_seconds = 0.0;
    double halo_seconds = 0.0;
    double global_reduction_seconds = 0.0;
    double output_seconds = 0.0;

    PerformanceCounters& operator+=(const PerformanceCounters& other) {
        linear_solves += other.linear_solves;
        krylov_iterations += other.krylov_iterations;
        sparse_matvecs += other.sparse_matvecs;
        halo_exchanges += other.halo_exchanges;
        halo_bytes += other.halo_bytes;
        global_reductions += other.global_reductions;
        equation_assemblies += other.equation_assemblies;
        preconditioner_setups += other.preconditioner_setups;
        preconditioner_applications += other.preconditioner_applications;
        output_writes += other.output_writes;
        elapsed_seconds += other.elapsed_seconds;
        assembly_seconds += other.assembly_seconds;
        preconditioner_seconds += other.preconditioner_seconds;
        preconditioner_apply_seconds += other.preconditioner_apply_seconds;
        linear_solve_seconds += other.linear_solve_seconds;
        sparse_matvec_seconds += other.sparse_matvec_seconds;
        halo_seconds += other.halo_seconds;
        global_reduction_seconds += other.global_reduction_seconds;
        output_seconds += other.output_seconds;
        return *this;
    }
};

struct LinearSolverConfig {
    LinearSolverType solver = LinearSolverType::BiCGSTAB;
    PreconditionerType preconditioner = PreconditionerType::ILUT;
    double absolute_tolerance = 1e-12;
    double relative_tolerance = 1e-8;
    int max_iterations = 1000;
    bool warm_start = false;
    // ILUT 强度；默认值保持原有后端行为。非默认值用于高纵横比/强对流网格。
    double ilut_drop_tolerance = 1e-3;
    int ilut_fill_factor = 2;
    // AMG 只作为 Krylov 预条件器；这些参数只改变计算后端，不改变方程 API。
    int amg_max_levels = 12;
    int amg_coarse_size = 48;
    int amg_smoothing_steps = 2;
    // AMG 仅作 Krylov 预条件器时，允许复用前几次方程的层级和粗层分解。
    // 线性算子始终使用当前矩阵；1 表示每次更新，保持最保守的数值路径。
    int amg_refresh_interval = 1;

    void validate() const;
};

inline bool operator==(const LinearSolverConfig& a, const LinearSolverConfig& b) {
    return a.solver==b.solver && a.preconditioner==b.preconditioner &&
        a.absolute_tolerance==b.absolute_tolerance && a.relative_tolerance==b.relative_tolerance &&
        a.max_iterations==b.max_iterations && a.warm_start==b.warm_start &&
        a.ilut_drop_tolerance==b.ilut_drop_tolerance && a.ilut_fill_factor==b.ilut_fill_factor &&
        a.amg_max_levels==b.amg_max_levels && a.amg_coarse_size==b.amg_coarse_size &&
        a.amg_smoothing_steps==b.amg_smoothing_steps && a.amg_refresh_interval==b.amg_refresh_interval;
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
