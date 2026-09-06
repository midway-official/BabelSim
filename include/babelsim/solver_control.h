#pragma once

#include <cmath>

namespace babelsim {

// 线性系统的配置和结果属于数值控制，不暴露任何稀疏矩阵或后端实现。
enum class LinearSolverType {
    ConjugateGradient,
    BiCGSTAB,
    GMRES,
    AlgebraicMultigrid,
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

struct LinearSolverConfig {
    LinearSolverType solver = LinearSolverType::BiCGSTAB;
    PreconditionerType preconditioner = PreconditionerType::ILUT;
    double absolute_tolerance = 1e-12;
    double relative_tolerance = 1e-8;
    int max_iterations = 1000;
    bool warm_start = false;
    // GMRES 的每轮 Krylov 子空间维数。较大值通常减少重启，但增加局部正交化和缓存占用。
    int gmres_restart = 30;
    // AMG 使用固定聚合与 V-cycle；这些参数只改变计算后端，不改变方程 API。
    int amg_max_levels = 12;
    int amg_coarse_size = 48;
    int amg_smoothing_steps = 2;
    // AMG 仅作 Krylov 预条件器时，允许复用前几次方程的层级和粗层分解。
    // 线性算子始终使用当前矩阵；1 表示每次更新，保持最保守的数值路径。
    int amg_refresh_interval = 1;

    void validate() const;
};

struct SolveResult {
    // 矢量方程的公开结果：各分量必须全部收敛；迭代数求和，绝对残差合成 L2，
    // 相对残差取最差分量。SIMPLE 内部仍可取得各分量结果用于详细诊断。
    SolveStatus status = SolveStatus::NumericalFailure;
    int iterations = 0;
    double initial_residual = 0.0;
    double final_residual = 0.0;
    double relative_residual = 0.0;

    bool converged() const { return status == SolveStatus::Converged; }
    bool healthy() const {
        return status != SolveStatus::NumericalFailure &&
            std::isfinite(initial_residual) &&
            std::isfinite(final_residual) &&
            std::isfinite(relative_residual);
    }
};

}  // babelsim 命名空间
