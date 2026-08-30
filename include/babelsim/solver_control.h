#pragma once

namespace babelsim {

// 线性系统的配置和结果属于数值控制，不暴露任何稀疏矩阵或后端实现。
enum class LinearSolverType {
    ConjugateGradient,
    BiCGSTAB,
};

enum class PreconditionerType {
    IncompleteCholesky,
    ILUT,
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

    void validate() const;
};

struct SolveResult {
    SolveStatus status = SolveStatus::NumericalFailure;
    int iterations = 0;
    double initial_residual = 0.0;
    double final_residual = 0.0;
    double relative_residual = 0.0;

    bool converged() const { return status == SolveStatus::Converged; }
};

}  // babelsim 命名空间
