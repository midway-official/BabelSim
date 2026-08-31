#pragma once

#include "babelsim/runtime.h"
#include <utility>

namespace babelsim {
// 不可压缩牛顿流体的常物性模型。若未来引入变黏度模型，仍由该层提供等价物性场，
// 而不是让 SIMPLE 接触底层存储或 Case 解析。
struct FluidProperties {
    double density = 1.0;
    double dynamic_viscosity = 1e-3;

    void validate() const;
};

// SIMPLE 只拥有算法控制；空间方法、时间方法和线性后端属于 RunTime 配置。
struct SimpleControl {
    int max_iterations = 1000;
    int non_orthogonal_corrections = 1;
    double velocity_relaxation = 0.7;
    double pressure_relaxation = 0.3;
    double continuity_tolerance = 1e-8;
    double velocity_tolerance = 1e-7;
    LinearSolverConfig velocity_solver{};
    LinearSolverConfig pressure_solver{
        LinearSolverType::ConjugateGradient,
        PreconditionerType::IncompleteCholesky,
        1e-12,
        1e-8,
        1000,
        false,
    };

    void validate() const;
};

inline RuntimeControl simpleRunTimeControl(
    RuntimeControl result,
    const SimpleControl& control)
{
    result.scalar_solver = control.pressure_solver;
    result.vector_solver = control.velocity_solver;
    return result;
}

inline RuntimeControl simpleRunTimeControl(
    const Methods& methods,
    const SimpleControl& control)
{
    RuntimeControl result;
    result.methods = methods;
    return simpleRunTimeControl(std::move(result), control);
}

// SIMPLE 专用的离散入口。它们把欠松弛、rAU 主对角提取和压力参考点固定在
// 数值层；SimpleSolver 只表达算法步骤，而普通 PDE Solver 仍只使用 solve()。
namespace simple {
void initializeFaceFlux(const VectorField& velocity, ScalarField& face_flux);
std::array<SolveResult, 3> solveMomentumEquation(
    const VectorEquationDefinition& equation,
    double relaxation,
    ScalarField& rAU);
SolveResult solvePressureCorrectionEquation(
    const ScalarEquationDefinition& equation,
    bool fix_reference);
}  // simple 命名空间

struct IncompressibleFields {
    explicit IncompressibleFields(const Mesh& mesh)
        : velocity(mesh, FieldLocation::Cell, "U"),
          pressure(mesh, FieldLocation::Cell, "p"),
          face_flux(mesh, FieldLocation::Face, "phi")
    {}

    VectorField velocity;
    ScalarField pressure;
    ScalarField face_flux;
};

// 通用 RunTime 计算的面通量平衡在不可压缩流中就是连续性度量。保留熟悉的名称，
// 但不把全局归约和数据分区暴露给 SIMPLE。
using ContinuityMetrics = FluxBalance;

struct SimpleIterationResult {
    std::array<SolveResult, 3> velocity;
    SolveResult pressure;
    ContinuityMetrics continuity;
    double relative_velocity_change = 0.0;
    double relative_pressure_correction = 0.0;
    bool healthy = false;
    bool linear_converged = false;
    bool converged = false;
};

}  // babelsim 命名空间
