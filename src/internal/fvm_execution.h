#pragma once

#include "babelsim/methods.h"

#include "babelsim/solver.h"
#include "babelsim/parallel.h"
#include "internal/equation_control.h"

#include <array>
#include <memory>

namespace babelsim::detail {
// FVM 执行后端：解释数学描述、同步局部场、装配与求解，独占离散工作区和时间历史。
// RunTime 仅负责运行域与物理时间推进，不再理解 Term、LDU、Eigen 或矩阵系数。
class FvmExecution {
public:
    FvmExecution(const Mesh&, const Methods&, const LinearSolverConfig&, const LinearSolverConfig&,
                 ParallelContext, double delta_t);
    ~FvmExecution();
    FvmExecution(const FvmExecution&) = delete;
    FvmExecution& operator=(const FvmExecution&) = delete;
    void beginStep(double delta_t);
    SolveResult solve(const ScalarEquationDefinition& equation, EquationControl control);
    std::array<SolveResult, 3> solve(
        const VectorEquationDefinition& equation,
        VectorEquationControl control);
    void evaluate(math::ScalarGradient operation, VectorField& result);
    void evaluate(math::NormalGradient operation, ScalarField& result);
    void evaluate(math::ScalarDiffusionFlux operation, ScalarField& result);
    void evaluate(math::VectorGradient operation, TensorField& result);
    void evaluate(math::FaceFlux operation, ScalarField& result);
    void evaluate(math::FaceDivergence operation, ScalarField& result);
    void evaluate(math::VectorDivergence operation, ScalarField& result);
    void evaluate(math::ScalarConvection operation, ScalarField& result);
    void evaluate(math::VectorConvection operation, VectorField& result);
    void evaluate(math::ScalarInterpolation operation, ScalarField& result);
    void evaluate(math::VectorInterpolation operation, VectorField& result);
    void evaluate(math::ScalarReconstruction operation, ScalarField& result);
    void evaluate(math::VectorReconstruction operation, VectorField& result);
    void evaluate(math::ScalarLaplacian operation, ScalarField& result);
    void subtract(
        const ScalarField& coefficient,
        math::ScalarGradient operation,
        VectorField& target);
    void add(math::FaceFlux operation, ScalarField& target, math::FaceRegion region);
    void subtract(math::ScalarDiffusionFlux operation, ScalarField& target, math::FaceRegion region);
    double relativeChange(const VectorField& current, const VectorField& previous) const;
    double relativeChange(const ScalarField& current, const ScalarField& previous) const;
    double relativeMagnitude(const ScalarField& value, const ScalarField& reference) const;
    FluxBalance fluxBalance(const ScalarField& face_flux) const;
    bool all(bool local_condition) const;


private:
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
};
}  // babelsim::detail 命名空间
