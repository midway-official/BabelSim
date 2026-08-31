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
    void evaluate(fvc::ScalarGradient operation, VectorField& result);
    void evaluate(fvc::NormalGradient operation, ScalarField& result);
    void evaluate(fvc::ScalarDiffusionFlux operation, ScalarField& result);
    void evaluate(fvc::VectorGradient operation, TensorField& result);
    void evaluate(fvc::FaceFlux operation, ScalarField& result);
    void evaluate(fvc::FaceDivergence operation, ScalarField& result);
    void evaluate(fvc::VectorDivergence operation, ScalarField& result);
    void evaluate(fvc::ScalarConvection operation, ScalarField& result);
    void evaluate(fvc::VectorConvection operation, VectorField& result);
    void evaluate(fvc::ScalarInterpolation operation, ScalarField& result);
    void evaluate(fvc::VectorInterpolation operation, VectorField& result);
    void evaluate(fvc::ScalarReconstruction operation, ScalarField& result);
    void evaluate(fvc::VectorReconstruction operation, VectorField& result);
    void evaluate(fvc::ScalarLaplacian operation, ScalarField& result);
    void subtract(
        const ScalarField& coefficient,
        fvc::ScalarGradient operation,
        VectorField& target);
    void add(fvc::FaceFlux operation, ScalarField& target, fvc::FaceRegion region);
    void subtract(fvc::ScalarDiffusionFlux operation, ScalarField& target, fvc::FaceRegion region);
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
