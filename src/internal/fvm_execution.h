#pragma once

#include "babelsim/math.h"
#include "babelsim/methods.h"
#include "babelsim/solver.h"

#include <memory>

namespace babelsim::detail {
class ComputeBackend;

// FVM 数值执行层：解释数学描述并调用通用算子；计算后端负责同步、归约与装配。
// 该类独占离散工作区，不依赖 MPI、Eigen 或具体稀疏矩阵实现；方程装配与求解由
// equ:: 显式对象承担，执行层只提供后端、网格与显式算法求值。
class FvmExecution {
public:
    FvmExecution(
        const Mesh&, const Methods&, std::unique_ptr<ComputeBackend>);
    ~FvmExecution();
    FvmExecution(const FvmExecution&) = delete;
    FvmExecution& operator=(const FvmExecution&) = delete;
    // Internal bridge for the immediate procedural assembler.
    ComputeBackend& backend();
    const Mesh& mesh() const;
    void evaluate(math::ScalarGradient operation, VectorField& result);
    void evaluate(math::NormalGradient operation, ScalarField& result);
    void evaluate(math::ScalarDiffusionFlux operation, ScalarField& result);
    void evaluate(math::VectorGradient operation, TensorField& result);
    void evaluate(math::FaceFlux operation, ScalarField& result);
    void evaluate(math::FaceDivergence operation, ScalarField& result);
    void evaluate(math::VectorDivergence operation, ScalarField& result);
    void evaluate(math::TensorDivergence operation, VectorField& result);
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
    PerformanceCounters performance() const;


private:
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
};
FvmExecution& execution();
}  // babelsim::detail 命名空间
