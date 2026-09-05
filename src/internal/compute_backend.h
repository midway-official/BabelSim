#pragma once

#include "babelsim/discrete_equation.h"
#include "babelsim/field.h"
#include "babelsim/solver_control.h"

#include <array>
#include <memory>

namespace babelsim {
struct ParallelContext;

namespace detail {

// 数值前端与计算后端之间的唯一执行契约。调用粒度是一整个同步、归约或线性方程，
// 不进入 cell/face 热循环；因此可以替换 MPI/线性代数实现而不改变 FVM 数值代码。
class ComputeBackend {
public:
    virtual ~ComputeBackend() = default;

    virtual void synchronize(ScalarField& field) = 0;
    virtual void synchronize(VectorField& field) = 0;
    virtual void synchronize(TensorField& field) = 0;

    virtual void sum(const double* local, double* global, int count) const = 0;
    virtual void maximum(const double* local, double* global, int count) const = 0;
    virtual bool all(bool local_condition) const = 0;

    virtual SolveResult solve(
        const ScalarDiscreteEquation& equation, ScalarField& unknown) = 0;
    virtual std::array<SolveResult, 3> solve(
        const VectorDiscreteEquation& equation, VectorField& unknown) = 0;
};

// 默认实现由一个独立翻译单元提供。构建时替换该实现即可更换计算后端，
// RunTime、FVM 和 Physics 均不需要修改。
std::unique_ptr<ComputeBackend> makeComputeBackend(
    const Mesh& mesh,
    const LinearSolverConfig& scalar_solver,
    const LinearSolverConfig& vector_solver,
    ParallelContext parallel);

}  // detail 命名空间
}  // babelsim 命名空间
