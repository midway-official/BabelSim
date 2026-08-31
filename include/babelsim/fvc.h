#pragma once

#include "babelsim/field.h"

namespace babelsim::fvc {

// fvc 描述“立即计算为场”的显式有限体积运算。描述对象只保存 Field 引用；
// evaluate()/subtract() 才执行：所有参与进程以相同顺序调用，后端同步所有输入、
// 选择离散方法，并同步写出的结果。调用者不负责 halo；同位输入与结果不得别名。
// 单面/单元局部核不属于公开 fvc API。
struct ScalarGradient {
    const ScalarField& field;
};

// 每个面的外法向梯度；执行时自动重构单元梯度、同步输入并修正非正交性。
struct NormalGradient {
    const ScalarField& field;
};

struct VectorGradient {
    const VectorField& field;
};

struct FaceFlux {
    const VectorField& velocity;
};

// 标量扩散面通量：coefficient * Sf·grad(field)。coefficient 可位于 cell 或 face，
// FVM 执行层负责同步、插值和梯度工作区。
struct ScalarDiffusionFlux {
    const ScalarField& coefficient;
    const ScalarField& field;
};

struct FaceDivergence {
    const ScalarField& flux;
};

struct VectorDivergence {
    const VectorField& field;
};

struct ScalarConvection {
    const ScalarField& flux;
    const ScalarField& field;
};

struct VectorConvection {
    const ScalarField& flux;
    const VectorField& field;
};

struct ScalarInterpolation {
    const ScalarField& field;
};

struct VectorInterpolation {
    const VectorField& field;
};

struct ScalarReconstruction {
    const ScalarField& field;
    const VectorField& gradient;
};

struct VectorReconstruction {
    const VectorField& field;
    const TensorField& gradient;
};

struct ScalarLaplacian {
    const ScalarField& field;
    double coefficient = 1.0;
    const ScalarField* coefficient_field = nullptr;
};

inline ScalarGradient grad(const ScalarField& field) { return {field}; }
inline NormalGradient normalGradient(const ScalarField& field) { return {field}; }
inline VectorGradient grad(const VectorField& field) { return {field}; }
inline FaceFlux flux(const VectorField& velocity) { return {velocity}; }
inline ScalarDiffusionFlux flux(
    const ScalarField& coefficient,
    const ScalarField& field)
{
    return {coefficient, field};
}
inline FaceDivergence div(const ScalarField& flux) { return {flux}; }
inline VectorDivergence div(const VectorField& field) { return {field}; }
inline ScalarConvection div(const ScalarField& flux, const ScalarField& field) {
    return {flux, field};
}
inline VectorConvection div(const ScalarField& flux, const VectorField& field) {
    return {flux, field};
}
inline ScalarInterpolation interpolate(const ScalarField& field) { return {field}; }
inline VectorInterpolation interpolate(const VectorField& field) { return {field}; }
inline ScalarReconstruction reconstruct(
    const ScalarField& field,
    const VectorField& gradient)
{
    return {field, gradient};
}
inline VectorReconstruction reconstruct(
    const VectorField& field,
    const TensorField& gradient)
{
    return {field, gradient};
}
inline ScalarLaplacian laplacian(const ScalarField& field) { return {field}; }
inline ScalarLaplacian laplacian(double diffusivity, const ScalarField& field) {
    return {field, diffusivity, nullptr};
}
inline ScalarLaplacian laplacian(
    const ScalarField& diffusivity,
    const ScalarField& field)
{
    return {field, 1.0, &diffusivity};
}

void evaluate(ScalarGradient operation, VectorField& result);
void evaluate(NormalGradient operation, ScalarField& result);
void evaluate(ScalarDiffusionFlux operation, ScalarField& result);
void evaluate(VectorGradient operation, TensorField& result);
void evaluate(FaceFlux operation, ScalarField& result);
void evaluate(FaceDivergence operation, ScalarField& result);
void evaluate(VectorDivergence operation, ScalarField& result);
void evaluate(ScalarConvection operation, ScalarField& result);
void evaluate(VectorConvection operation, VectorField& result);
void evaluate(ScalarInterpolation operation, ScalarField& result);
void evaluate(VectorInterpolation operation, VectorField& result);
void evaluate(ScalarReconstruction operation, ScalarField& result);
void evaluate(VectorReconstruction operation, VectorField& result);
void evaluate(ScalarLaplacian operation, ScalarField& result);
void subtract(
    const ScalarField& coefficient,
    ScalarGradient operation,
    VectorField& target);
void subtract(ScalarDiffusionFlux operation, ScalarField& target);

}  // babelsim::fvc 命名空间
