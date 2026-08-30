#pragma once

#include "babelsim/field.h"
#include "babelsim/methods.h"

namespace babelsim::fvc {

// fvc 描述“立即计算为场”的显式有限体积运算。描述对象只保存 Field 引用；
// RunTime::evaluate() 才会选择方法、同步 halo 并写入调用者提供的结果场。
struct ScalarGradient {
    const ScalarField& field;
};

struct VectorGradient {
    const VectorField& field;
};

struct FaceFlux {
    const VectorField& velocity;
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

// 面法向梯度的面积积分，是扩散和压力修正共享的显式面量。
double integratedNormalGradient(
    const ScalarField& field,
    const VectorField& gradient,
    Index face,
    DiffusionMethod method = DiffusionMethod::Corrected);
Vec3 integratedNormalGradient(
    const VectorField& field,
    const TensorField& gradient,
    Index face,
    DiffusionMethod method = DiffusionMethod::Corrected);

inline ScalarGradient grad(const ScalarField& field) { return {field}; }
inline VectorGradient grad(const VectorField& field) { return {field}; }
inline FaceFlux flux(const VectorField& velocity) { return {velocity}; }
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

}  // babelsim::fvc 命名空间
