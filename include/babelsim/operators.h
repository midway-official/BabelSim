#pragma once

#include "babelsim/equation.h"
#include "babelsim/field.h"
#include "babelsim/methods.h"

namespace babelsim {

// 算子是局部有限体积核。分区运行时，调用方在调用算子前用 HaloExchange 同步每个输入
// 单元/面场；算子不会读取隐藏的全局场，也不会执行 MPI 集体通信。

void interpolate(
    const ScalarField& cell,
    ScalarField& face,
    InterpolationMethod method = InterpolationMethod::Linear);
void interpolate(
    const VectorField& cell,
    VectorField& face,
    InterpolationMethod method = InterpolationMethod::Linear);

// 先线性插值到 owner-neighbour 交点，再由该交点的梯度修正重构到几何面中心。
void reconstruct(
    const ScalarField& cell,
    const VectorField& cell_gradient,
    ScalarField& face);
void reconstruct(
    const VectorField& cell,
    const TensorField& cell_gradient,
    VectorField& face);

void gradient(
    const ScalarField& scalar,
    VectorField& result,
    GradientMethod method = GradientMethod::GreenGauss);
void gradient(
    const VectorField& vector,
    TensorField& result,
    GradientMethod method = GradientMethod::GreenGauss);

void flux(
    const VectorField& velocity,
    ScalarField& face_flux,
    InterpolationMethod method = InterpolationMethod::Linear);

void divergence(
    const VectorField& vector,
    ScalarField& result,
    InterpolationMethod method = InterpolationMethod::Linear);
void divergence(const ScalarField& face_flux, ScalarField& result);

void laplacian(
    const ScalarField& scalar,
    ScalarField& result,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);

void addConvection(
    ScalarEquation& equation,
    const ScalarField& face_flux,
    const ScalarField& transported,
    ConvectionMethod method = ConvectionMethod::Upwind);
void addConvection(
    VectorEquation& equation,
    const ScalarField& face_flux,
    const VectorField& transported,
    ConvectionMethod method = ConvectionMethod::Upwind);

void addDiffusion(
    ScalarEquation& equation,
    double diffusivity,
    const ScalarField& scalar,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);
void addDiffusion(
    ScalarEquation& equation,
    const ScalarField& face_diffusivity,
    const ScalarField& scalar,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);
void addDiffusion(
    VectorEquation& equation,
    double diffusivity,
    const VectorField& vector,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);
void addDiffusion(
    VectorEquation& equation,
    const ScalarField& face_diffusivity,
    const VectorField& vector,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);

void addTimeDerivative(
    ScalarEquation& equation,
    const ScalarField& previous,
    double dt,
    double density = 1.0,
    TimeMethod method = TimeMethod::Euler,
    const ScalarField* older = nullptr);
void addTimeDerivative(
    VectorEquation& equation,
    const VectorField& previous,
    double dt,
    double density = 1.0,
    TimeMethod method = TimeMethod::Euler,
    const VectorField* older = nullptr);

}  // babelsim 命名空间
