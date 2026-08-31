#pragma once

#include "babelsim/discrete_equation.h"
#include "babelsim/field.h"
#include "babelsim/methods.h"

namespace babelsim {

// 算子是局部有限体积核。分区运行时，调用方在调用算子前用 HaloExchange 同步每个输入
// 单元/面场；算子不会读取隐藏的全局场，也不会执行 MPI 集体通信。

void interpolate(
    const ScalarField& cell,
    ScalarField& face,
    InterpolationMethod method = InterpolationMethod::Corrected,
    GradientMethod gradient_method = GradientMethod::LeastSquares);
void interpolate(
    const VectorField& cell,
    VectorField& face,
    InterpolationMethod method = InterpolationMethod::Corrected,
    GradientMethod gradient_method = GradientMethod::LeastSquares);

// 先线性插值到 owner 与 neighbour 连线交点，再由该交点的梯度修正重构到几何面中心。
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

// 返回面积积分后的面法向梯度。正交差分保持隐式所需的紧凑计算模板，非正交部分由
// 已重构的单元梯度显式补偿；扩散和压力修正共享这一几何定义。
double integratedNormalGradient(
    const ScalarField& scalar,
    const VectorField& scalar_gradient,
    Index face,
    DiffusionMethod method = DiffusionMethod::Corrected);
Vec3 integratedNormalGradient(
    const VectorField& vector,
    const TensorField& vector_gradient,
    Index face,
    DiffusionMethod method = DiffusionMethod::Corrected);

// cell 矢量按所选方法插值；face 矢量直接与面面积向量点乘。
void flux(
    const VectorField& velocity,
    ScalarField& face_flux,
    InterpolationMethod method = InterpolationMethod::Corrected,
    GradientMethod gradient_method = GradientMethod::LeastSquares);

void diffusionFlux(
    const ScalarField& face_diffusivity,
    const ScalarField& scalar,
    const VectorField& scalar_gradient,
    ScalarField& face_flux,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);

void divergence(
    const VectorField& vector,
    ScalarField& result,
    InterpolationMethod method = InterpolationMethod::Corrected,
    GradientMethod gradient_method = GradientMethod::LeastSquares);
void divergence(const ScalarField& face_flux, ScalarField& result);
void convection(
    const ScalarField& face_flux,
    const ScalarField& transported,
    ScalarField& result,
    ConvectionMethod method = ConvectionMethod::Upwind,
    InterpolationMethod interpolation_method = InterpolationMethod::Corrected,
    GradientMethod gradient_method = GradientMethod::LeastSquares);
void convection(
    const ScalarField& face_flux,
    const VectorField& transported,
    VectorField& result,
    ConvectionMethod method = ConvectionMethod::Upwind,
    InterpolationMethod interpolation_method = InterpolationMethod::Corrected,
    GradientMethod gradient_method = GradientMethod::LeastSquares);

void laplacian(
    const ScalarField& scalar,
    ScalarField& result,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);
void laplacian(
    double diffusivity,
    const ScalarField& scalar,
    ScalarField& result,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);
void laplacian(
    const ScalarField& face_diffusivity,
    const ScalarField& scalar,
    ScalarField& result,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);

void addConvection(
    ScalarDiscreteEquation& equation,
    const ScalarField& face_flux,
    const ScalarField& transported,
    ConvectionMethod method = ConvectionMethod::Upwind,
    InterpolationMethod interpolation_method = InterpolationMethod::Corrected,
    GradientMethod gradient_method = GradientMethod::LeastSquares,
    double flux_scale = 1.0);
void addConvection(
    VectorDiscreteEquation& equation,
    const ScalarField& face_flux,
    const VectorField& transported,
    ConvectionMethod method = ConvectionMethod::Upwind,
    InterpolationMethod interpolation_method = InterpolationMethod::Corrected,
    GradientMethod gradient_method = GradientMethod::LeastSquares,
    double flux_scale = 1.0);

void addDiffusion(
    ScalarDiscreteEquation& equation,
    double diffusivity,
    const ScalarField& scalar,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);
void addDiffusion(
    ScalarDiscreteEquation& equation,
    const ScalarField& face_diffusivity,
    const ScalarField& scalar,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);
void addDiffusion(
    VectorDiscreteEquation& equation,
    double diffusivity,
    const VectorField& vector,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);
void addDiffusion(
    VectorDiscreteEquation& equation,
    const ScalarField& face_diffusivity,
    const VectorField& vector,
    GradientMethod gradient_method = GradientMethod::GreenGauss,
    DiffusionMethod diffusion_method = DiffusionMethod::Corrected);

void addTimeDerivative(
    ScalarDiscreteEquation& equation,
    const ScalarField& previous,
    double dt,
    double density = 1.0,
    TimeMethod method = TimeMethod::Euler,
    const ScalarField* older = nullptr);
void addTimeDerivative(
    VectorDiscreteEquation& equation,
    const VectorField& previous,
    double dt,
    double density = 1.0,
    TimeMethod method = TimeMethod::Euler,
    const VectorField* older = nullptr);
void addTimeDerivative(
    ScalarDiscreteEquation& equation,
    const ScalarField& previous,
    double dt,
    const ScalarField& volumetric_capacity,
    TimeMethod method = TimeMethod::Euler,
    const ScalarField* older = nullptr);
void addTimeDerivative(
    VectorDiscreteEquation& equation,
    const VectorField& previous,
    double dt,
    const ScalarField& density,
    TimeMethod method = TimeMethod::Euler,
    const VectorField* older = nullptr);

}  // babelsim 命名空间
