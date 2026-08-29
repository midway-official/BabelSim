#pragma once

#include "babelsim/equation.h"
#include "babelsim/field.h"
#include "babelsim/methods.h"

namespace babelsim {

// Operators are local finite-volume kernels. In a decomposed run, callers
// exchange each input cell/face field with HaloExchange before invoking a
// kernel; no operator reads a hidden global field or performs an MPI collective.

void interpolate(
    const ScalarField& cell,
    ScalarField& face,
    InterpolationMethod method = InterpolationMethod::Linear);
void interpolate(
    const VectorField& cell,
    VectorField& face,
    InterpolationMethod method = InterpolationMethod::Linear);

// Linear interpolation to the owner-neighbour intersection followed by a
// gradient correction from that intersection to the geometric face centre.
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

}  // namespace babelsim
