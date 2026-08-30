#include "babelsim/incompressible_operators.h"

#include "test_util.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace babelsim;

namespace {

Mesh affineNonOrthogonalMesh() {
    constexpr Index n = 5;
    std::vector<Vec3> points;
    points.reserve(static_cast<std::size_t>((n + 1) * (n + 1) * (n + 1)));
    for (Index k = 0; k <= n; ++k) {
        for (Index j = 0; j <= n; ++j) {
            for (Index i = 0; i <= n; ++i) {
                points.push_back({
                    static_cast<double>(i) + 0.30 * j + 0.12 * k,
                    static_cast<double>(j) + 0.18 * k,
                    static_cast<double>(k),
                });
            }
        }
    }
    return Mesh::structured({n, n, n}, std::move(points));
}

double affinePressure(const Vec3& point) {
    return 1.7 * point.x - 0.8 * point.y + 0.45 * point.z + 2.0;
}

RuntimeControl runtimeControl(
    InterpolationMethod interpolation,
    GradientMethod gradient,
    DiffusionMethod diffusion)
{
    RuntimeControl control;
    control.methods.interpolation = interpolation;
    control.methods.gradient = gradient;
    control.methods.diffusion = diffusion;
    return control;
}

}  // 匿名命名空间

int main() {
    const Mesh mesh = Mesh::cartesian({2, 2, 1}, {0, 0, 0}, {1, 1, 1});
    RunTime run_time = RunTime::forMesh(
        mesh, runtimeControl(
                  InterpolationMethod::Linear, GradientMethod::GreenGauss,
                  DiffusionMethod::Orthogonal));
    VectorField velocity(mesh, FieldLocation::Cell, "U");
    ScalarField pressure(mesh, FieldLocation::Cell, "p");
    ScalarField mobility(mesh, FieldLocation::Cell, "rAU", 0.25);
    VectorField gradient(mesh, FieldLocation::Cell, "gradP");
    ScalarField flux(mesh, FieldLocation::Face, "phi");
    for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        velocity.boundary(patch) = zeroGradient();
        pressure.boundary(patch) = zeroGradient();
    }

    MomentumInterpolationWorkspace workspace(mesh);
    MomentumInterpolation::apply(
        run_time, mesh, velocity, pressure, mobility, gradient, flux, workspace,
        InterpolationMethod::Linear, GradientMethod::GreenGauss,
        DiffusionMethod::Orthogonal);
    for (double value : flux.values()) require(near(value, 0.0), "zero state gained flux");

    ScalarField correction(mesh, FieldLocation::Cell, "pPrime");
    ScalarField face_mobility(mesh, FieldLocation::Face, "rAUFace");
    ScalarField divergence(mesh, FieldLocation::Cell, "divPhi");
    const SolveResult pressure_result = PressureCorrection::solve(
        run_time, correction, mesh, flux, mobility, pressure, face_mobility,
        divergence, false);
    require(
        pressure_result.status != SolveStatus::NumericalFailure,
        "pressure correction solve failed");
    VectorField correction_gradient(mesh, FieldLocation::Cell, "gradPPrime");
    PressureCorrection::apply(
        run_time, mesh, 0.3, pressure, velocity, flux, correction, mobility,
        correction_gradient, DiffusionMethod::Orthogonal);
    for (Index cell : mesh.owned_cells) {
        require(near(pressure[cell], 0.0) && near(velocity[cell], {}),
                "zero pressure correction changed the field");
    }

    // 三维仿射非正交网格上，修正的 Rhie-Chow 项必须和真实面梯度相抵消；
    // 这同时验证专用算子内部复用了 fvc 插值而没有退回正交近似。
    const Mesh nonorthogonal = affineNonOrthogonalMesh();
    const Index centre = nonorthogonal.cellId(2, 2, 2);
    const Vec3 exact_gradient{1.7, -0.8, 0.45};
    VectorField zero_velocity(nonorthogonal, FieldLocation::Cell, "U", {});
    ScalarField affine_pressure(nonorthogonal, FieldLocation::Cell, "p");
    ScalarField affine_mobility(nonorthogonal, FieldLocation::Cell, "rAU", 0.25);
    VectorField affine_gradient(
        nonorthogonal, FieldLocation::Cell, "gradP", exact_gradient);
    for (Index cell = 0; cell < nonorthogonal.cellCount(); ++cell) {
        affine_pressure[cell] = affinePressure(
            nonorthogonal.cell_centres[static_cast<std::size_t>(cell)]);
    }
    for (Index patch = 0;
         patch < static_cast<Index>(nonorthogonal.patches.size()); ++patch) {
        zero_velocity.boundary(patch) = zeroGradient();
        affine_pressure.boundary(patch) = zeroGradient();
    }
    RunTime corrected_run_time = RunTime::forMesh(
        nonorthogonal, runtimeControl(
                           InterpolationMethod::Corrected,
                           GradientMethod::LeastSquares,
                           DiffusionMethod::Corrected));
    RunTime orthogonal_run_time = RunTime::forMesh(
        nonorthogonal, runtimeControl(
                           InterpolationMethod::Linear,
                           GradientMethod::LeastSquares,
                           DiffusionMethod::Orthogonal));
    ScalarField corrected_flux(nonorthogonal, FieldLocation::Face, "correctedPhi");
    ScalarField orthogonal_flux(nonorthogonal, FieldLocation::Face, "orthogonalPhi");
    MomentumInterpolationWorkspace corrected_workspace(nonorthogonal);
    MomentumInterpolationWorkspace orthogonal_workspace(nonorthogonal);
    MomentumInterpolation::apply(
        corrected_run_time, nonorthogonal, zero_velocity, affine_pressure,
        affine_mobility, affine_gradient, corrected_flux, corrected_workspace,
        InterpolationMethod::Corrected, GradientMethod::LeastSquares,
        DiffusionMethod::Corrected);
    MomentumInterpolation::apply(
        orthogonal_run_time, nonorthogonal, zero_velocity, affine_pressure,
        affine_mobility, affine_gradient, orthogonal_flux, orthogonal_workspace,
        InterpolationMethod::Linear, GradientMethod::LeastSquares,
        DiffusionMethod::Orthogonal);
    double maximum_corrected_flux = 0.0;
    double maximum_orthogonal_flux = 0.0;
    for (Index face : nonorthogonal.cell_faces[static_cast<std::size_t>(centre)]) {
        maximum_corrected_flux = std::max(
            maximum_corrected_flux, std::abs(corrected_flux[face]));
        maximum_orthogonal_flux = std::max(
            maximum_orthogonal_flux, std::abs(orthogonal_flux[face]));
    }
    require(maximum_corrected_flux < 1e-12,
            "corrected momentum interpolation is not affine exact");
    require(maximum_orthogonal_flux > 1e-3,
            "momentum interpolation test does not exercise non-orthogonality");

    ScalarField affine_correction(nonorthogonal, FieldLocation::Cell, "pPrime");
    ScalarField pressure_after(nonorthogonal, FieldLocation::Cell, "pressureAfter");
    VectorField velocity_after(nonorthogonal, FieldLocation::Cell, "velocityAfter");
    ScalarField correction_flux(nonorthogonal, FieldLocation::Face, "correctionPhi");
    VectorField affine_correction_gradient(
        nonorthogonal, FieldLocation::Cell, "gradPPrime");
    for (Index cell = 0; cell < nonorthogonal.cellCount(); ++cell) {
        affine_correction[cell] = affinePressure(
            nonorthogonal.cell_centres[static_cast<std::size_t>(cell)]);
    }
    for (Index patch = 0;
         patch < static_cast<Index>(nonorthogonal.patches.size()); ++patch) {
        affine_correction.boundary(patch) = zeroGradient();
        pressure_after.boundary(patch) = zeroGradient();
        velocity_after.boundary(patch) = zeroGradient();
    }
    PressureCorrection::apply(
        corrected_run_time, nonorthogonal, 0.3, pressure_after, velocity_after,
        correction_flux, affine_correction, affine_mobility,
        affine_correction_gradient, DiffusionMethod::Corrected);
    for (Index face : nonorthogonal.cell_faces[static_cast<std::size_t>(centre)]) {
        const std::size_t index = static_cast<std::size_t>(face);
        const double expected = -0.25 * dot(
            exact_gradient, nonorthogonal.face_area_vectors[index]);
        require(near(correction_flux[face], expected, 1e-12),
                "pressure correction missed the non-orthogonal face gradient");
    }
    std::cout << "specialized_operator_test: fvm/fvc-backed SIMPLE operators passed\n";
}
