#include "babelsim/incompressible_operators.h"
#include "babelsim/operators.h"

#include "test_util.h"

#include <algorithm>
#include <array>
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

}  // 匿名命名空间

int main() {
    const Mesh mesh = Mesh::cartesian({2, 2, 1}, {0, 0, 0}, {1, 1, 1});
    VectorField velocity(mesh, FieldLocation::Cell, "U");
    ScalarField pressure(mesh, FieldLocation::Cell, "p");
    ScalarField mobility(mesh, FieldLocation::Cell, "rAU", 0.25);
    VectorField gradient(mesh, FieldLocation::Cell, "gradP");
    ScalarField flux(mesh, FieldLocation::Face, "phi");
    for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        velocity.setBoundary(patch, BoundaryCondition<Vec3>::zeroGradient());
        pressure.setBoundary(patch, BoundaryCondition<double>::zeroGradient());
    }

    MomentumInterpolation::apply(
        mesh, velocity, pressure, mobility, gradient, flux,
        InterpolationMethod::Linear, GradientMethod::GreenGauss,
        DiffusionMethod::Orthogonal);
    for (double value : flux.values()) require(near(value, 0.0), "zero state gained flux");

    ScalarEquation equation(mesh);
    ScalarField correction(mesh, FieldLocation::Cell, "pPrime");
    PressureCorrection::assemble(
        equation, correction, mesh, flux, mobility, pressure,
        nullptr, DiffusionMethod::Orthogonal, false, ParallelContext{});
    for (Index cell : mesh.owned_cells) {
        require(equation.diagonal[static_cast<std::size_t>(cell)] > 0.0,
                "pressure correction has no positive diagonal");
    }
    VectorField correction_gradient(mesh, FieldLocation::Cell, "gradPPrime");
    PressureCorrection::apply(
        mesh, 0.3, pressure, velocity, flux, correction, correction_gradient, mobility,
        DiffusionMethod::Orthogonal);
    for (Index cell : mesh.owned_cells) {
        require(near(pressure[cell], 0.0) && near(velocity[cell], {}),
                "zero pressure correction changed the field");
    }

    // 在三维仿射非正交网格上，Rhie-Chow 的插值压力响应必须与直接面压力梯度
    // 对仿射压力场精确抵消；正交近似则应留下可测的交叉扩散误差。
    const Mesh nonorthogonal = affineNonOrthogonalMesh();
    const Index centre = nonorthogonal.cellId(2, 2, 2);
    const Vec3 exact_gradient{1.7, -0.8, 0.45};
    VectorField zero_velocity(
        nonorthogonal, FieldLocation::Cell, "U", {});
    ScalarField affine_pressure(
        nonorthogonal, FieldLocation::Cell, "p");
    ScalarField affine_mobility(
        nonorthogonal, FieldLocation::Cell, "rAU", 0.25);
    VectorField affine_gradient(
        nonorthogonal, FieldLocation::Cell, "gradP", exact_gradient);
    for (Index cell = 0; cell < nonorthogonal.cellCount(); ++cell) {
        affine_pressure[cell] = affinePressure(
            nonorthogonal.cell_centres[static_cast<std::size_t>(cell)]);
    }
    for (Index patch = 0;
         patch < static_cast<Index>(nonorthogonal.patches.size()); ++patch) {
        zero_velocity.setBoundary(
            patch, BoundaryCondition<Vec3>::zeroGradient());
        affine_pressure.setBoundary(
            patch, BoundaryCondition<double>::zeroGradient());
    }
    ScalarField corrected_flux(
        nonorthogonal, FieldLocation::Face, "correctedPhi");
    ScalarField orthogonal_flux(
        nonorthogonal, FieldLocation::Face, "orthogonalPhi");
    MomentumInterpolation::apply(
        nonorthogonal, zero_velocity, affine_pressure, affine_mobility,
        affine_gradient, corrected_flux, InterpolationMethod::Corrected,
        GradientMethod::LeastSquares, DiffusionMethod::Corrected);
    MomentumInterpolation::apply(
        nonorthogonal, zero_velocity, affine_pressure, affine_mobility,
        affine_gradient, orthogonal_flux, InterpolationMethod::Linear,
        GradientMethod::LeastSquares, DiffusionMethod::Orthogonal);
    double maximum_corrected_flux = 0.0;
    double maximum_orthogonal_flux = 0.0;
    for (Index face :
         nonorthogonal.cell_faces[static_cast<std::size_t>(centre)]) {
        maximum_corrected_flux = std::max(
            maximum_corrected_flux, std::abs(corrected_flux[face]));
        maximum_orthogonal_flux = std::max(
            maximum_orthogonal_flux, std::abs(orthogonal_flux[face]));
    }
    require(
        maximum_corrected_flux < 1e-12,
        "corrected momentum interpolation is not affine exact");
    require(
        maximum_orthogonal_flux > 1e-3,
        "momentum interpolation test does not exercise non-orthogonality");

    ScalarField affine_correction(
        nonorthogonal, FieldLocation::Cell, "pPrime");
    VectorField affine_correction_gradient(
        nonorthogonal, FieldLocation::Cell, "gradPPrime", exact_gradient);
    ScalarField pressure_after(
        nonorthogonal, FieldLocation::Cell, "pressureAfter");
    VectorField velocity_after(
        nonorthogonal, FieldLocation::Cell, "velocityAfter");
    ScalarField correction_flux(
        nonorthogonal, FieldLocation::Face, "correctionPhi");
    for (Index cell = 0; cell < nonorthogonal.cellCount(); ++cell) {
        affine_correction[cell] = affinePressure(
            nonorthogonal.cell_centres[static_cast<std::size_t>(cell)]);
    }
    for (Index patch = 0;
         patch < static_cast<Index>(nonorthogonal.patches.size()); ++patch) {
        affine_correction.setBoundary(
            patch, BoundaryCondition<double>::zeroGradient());
        pressure_after.setBoundary(
            patch, BoundaryCondition<double>::zeroGradient());
        velocity_after.setBoundary(
            patch, BoundaryCondition<Vec3>::zeroGradient());
    }
    PressureCorrection::apply(
        nonorthogonal, 0.3, pressure_after, velocity_after, correction_flux,
        affine_correction, affine_correction_gradient, affine_mobility,
        DiffusionMethod::Corrected);
    for (Index face :
         nonorthogonal.cell_faces[static_cast<std::size_t>(centre)]) {
        const auto f = static_cast<std::size_t>(face);
        const double expected = -0.25 * dot(
            exact_gradient, nonorthogonal.face_area_vectors[f]);
        require(
            near(correction_flux[face], expected, 1e-12),
            "pressure correction missed the non-orthogonal face gradient");
    }

    ScalarField zero_flux(
        nonorthogonal, FieldLocation::Face, "zeroPhi");
    ScalarEquation nonorthogonal_equation(nonorthogonal);
    PressureCorrection::assemble(
        nonorthogonal_equation, affine_correction, nonorthogonal, zero_flux,
        affine_mobility, affine_pressure, &affine_correction_gradient,
        DiffusionMethod::Corrected, true, ParallelContext{});
    double row_residual =
        nonorthogonal_equation.diagonal[static_cast<std::size_t>(centre)] *
            affine_correction[centre] -
        nonorthogonal_equation.source[static_cast<std::size_t>(centre)];
    for (Index face :
         nonorthogonal.cell_faces[static_cast<std::size_t>(centre)]) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = nonorthogonal.face_owner[f];
        const Index neighbour = nonorthogonal.face_neighbour[f];
        row_residual += owner == centre
            ? nonorthogonal_equation.upper[f] * affine_correction[neighbour]
            : nonorthogonal_equation.lower[f] * affine_correction[owner];
    }
    require(
        std::abs(row_residual) < 1e-12,
        "corrected pressure equation does not preserve an affine solution");
    std::cout << "specialized_operator_test: momentum interpolation and pressure correction passed\n";
}
