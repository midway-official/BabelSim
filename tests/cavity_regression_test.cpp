#include "internal/field_access.h"
#include "support/simple_reference.h"
#include "babelsim/runtime.h"

#include "test_util.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace babelsim;

int main() {
    constexpr Index n = 12;
    auto patches = boxPatches();
    for (Index side : {0, 1, 2, 3}) {
        patches[static_cast<std::size_t>(side)].kind = PatchKind::Wall;
    }
    patches[static_cast<std::size_t>(4)].kind = PatchKind::Symmetry;
    patches[static_cast<std::size_t>(5)].kind = PatchKind::Symmetry;
    const Mesh mesh = makeHexBox(
        {n, n, 1}, {0, 0, 0}, {1, 1, 1}, patches);
    IncompressibleFields fields(mesh);
    for (Index side : {0, 1, 2}) {
        fields.velocity.setBoundary(
            static_cast<Index>(side),
            BoundaryCondition<Vec3>::fixedValue({}));
    }
    fields.velocity.setBoundary(
        static_cast<Index>(3),
        BoundaryCondition<Vec3>::fixedValue({1.0, 0.0, 0.0}));
    for (Index side : {4, 5}) {
        fields.velocity.setBoundary(
            static_cast<Index>(side),
            BoundaryCondition<Vec3>::symmetry());
        fields.pressure.setBoundary(
            static_cast<Index>(side),
            BoundaryCondition<double>::symmetry());
    }

    RuntimeControl run_control;
    SimpleControl control;
    control.max_iterations = 3000;
    control.velocity_relaxation = 0.5;
    control.pressure_relaxation = 0.3;
    control.continuity_tolerance = 1e-8;
    control.velocity_tolerance = 1e-6;
    control.momentum_equation.spatial.gradient = GradientMethod::GreenGauss;
    control.momentum_equation.spatial.convection = ConvectionMethod::Upwind;
    control.momentum_equation.spatial.diffusion = DiffusionMethod::Orthogonal;
    control.pressure_equation.spatial = control.momentum_equation.spatial;
    control.pressure_equation.linear.solver = LinearSolverType::ConjugateGradient;
    control.pressure_equation.linear.preconditioner = PreconditionerType::IncompleteCholesky;
    control.momentum_equation.linear.absolute_tolerance = 1e-16;
    control.momentum_equation.linear.relative_tolerance = 1e-10;
    control.pressure_equation.linear.absolute_tolerance = 1e-16;
    control.pressure_equation.linear.relative_tolerance = 1e-10;

    RunTime run_time = RunTime::forMesh(mesh, run_control);
    int iterations = 0;
    const auto result = solveIncompressible(fields, {1.0, 0.01}, control, &iterations);
    require(result.healthy, "SIMPLE produced a numerical failure");
    require(result.converged, "closed-cavity SIMPLE did not converge");
    const Index centre = hexCellIndex(n / 2 - 1, n / 2 - 1, 0, n, n);
    require(detail::fieldData(fields.velocity)[centre].x < -0.05, "cavity primary vortex is missing");

    const Vec3 upper_left = detail::fieldData(fields.velocity)[hexCellIndex(0, n - 3, 0, n, n)];
    const Vec3 upper_right = detail::fieldData(fields.velocity)[hexCellIndex(n - 1, n - 3, 0, n, n)];
    require(
        upper_left.y > 0.0 && upper_right.y < 0.0,
        "cavity circulation direction is incorrect");
    std::cout << "cavity_regression_test: iterations=" << iterations
              << " mass=" << result.continuity.relative
              << " dU=" << result.relative_velocity_change
              << " centreU=" << detail::fieldData(fields.velocity)[centre].x
              << " pLin=" << result.pressure.relative_residual << '\n';
}
