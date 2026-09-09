#include "internal/field_access.h"
#include "physics/simple/algorithm.h"
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
    run_control.scalar_solver.solver = LinearSolverType::ConjugateGradient;
    run_control.scalar_solver.preconditioner = PreconditionerType::IncompleteCholesky;
    Methods& methods = run_control.methods;
    methods.gradient = GradientMethod::GreenGauss;
    methods.convection = ConvectionMethod::Upwind;
    methods.diffusion = DiffusionMethod::Orthogonal;
    SimpleControl control;
    control.max_iterations = 3000;
    control.velocity_relaxation = 0.5;
    control.pressure_relaxation = 0.3;
    control.continuity_tolerance = 1e-8;
    control.velocity_tolerance = 1e-6;
    run_control.vector_solver.absolute_tolerance = 1e-16;
    run_control.vector_solver.relative_tolerance = 1e-10;
    run_control.scalar_solver.absolute_tolerance = 1e-16;
    run_control.scalar_solver.relative_tolerance = 1e-10;

    RunTime run_time = RunTime::forMesh(mesh, run_control);
    SteadySimpleAlgorithm solver(fields, {1.0, 0.01}, control);
    SimpleIterationResult result;
    int iterations = 0;
    for (int iteration = 1; iteration <= control.max_iterations; ++iteration) {
        result = solver.iterate();
        iterations = iteration;
        require(result.healthy, "closed-cavity SIMPLE became unhealthy");
        if (result.converged) {
            break;
        }
    }
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
