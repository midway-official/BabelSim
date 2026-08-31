#include "internal/field_access.h"
#include "babelsim/simple.h"
#include "babelsim/simple_control.h"
#include "babelsim/runtime.h"

#include "test_util.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace babelsim;

int main() {
    constexpr Index n = 12;
    auto patches = defaultPatches();
    for (Side side : {Side::XMin, Side::XMax, Side::YMin, Side::YMax}) {
        patches[static_cast<std::size_t>(side)].kind = PatchKind::Wall;
    }
    patches[static_cast<std::size_t>(Side::ZMin)].kind = PatchKind::Symmetry;
    patches[static_cast<std::size_t>(Side::ZMax)].kind = PatchKind::Symmetry;
    const Mesh mesh = Mesh::cartesian(
        {n, n, 1}, {0, 0, 0}, {1, 1, 1}, patches);
    IncompressibleFields fields(mesh);
    for (Side side : {Side::XMin, Side::XMax, Side::YMin}) {
        fields.velocity.setBoundary(
            static_cast<Index>(side),
            BoundaryCondition<Vec3>::fixedValue({}));
    }
    fields.velocity.setBoundary(
        static_cast<Index>(Side::YMax),
        BoundaryCondition<Vec3>::fixedValue({1.0, 0.0, 0.0}));
    for (Side side : {Side::ZMin, Side::ZMax}) {
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
    SimpleSolver solver(run_time, fields, {1.0, 0.01}, control);
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
    const Index centre = mesh.cellId(n / 2 - 1, n / 2 - 1, 0);
    require(detail::fieldData(fields.velocity)[centre].x < -0.05, "cavity primary vortex is missing");

    const Vec3 upper_left = detail::fieldData(fields.velocity)[mesh.cellId(0, n - 3, 0)];
    const Vec3 upper_right = detail::fieldData(fields.velocity)[mesh.cellId(n - 1, n - 3, 0)];
    require(
        upper_left.y > 0.0 && upper_right.y < 0.0,
        "cavity circulation direction is incorrect");
    std::cout << "cavity_regression_test: iterations=" << iterations
              << " mass=" << result.continuity.relative
              << " dU=" << result.relative_velocity_change
              << " centreU=" << detail::fieldData(fields.velocity)[centre].x
              << " pLin=" << result.pressure.relative_residual << '\n';
}
