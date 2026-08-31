#include "internal/mesh_access.h"
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
    constexpr Index n = 6;
    auto patches = defaultPatches();
    for (auto& patch : patches) {
        patch.kind = PatchKind::Wall;
    }
    const Mesh mesh = Mesh::cartesian(
        {n, n, n}, {0, 0, 0}, {1, 1, 1}, patches);
    IncompressibleFields fields(mesh);
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        fields.velocity.setBoundary(
            patch, BoundaryCondition<Vec3>::fixedValue({}));
    }
    fields.velocity.setBoundary(
        static_cast<Index>(Side::YMax),
        BoundaryCondition<Vec3>::fixedValue({1.0, 0.0, 0.0}));

    RuntimeControl run_control;
    run_control.scalar_solver.solver = LinearSolverType::ConjugateGradient;
    run_control.scalar_solver.preconditioner = PreconditionerType::IncompleteCholesky;
    Methods& methods = run_control.methods;
    methods.gradient = GradientMethod::GreenGauss;
    methods.convection = ConvectionMethod::Upwind;
    methods.diffusion = DiffusionMethod::Orthogonal;
    SimpleControl control;
    control.max_iterations = 2000;
    control.velocity_relaxation = 0.5;
    control.pressure_relaxation = 0.3;
    control.continuity_tolerance = 1e-8;
    control.velocity_tolerance = 1e-6;
    run_control.vector_solver.absolute_tolerance = 1e-15;
    run_control.vector_solver.relative_tolerance = 1e-9;
    run_control.scalar_solver.absolute_tolerance = 1e-15;
    run_control.scalar_solver.relative_tolerance = 1e-9;

    RunTime run_time = RunTime::forMesh(mesh, run_control);
    SimpleSolver solver(run_time, fields, {1.0, 0.01}, control);
    SimpleIterationResult result;
    int iterations = 0;
    for (int iteration = 1; iteration <= control.max_iterations; ++iteration) {
        result = solver.iterate();
        iterations = iteration;
        require(result.healthy, "3D cavity SIMPLE became unhealthy");
        if (result.converged) {
            break;
        }
    }
    require(result.converged, "3D cavity SIMPLE did not converge");

    const Index lower_centre = mesh.cellId(n / 2 - 1, n / 2 - 1, n / 2 - 1);
    const Index upper_centre = mesh.cellId(n / 2 - 1, n / 2 - 1, n / 2);
    const Vec3 lower = detail::fieldData(fields.velocity)[lower_centre];
    const Vec3 upper = detail::fieldData(fields.velocity)[upper_centre];
    require(
        0.5 * (lower.x + upper.x) < -0.02,
        "3D cavity primary vortex is missing");
    require(
        near(lower.x, upper.x, 2e-5) && near(lower.y, upper.y, 2e-5) &&
            near(lower.z, -upper.z, 2e-5),
        "3D cavity violates centre-plane symmetry");

    double maximum_spanwise_velocity = 0.0;
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        maximum_spanwise_velocity = std::max(
            maximum_spanwise_velocity, std::abs(detail::fieldData(fields.velocity)[cell].z));
    }
    require(
        maximum_spanwise_velocity > 1e-4,
        "3D cavity collapsed to a duplicated 2D solution");
    require(
        result.continuity.relative <= control.continuity_tolerance,
        "3D cavity is not mass conservative");

    std::cout << "cavity_3d_test: cells=" << mesh.cellCount()
              << " iterations=" << iterations
              << " mass=" << result.continuity.relative
              << " centreU=" << 0.5 * (lower.x + upper.x)
              << " maxW=" << maximum_spanwise_velocity << '\n';
}
