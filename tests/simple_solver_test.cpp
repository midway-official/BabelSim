#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "physics/simple/algorithm.h"
#include "babelsim/runtime.h"
#include "babelsim/mesh_io.h"

#include "test_util.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace babelsim;

namespace {

void configureChannelBoundaries(IncompressibleFields& fields) {
    const Mesh& mesh = fields.velocity.mesh();
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        switch (detail::meshData(mesh).patches[static_cast<std::size_t>(patch)].kind) {
            case PatchKind::Inlet:
                fields.velocity.setBoundary(
                    patch, BoundaryCondition<Vec3>::fixedValue({1.0, 0.0, 0.0}));
                fields.pressure.setBoundary(
                    patch, BoundaryCondition<double>::zeroGradient());
                break;
            case PatchKind::Outlet:
                fields.velocity.setBoundary(
                    patch, BoundaryCondition<Vec3>::zeroGradient());
                fields.pressure.setBoundary(
                    patch, BoundaryCondition<double>::fixedValue(0.0));
                break;
            case PatchKind::Wall:
                fields.velocity.setBoundary(
                    patch, BoundaryCondition<Vec3>::fixedValue({}));
                fields.pressure.setBoundary(
                    patch, BoundaryCondition<double>::zeroGradient());
                break;
            case PatchKind::Symmetry:
                fields.velocity.setBoundary(
                    patch, BoundaryCondition<Vec3>::symmetry());
                fields.pressure.setBoundary(
                    patch, BoundaryCondition<double>::symmetry());
                break;
            case PatchKind::Generic:
            case PatchKind::Processor:
                fields.velocity.setBoundary(
                    patch, BoundaryCondition<Vec3>::zeroGradient());
                fields.pressure.setBoundary(
                    patch, BoundaryCondition<double>::zeroGradient());
                break;
        }
    }
}

}  // 匿名命名空间

int main() {
    const Mesh mesh = readMeshFile("tests/data/babelsim_channel.mesh");
    IncompressibleFields fields(mesh);
    configureChannelBoundaries(fields);

    bool rejected = false;
    try { SteadySimpleAlgorithm without_run(fields, {1.0, 0.01}, {}); }
    catch (const std::logic_error&) { rejected = true; }
    require(rejected, "SIMPLE accepted fields without an active execution domain");

    RuntimeControl run_control;
    run_control.scalar_solver.solver = LinearSolverType::ConjugateGradient;
    run_control.scalar_solver.preconditioner = PreconditionerType::IncompleteCholesky;
    Methods& methods = run_control.methods;
    methods.gradient = GradientMethod::GreenGauss;
    methods.convection = ConvectionMethod::Upwind;
    methods.diffusion = DiffusionMethod::Orthogonal;
    SimpleControl control;
    control.max_iterations = 500;
    control.velocity_relaxation = 0.5;
    control.pressure_relaxation = 0.3;
    control.continuity_tolerance = 1e-9;
    control.velocity_tolerance = 1e-7;
    run_control.vector_solver.absolute_tolerance = 1e-13;
    run_control.vector_solver.relative_tolerance = 1e-10;
    run_control.scalar_solver.absolute_tolerance = 1e-13;
    run_control.scalar_solver.relative_tolerance = 1e-10;

    RunTime run_time = RunTime::forMesh(mesh, run_control);
    const Mesh other_mesh = Mesh::cartesian({2, 2, 1}, {}, {1, 1, 1});
    IncompressibleFields other_fields(other_mesh);
    rejected = false;
    try { SteadySimpleAlgorithm wrong_mesh(other_fields, {1.0, 0.01}, control); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "SIMPLE accepted fields from a different execution domain");
    // rho=2 同时验证 eqn::div(rho, phi, U) 的常数通量缩放路径；保持相同运动
    // 黏度以维持该回归的 Reynolds 数。
    SteadySimpleAlgorithm solver(fields, {2.0, 0.2}, control);
    SimpleIterationResult result;
    int iterations = 0;
    for (int iteration = 1; iteration <= control.max_iterations; ++iteration) {
        result = solver.iterate();
        iterations = iteration;
        require(result.healthy, "SIMPLE produced a numerical failure");
        if (result.converged) {
            break;
        }
    }
    require(result.converged, "SIMPLE did not converge on the native channel");
    require(
        result.continuity.relative <= control.continuity_tolerance,
        "SIMPLE converged without satisfying continuity");

    double maximum_z_velocity = 0.0;
    for (const Vec3& velocity : detail::fieldValues(fields.velocity)) {
        maximum_z_velocity = std::max(maximum_z_velocity, std::abs(velocity.z));
        require(isFinite(velocity), "SIMPLE velocity contains a non-finite value");
    }
    require(maximum_z_velocity < 1e-13, "nz=1 SIMPLE generated a z velocity");

    std::cout << "simple_solver_test: iterations=" << iterations
              << " mass=" << result.continuity.relative
              << " dU=" << result.relative_velocity_change << '\n';
}
