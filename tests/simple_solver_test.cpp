#include "babelsim/incompressible.h"
#include "babelsim/mesh_io.h"

#include "test_util.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace babelsim;

namespace {

void configureChannelBoundaries(IncompressibleFields& fields) {
    const Mesh& mesh = fields.velocity.mesh();
    for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        switch (mesh.patches[static_cast<std::size_t>(patch)].kind) {
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

    Methods methods;
    methods.gradient = GradientMethod::GreenGauss;
    methods.convection = ConvectionMethod::Upwind;
    methods.diffusion = DiffusionMethod::Orthogonal;
    SimpleControl control;
    control.max_iterations = 500;
    control.velocity_relaxation = 0.5;
    control.pressure_relaxation = 0.3;
    control.continuity_tolerance = 1e-9;
    control.velocity_tolerance = 1e-7;
    control.velocity_solver.absolute_tolerance = 1e-13;
    control.velocity_solver.relative_tolerance = 1e-10;
    control.pressure_solver.absolute_tolerance = 1e-13;
    control.pressure_solver.relative_tolerance = 1e-10;

    RunTime run_time = RunTime::forMesh(mesh, simpleRunTimeControl(methods, control));
    // rho=2 同时验证 fvm::div(rho, phi, U) 的常数通量缩放路径；保持相同运动
    // 黏度以维持该回归的 Reynolds 数。
    SimpleSolver solver(run_time, fields, {2.0, 0.2}, control);
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
    for (const Vec3& velocity : fields.velocity.values()) {
        maximum_z_velocity = std::max(maximum_z_velocity, std::abs(velocity.z));
        require(isFinite(velocity), "SIMPLE velocity contains a non-finite value");
    }
    require(maximum_z_velocity < 1e-13, "nz=1 SIMPLE generated a z velocity");

    std::cout << "simple_solver_test: iterations=" << iterations
              << " mass=" << result.continuity.relative
              << " dU=" << result.relative_velocity_change << '\n';
}
