#include "babelsim/thermal.h"

#include "babelsim/fvm.h"

#include "test_util.h"

#include <cmath>
#include <iostream>

using namespace babelsim;

int main() {
    const Mesh mesh = Mesh::cartesian({1, 1, 1}, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    ScalarField temperature(mesh, FieldLocation::Cell, "T", 1.0);
    for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        temperature.setBoundary(patch, BoundaryCondition<double>::fixedValue(0.0));
    }

    RuntimeControl control;
    control.methods.time = TimeMethod::Euler;
    control.methods.gradient = GradientMethod::GreenGauss;
    control.methods.diffusion = DiffusionMethod::Orthogonal;
    control.time = {0.0, 0.1, 0.1};
    control.scalar_solver.absolute_tolerance = 1e-14;
    control.scalar_solver.relative_tolerance = 1e-12;
    RunTime run_time = RunTime::forMesh(mesh, control);

    const HeatResult result = solveTransientHeat(
        run_time, temperature, {1.0, 1.0, 1.0});
    require(result.converged, "heat equation did not converge");
    require(result.steps == 1, "heat run used an unexpected number of time steps");
    require(
        near(temperature[0], 10.0 / 22.0, 1e-12),
        "implicit heat equation does not match the one-cell FVM result");
    require(run_time.step() == 1, "runtime did not advance exactly one step");

    // 系数也可以是 cell Field。RunTime 自动将扩散系数插值到面，同时保持 ddt
    // 系数在单元上；Solver 不需要编写插值、矩阵组装或并行同步。
    ScalarField variable_temperature(mesh, FieldLocation::Cell, "Tv", 1.0);
    ScalarField heat_capacity(mesh, FieldLocation::Cell, "rhoCp", 2.0);
    ScalarField conductivity(mesh, FieldLocation::Cell, "k", 1.0);
    for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        variable_temperature.boundary(patch) = fixedValue(0.0);
    }
    RunTime variable_run_time = RunTime::forMesh(mesh, control);
    require(variable_run_time.loop(), "variable-coefficient heat step was not started");
    const SolveResult variable_result = solve(
        variable_run_time,
        fvm::ddt(heat_capacity, variable_temperature) ==
            fvm::laplacian(conductivity, variable_temperature));
    require(variable_result.converged(), "variable-coefficient heat solve did not converge");
    require(
        near(variable_temperature[0], 20.0 / 32.0, 1e-12),
        "field time/diffusion coefficients do not match the one-cell FVM result");

    std::cout << "heat_solver_test: T=" << temperature[0]
              << " residual=" << result.linear.relative_residual << '\n';
}
