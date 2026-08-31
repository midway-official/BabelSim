#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/runtime.h"

#include "babelsim/fvm.h"

#include "test_util.h"

#include <cmath>
#include <iostream>

using namespace babelsim;

int main() {
    const Mesh mesh = Mesh::cartesian({1, 1, 1}, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    ScalarField temperature(mesh, FieldLocation::Cell, "T", 1.0);
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
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

    require(run_time.loop(), "heat run did not start its time step");
    const SolveResult result = solve(
        fvm::ddt(temperature) == fvm::laplacian(1.0, temperature));
    require(result.converged(), "heat equation did not converge");
    require(!run_time.loop(), "heat run used an unexpected number of time steps");
    require(
        near(detail::fieldData(temperature)[0], 10.0 / 22.0, 1e-12),
        "implicit heat equation does not match the one-cell FVM result");
    require(run_time.step() == 1, "runtime did not advance exactly one step");

    // 系数也可以是 cell Field。RunTime 自动将扩散系数插值到面，同时保持 ddt
    // 系数在单元上；Solver 不需要编写插值、矩阵组装或并行同步。
    ScalarField variable_temperature(mesh, FieldLocation::Cell, "Tv", 1.0);
    ScalarField heat_capacity(mesh, FieldLocation::Cell, "rhoCp", 2.0);
    ScalarField conductivity(mesh, FieldLocation::Cell, "k", 1.0);
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        variable_temperature.boundary(patch) = fixedValue(0.0);
    }
    const SolveResult variable_result = solve(
        fvm::ddt(heat_capacity, variable_temperature) ==
            fvm::laplacian(conductivity, variable_temperature));
    require(variable_result.converged(), "variable-coefficient heat solve did not converge");
    require(
        near(detail::fieldData(variable_temperature)[0], 20.0 / 32.0, 1e-12),
        "field time/diffusion coefficients do not match the one-cell FVM result");

    // 同一 API 还接受材料场和热源场。温度相关材料只需在每步前更新 conductivity，
    // 不需要改动 FVM、运行时或 MPI 层。
    ScalarField field_temperature(mesh, FieldLocation::Cell, "Tf", 1.0);
    ScalarField field_source(mesh, FieldLocation::Cell, "Q", 0.0);
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        field_temperature.boundary(patch) = fixedValue(0.0);
    }
    const SolveResult field_result = solve(
        fvm::ddt(heat_capacity, field_temperature) ==
            fvm::laplacian(conductivity, field_temperature) + fvm::source(field_source));
    require(field_result.converged(), "Field-material heat step did not converge");

    std::cout << "heat_solver_test: T=" << detail::fieldData(temperature)[0]
              << " residual=" << result.relative_residual << '\n';
}
