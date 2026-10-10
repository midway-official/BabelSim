#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/equ.h"
#include "babelsim/runtime.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

namespace {

void checkHeatEquation() {
    const Mesh mesh = makeHexBox({1, 1, 1}, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    ScalarField temperature(mesh, FieldLocation::Cell, "T", 1.0);
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        temperature.setBoundary(patch, BoundaryCondition<double>::fixedValue(0.0));
    }

    RuntimeControl control;
    control.methods.time = TimeMethod::Euler;
    control.time = {0.0, 0.1, 0.1};
    RunTime run_time = RunTime::forMesh(mesh, control);

    require(run_time.loop(), "heat run did not start its time step");
    // 时间层由求解器唯一显式推进：每物理步恰好 save 一次。
    time::History<double> history = time::history(temperature);
    history.save(temperature, run_time.deltaT());
    auto equation = testEquation(temperature);
    equ::ddt(equation, 1.0, history);
    equ::laplacian(equation, 1.0, -1.0);
    const SolveResult result = equ::solve(equation);
    require(result.converged(), "heat equation did not converge");
    require(!run_time.loop(), "heat run used an unexpected number of time steps");
    require(
        near(detail::fieldData(temperature)[0], 10.0 / 22.0, 1e-12),
        "implicit heat equation does not match the one-cell FVM result");
    require(run_time.step() == 1, "runtime did not advance exactly one step");

    // 系数也可以是 cell Field：扩散系数由算子插值到面，ddt 系数保持在单元上；
    // Solver 不需要编写插值、矩阵组装或并行同步。
    ScalarField variable_temperature(mesh, FieldLocation::Cell, "Tv", 1.0);
    ScalarField heat_capacity(mesh, FieldLocation::Cell, "rhoCp", 2.0);
    ScalarField conductivity(mesh, FieldLocation::Cell, "k", 1.0);
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        variable_temperature.boundary(patch) = fixedValue(0.0);
    }
    time::History<double> variable_history = time::history(variable_temperature);
    variable_history.save(variable_temperature, run_time.deltaT());
    auto variable_equation = testEquation(variable_temperature);
    equ::ddt(variable_equation, heat_capacity, variable_history);
    equ::laplacian(variable_equation, conductivity, -1.0);
    const SolveResult variable_result = equ::solve(variable_equation);
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
    time::History<double> field_history = time::history(field_temperature);
    field_history.save(field_temperature, run_time.deltaT());
    auto field_equation = testEquation(field_temperature);
    equ::ddt(field_equation, heat_capacity, field_history);
    equ::laplacian(field_equation, conductivity, -1.0);
    equ::source(field_equation, field_source);
    const SolveResult field_result = equ::solve(field_equation);
    require(field_result.converged(), "Field-material heat step did not converge");
}

void checkTransportEquation() {
    const Mesh mesh = makeHexBox({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
    ScalarField concentration(mesh, FieldLocation::Cell, "C", 0.0);
    ScalarField flux(mesh, FieldLocation::Face, "phi", 0.0);
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        concentration.boundary(patch) = zeroGradient();
    }

    RuntimeControl control;
    control.methods.time = TimeMethod::Euler;
    control.time = {0.0, 0.1, 0.1};
    RunTime run_time = RunTime::forMesh(mesh, control);
    require(run_time.loop(), "transport run did not start its time step");
    time::History<double> history = time::history(concentration);
    history.save(concentration, run_time.deltaT());
    // ddt(2, C) + div(phi, C) == laplacian(0, C) + source(6)
    auto equation = testEquation(concentration);
    equ::ddt(equation, 2.0, history);
    equ::div(equation, flux);
    equ::laplacian(equation, 0.0, -1.0);
    equ::source(equation, 6.0);
    const SolveResult result = equ::solve(equation);
    require(result.converged(), "transport equation did not converge");
    require(!run_time.loop() && run_time.step() == 1, "transport time step count changed");
    require(near(detail::fieldData(concentration)[0], 0.3, 1e-12), "transport source update is incorrect");
}

}  // namespace

int main() {
    checkHeatEquation();
    checkTransportEquation();
    std::cout << "scalar_equation_test: heat coefficients and transport assembly passed\n";
}
