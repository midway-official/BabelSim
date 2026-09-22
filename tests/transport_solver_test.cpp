#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/equ.h"
#include "babelsim/runtime.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
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
    std::cout << "transport_solver_test: C=" << detail::fieldData(concentration)[0] << '\n';
}
