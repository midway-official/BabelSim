#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/runtime.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    const Mesh mesh = Mesh::cartesian({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
    ScalarField concentration(mesh, FieldLocation::Cell, "C", 0.0);
    ScalarField flux(mesh, FieldLocation::Face, "phi", 0.0);
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        concentration.boundary(patch) = zeroGradient();
    }

    RuntimeControl control;
    control.methods.time = TimeMethod::Euler;
    control.time = {0.0, 0.1, 0.1};
    control.scalar_solver.absolute_tolerance = 1e-14;
    RunTime run_time = RunTime::forMesh(mesh, control);
    require(run_time.loop(), "transport run did not start its time step");
    const SolveResult result = solve(
        eqn::ddt(2.0, concentration) + eqn::div(flux, concentration) ==
            eqn::laplacian(0.0, concentration) + eqn::source(6.0));
    require(result.converged(), "transport equation did not converge");
    require(!run_time.loop() && run_time.step() == 1, "transport time step count changed");
    require(near(detail::fieldData(concentration)[0], 0.3, 1e-12), "transport source update is incorrect");
    std::cout << "transport_solver_test: C=" << detail::fieldData(concentration)[0] << '\n';
}
