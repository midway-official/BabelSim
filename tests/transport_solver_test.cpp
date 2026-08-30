#include "babelsim/transport.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    const Mesh mesh = Mesh::cartesian({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
    ScalarField concentration(mesh, FieldLocation::Cell, "C", 0.0);
    ScalarField flux(mesh, FieldLocation::Face, "phi", 0.0);
    for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        concentration.boundary(patch) = zeroGradient();
    }

    RuntimeControl control;
    control.methods.time = TimeMethod::Euler;
    control.time = {0.0, 0.1, 0.1};
    control.scalar_solver.absolute_tolerance = 1e-14;
    RunTime run_time = RunTime::forMesh(mesh, control);
    const ScalarTransportResult result = solveTransientScalarTransport(
        run_time, concentration, flux, 2.0, 0.0, 6.0);

    require(result.converged && result.steps == 1, "transport equation did not converge");
    require(near(concentration[0], 0.3, 1e-12), "transport source update is incorrect");
    std::cout << "transport_solver_test: C=" << concentration[0] << '\n';
}
