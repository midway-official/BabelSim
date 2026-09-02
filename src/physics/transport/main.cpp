#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/solver.h"

namespace babelsim {

int runTransport(Case& problem) {
    ScalarField& C = problem.scalarField("C");
    VectorField& U = problem.vectorField("U");
    ScalarField& phi = problem.faceFlux("phi", U);
    const double storage = problem.physics().positive("storage");
    const double D = problem.physics().nonnegative("diffusivity");
    const double S = problem.physics().number("source");

    while (problem.loop()) {
        if (!solve(eqn::ddt(storage, C) + eqn::div(phi, C) ==
                   eqn::laplacian(D, C) + eqn::source(S)).converged()) return 2;
    }
    return 0;
}

const SolverRegistration transport("transport", runTransport);

}  // babelsim 命名空间
