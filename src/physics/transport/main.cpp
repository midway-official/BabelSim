#include "babelsim/case.h"
#include "babelsim/solver.h"

namespace babelsim {

int runTransport(Case& problem) {
    ScalarField& C = problem.scalarField("C");
    VectorField& U = problem.vectorField("U");
    ScalarField& phi = problem.faceFlux("phi", U);
    const double storage = problem.properties().positive("storage");
    const double D = problem.properties().nonnegative("diffusivity");
    const double S = problem.properties().number("source");

    while (problem.loop()) {
        if (!solve(fvm::ddt(storage, C) + fvm::div(phi, C) ==
                   fvm::laplacian(D, C) + fvm::source(S)).converged()) return 2;
    }
    return 0;
}

}  // babelsim 命名空间
