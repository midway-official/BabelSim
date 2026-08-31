#include "babelsim/case.h"
#include "babelsim/solver.h"

namespace babelsim {

int runHeat(Case& problem) {
    ScalarField& T = problem.scalarField("T");
    const double rho = problem.properties().positive("density");
    const double cp = problem.properties().positive("heatCapacity");
    const double k = problem.properties().nonnegative("conductivity");
    const double Q = problem.properties().number("source");

    while (problem.loop()) {
        if (!solve(fvm::ddt(rho * cp, T) ==
                   fvm::laplacian(k, T) + fvm::source(Q)).converged()) return 2;
    }
    return 0;
}

}  // babelsim 命名空间
