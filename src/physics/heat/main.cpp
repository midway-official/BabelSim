#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/solver.h"

namespace babelsim {

int runHeat(Case& problem) {
    ScalarField& T = problem.scalarField("T");
    const double rho = problem.physics().positive("density");
    const double cp = problem.physics().positive("heatCapacity");
    const double k = problem.physics().nonnegative("conductivity");
    const double Q = problem.physics().number("source");

    while (problem.loop()) {
        if (!solve(eqn::ddt(rho * cp, T) ==
                   eqn::laplacian(k, T) + eqn::source(Q)).converged()) return 2;
    }
    return 0;
}

const SolverRegistration heat("heat", runHeat);

}  // babelsim 命名空间
