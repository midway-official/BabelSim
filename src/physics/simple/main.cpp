#include "babelsim/application.h"
#include "babelsim/case.h"
#include "algorithm.h"

namespace babelsim {

int runSimple(Case& problem) {
    SteadySimpleAlgorithm simple(problem);
    problem.start();
    while (simple.loop()) {
        simple.solveMomentum();
        simple.solvePressure();
        simple.correctVelocity();
        simple.correctFlux();
        simple.checkContinuity();
    }
    return simple.converged() ? 0 : 2;
}

const SolverRegistration simple("simple", runSimple);

}  // babelsim 命名空间
