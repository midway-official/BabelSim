#include "babelsim/application.h"
#include "babelsim/case.h"
#include "algorithm.h"

namespace babelsim {

int runTransientSimple(Case& problem) {
    TransientSimpleAlgorithm simple(problem);
    while (problem.loop()) {
        simple.beginTimeStep();
        while (simple.loop()) {
            simple.solveMomentum();
            simple.solvePressure();
            simple.correctVelocity();
            simple.correctFlux();
            simple.correctTurbulence();
            simple.checkContinuity();
        }
        if (!simple.converged()) return 2;
    }
    return 0;
}

const SolverRegistration transient_simple("transientSimple", runTransientSimple);

}  // babelsim 命名空间
