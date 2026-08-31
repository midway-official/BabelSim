#include "babelsim/case.h"
#include "babelsim/simple.h"

namespace babelsim {

int runSimple(Case& problem) {
    SimpleSolver simple(problem);
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

}  // babelsim 命名空间
