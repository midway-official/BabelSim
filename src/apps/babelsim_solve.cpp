#include "babelsim/application.h"

namespace babelsim {
int runHeat(Case&);
int runSimple(Case&);
int runTransport(Case&);
}

int main(int argc, char* argv[]) {
    const babelsim::SolverEntry solvers[] = {
        {"heat", babelsim::runHeat},
        {"simple", babelsim::runSimple},
        {"transport", babelsim::runTransport},
    };
    return babelsim::runApplication(argc, argv, solvers);
}
