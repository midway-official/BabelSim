#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/solver.h"
#include <iostream>

namespace babelsim {
int runTransport(Case& problem) {
    auto& C = problem.scalarField("C");
    auto& U = problem.vectorField("U");
    auto& phi = problem.faceFlux("phi", U);
    const auto& physical = problem.physics();
    const double storage = physical.positive("storage");
    const double D = physical.nonnegative("diffusivity");
    const double S = physical.number("source");

    loadMethods(problem);
    const auto linearOptions = readLinearControl(problem, C);
    const int writeInterval = readWriteInterval(problem);
    auto time = enableTime(problem);
    auto C_old = math::history(C);
    auto transportEquation = equ::createEquation(C);
    problem.validate();

    while (time.value() < time.end()) {
        advance(time);
        math::saveOld(C_old, C, time.dt());
        equ::reset(transportEquation);
        equ::ddt(transportEquation, storage, C_old);
        equ::div(transportEquation, phi);
        equ::laplacian(transportEquation, D, -1.0);
        equ::source(transportEquation, S);
        const auto result = equ::solve(transportEquation, C, linearOptions);

        if (primaryProcess())
            std::cout << "transport time=" << time.value()
                      << " residual=" << result.relative_residual << '\n';
        if (!result.converged()) return 2;
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return 0;
}
const SolverRegistration transport("transport", runTransport);
}
