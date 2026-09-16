#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/solver.h"
#include <iostream>

namespace babelsim {
int runHeat(Case& problem) {
    auto& T = problem.scalarField("T");
    const auto& physical = problem.physics();
    const double rho = physical.positive("density");
    const double cp = physical.positive("heatCapacity");
    const double k = physical.nonnegative("conductivity");
    const double Q = physical.number("source");

    loadMethods(problem);
    const auto linearOptions = readLinearControl(problem, T);
    const int writeInterval = readWriteInterval(problem);
    auto time = enableTime(problem);
    auto T_old = math::history(T);
    auto temperatureEquation = equ::createEquation(T);
    problem.output(T);
    problem.validate();

    while (time.value() < time.end()) {
        advance(time);
        math::saveOld(T_old, T, time.dt());
        equ::reset(temperatureEquation);
        equ::ddt(temperatureEquation, rho * cp, T_old);
        equ::laplacian(temperatureEquation, k, -1.0);
        equ::source(temperatureEquation, Q);
        const auto result = equ::solve(temperatureEquation, T, linearOptions);

        if (primaryProcess())
            std::cout << "heat time=" << time.value()
                      << " residual=" << result.relative_residual << '\n';
        if (!result.converged()) return 2;
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return 0;
}
const SolverRegistration heat("heat", runHeat);
}
