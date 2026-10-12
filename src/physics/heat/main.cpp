#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/builtin_descriptions.h"
#include "babelsim/time_discretization.h"
#include "babelsim/solver.h"
#include "babelsim/monitor.h"

namespace babelsim {
SolverResult runHeat(Case& problem) {
    const monitor::Reporter reporter("heat");
    auto& T = problem.scalarField("T");
    const auto parameters = problem.parameters(builtin::heat());
    const auto& physical = parameters.physics;
    const double rho = physical.number("density");
    const double cp = physical.number("heatCapacity");
    const double k = physical.number("conductivity");
    const double Q = physical.number("source");

    const int writeInterval = readWriteInterval(problem);
    auto time = time::start(problem);
    auto T_old = time::history(T);
    time::HistoryGroup histories;
    histories.add(T, T_old);
    auto temperatureEquation = equ::createEquation(problem, "temperature", T, {"diffusion"});
    problem.validate();

    while (time.value() < time.end()) {
        time.advance();
        histories.save(time);
        if (problem.methods().time == TimeMethod::BDF2 && T_old.levels() >= 2)
            time::predict(T, T_old, problem.methods().time);
        temperatureEquation.reset();
        equ::ddt(temperatureEquation, rho * cp, T_old);
        equ::laplacian(temperatureEquation, k, -1, "diffusion");
        equ::source(temperatureEquation, Q);
        const auto result = equ::solve(temperatureEquation);

        reporter.record({{"time", time.value()}, {"residual", result.relative_residual}});
        if (!result.converged()) return SolverResult{result.status};
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return SolverResult::completed();
}
const SolverRegistration heat("heat", runHeat, &builtin::heat());
}
