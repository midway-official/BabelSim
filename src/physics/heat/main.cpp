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
    const auto& spec = builtin::heat();
    const auto config = problem.read.config(spec);
    const auto method = problem.read.methods().time;
    const int writeInterval = readWriteInterval(problem);

    auto& T = problem.load.scalarField("T");
    auto temperatureEquation = equ::createEquation(problem, spec.equation("temperature"), T);
    time::HistoryGroup histories;
    auto& T_old = histories.addHistory(T);
    auto time = time::makeStepper(problem);
    problem.finalizeInitialization();

    while (time.value() < time.end()) {
        time.advance();
        histories.save(time);
        if (method == TimeMethod::BDF2 && T_old.levels() >= 2)
            time::predict(T, T_old, method);
        temperatureEquation.reset();
        equ::ddt(temperatureEquation, config.density * config.heatCapacity, T_old);
        equ::laplacian(temperatureEquation, config.conductivity, -1, "diffusion");
        equ::source(temperatureEquation, config.source);
        const auto result = equ::solve(temperatureEquation);

        reporter.record({{"time", time.value()}, {"residual", result.relative_residual}});
        if (!result.converged()) return SolverResult{result.status};
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return SolverResult::completed();
}
const SolverRegistration heat("heat", runHeat, &builtin::heat());
}
