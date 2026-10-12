#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/builtin_descriptions.h"
#include "babelsim/time_discretization.h"
#include "babelsim/time_step.h"
#include "babelsim/monitor.h"

namespace {
const babelsim::ConfigurationDescription<babelsim::builtin::HeatConfig>& description() {
    static const auto module = [] {
        auto result = babelsim::builtin::heat();
        result.name = "customHeat";
        result.explanation = "External heat solver using public time/history tools.";
        return result;
    }();
    return module;
}
babelsim::SolverResult runCustomHeat(babelsim::Case& problem) {
    using namespace babelsim;
    const auto& spec = description();
    const auto config = problem.read.config(spec);
    const auto method = problem.read.methods().time;
    auto& T = problem.load.scalarField("T");
    auto equation = equ::createEquation(problem, spec.equation("temperature"), T);
    time::HistoryGroup histories;
    auto& history = histories.addHistory(T);
    auto clock = time::makeStepper(problem);
    const int interval = readWriteInterval(problem);
    const monitor::Reporter reporter("customHeat");
    problem.finalizeInitialization();

    while (!clock.finished()) {
        auto step = time::beginStep(clock, histories);
        if (method == TimeMethod::BDF2 && history.levels() >= 2)
            time::predict(T, history, method);
        equation.reset();
        equ::ddt(equation, config.density * config.heatCapacity, history);
        equ::laplacian(equation, config.conductivity, -1, "diffusion");
        equ::source(equation, config.source);
        const auto result = equ::solve(equation);
        reporter.record({{"time", clock.value()}, {"residual", result.relative_residual}});
        if (!result.converged()) return SolverResult{result.status};
        step.accept();
        if (clock.step() % interval == 0 || clock.finished()) write(problem, clock);
    }
    return SolverResult::completed();
}

const babelsim::SolverRegistration registration("customHeat", runCustomHeat, &description());
} // namespace
