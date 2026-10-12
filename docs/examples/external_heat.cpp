#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/builtin_descriptions.h"
#include "babelsim/time_discretization.h"
#include "babelsim/time_step.h"
#include "babelsim/monitor.h"

namespace {
const babelsim::ModuleDescription& description() {
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
    auto& T = problem.scalarField("T");
    const auto parameters = problem.parameters(description());
    const auto& physical = parameters.physics;
    const double rho = physical.number("density");
    const double cp = physical.number("heatCapacity");
    const double k = physical.number("conductivity");
    const double Q = physical.number("source");
    auto equation = equ::createEquation(problem, "temperature", T, {"diffusion"});
    auto clock = time::start(problem);
    auto history = time::history(T);
    time::HistoryGroup histories;
    histories.add(T, history);
    const int interval = readWriteInterval(problem);
    const monitor::Reporter reporter("customHeat");
    problem.validate();

    while (!clock.finished()) {
        auto step = time::beginStep(clock, histories);
        if (problem.methods().time == TimeMethod::BDF2 && history.levels() >= 2)
            time::predict(T, history, problem.methods().time);
        equation.reset();
        equ::ddt(equation, rho * cp, history);
        equ::laplacian(equation, k, -1, "diffusion");
        equ::source(equation, Q);
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
