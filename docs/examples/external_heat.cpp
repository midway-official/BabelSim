#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/monitor.h"

namespace {
babelsim::SolverResult runCustomHeat(babelsim::Case& problem) {
    using namespace babelsim;
    auto& T = problem.scalarField("T");
    const auto& physical = problem.physics();
    const double rho = physical.positive("density");
    const double cp = physical.positive("heatCapacity");
    const double k = physical.nonnegative("conductivity");
    const double Q = physical.number("source");
    auto equation = equ::createEquation(problem, "temperature", T, {"diffusion"});
    auto clock = time::start(problem);
    auto history = time::history(T);
    const int interval = readWriteInterval(problem);
    const monitor::Reporter reporter("customHeat");
    problem.validate();

    while (!clock.finished()) {
        clock.advance();
        history.save(T, clock.dt());
        if (problem.methods().time == TimeMethod::BDF2 && history.levels() >= 2) {
            const double ratio = history.dt() / history.previousDt();
            T.assignScaled(1.0 + ratio, history.previous());
            T.addScaled(-ratio, history.older());
        }
        equation.reset();
        equ::ddt(equation, rho * cp, history);
        equ::laplacian(equation, k, -1, "diffusion");
        equ::source(equation, Q);
        const auto result = equ::solve(equation);
        reporter.record({{"time", clock.value()}, {"residual", result.relative_residual}});
        if (!result.converged()) return SolverResult{result.status};
        if (clock.step() % interval == 0 || clock.finished()) write(problem, clock);
    }
    return SolverResult::completed();
}

const babelsim::SolverRegistration registration("customHeat", runCustomHeat);
} // namespace
