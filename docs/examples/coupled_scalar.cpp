#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/math.h"
#include "babelsim/monitor.h"
#include "babelsim/solver.h"
#include "babelsim/time_step.h"
#include <algorithm>
#include <stdexcept>

namespace {
babelsim::SolverResult runCoupled(babelsim::Case& problem) {
    using namespace babelsim;
    if (problem.methods().time != TimeMethod::Euler)
        throw std::invalid_argument("coupledScalar tutorial requires time euler");
    auto& T = problem.scalarField("T");
    auto& C = problem.scalarField("C");
    auto& total = problem.createScalarField("total");
    problem.output(total);
    const double D = problem.physics().nonnegative("diffusivity");
    const double a = problem.physics().nonnegative("coupling");
    const int limit = problem.solution().integer("couplingIterations", 100, 1, 10000);
    const double tolerance = problem.solution().positive("couplingTolerance", 1e-10);
    auto temperature = equ::createEquation(problem, "temperature", T, {"diffusion"});
    auto concentration = equ::createEquation(problem, "concentration", C, {"diffusion"});
    auto oldT = time::history(T);
    auto oldC = time::history(C);
    time::HistoryGroup histories;
    histories.add(T, oldT);
    histories.add(C, oldC);
    histories.track(total);
    auto clock = time::start(problem);
    const int interval = readWriteInterval(problem);
    const monitor::Reporter reporter("coupledScalar");
    problem.validate();
    while (!clock.finished()) {
        auto step = time::beginStep(clock, histories);
        bool settled = false;
        for (int iteration = 0; iteration < limit; ++iteration) {
            const ScalarField previousT = T;
            const ScalarField previousC = C;
            temperature.reset();
            equ::ddt(temperature, 1.0, oldT);
            equ::laplacian(temperature, D, -1, "diffusion");
            equ::reaction(temperature, a);
            equ::source(temperature, C, a);
            const auto solveT = equ::solve(temperature);
            if (!solveT.converged()) return SolverResult{solveT.status};
            concentration.reset();
            equ::ddt(concentration, 1.0, oldC);
            equ::laplacian(concentration, D, -1, "diffusion");
            equ::reaction(concentration, a);
            equ::source(concentration, T, a);
            const auto solveC = equ::solve(concentration);
            if (!solveC.converged()) return SolverResult{solveC.status};
            const double change = std::max(diagnostics::relativeChange(T, previousT),
                                           diagnostics::relativeChange(C, previousC));
            if (change <= tolerance) {
                settled = true;
                break;
            }
        }
        reporter.record({{"time", clock.value()}, {"settled", settled ? 1.0 : 0.0}});
        if (!settled) return SolverResult::notConverged();
        total.assign(T + C);
        step.accept();
        if (clock.step() % interval == 0 || clock.finished()) write(problem, clock);
    }
    return SolverResult::completed();
}
const babelsim::SolverRegistration registration("coupledScalar", runCoupled);
} // namespace
