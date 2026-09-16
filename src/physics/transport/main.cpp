#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/solver.h"
#include "babelsim/monitor.h"

namespace babelsim {
SolverResult runTransport(Case& problem) {
    const monitor::Reporter reporter("transport");
    auto& C = problem.scalarField("C");
    auto& U = problem.vectorField("U");
    auto& phi = problem.createFaceField("phi");
    phi = math::flux(U);
    const auto& physical = problem.physics();
    const double storage = physical.positive("storage");
    const double D = physical.nonnegative("diffusivity");
    const double S = physical.number("source");

    const auto linearOptions = readLinearControl(problem, C);
    const int writeInterval = readWriteInterval(problem);
    auto time = time::start(problem);
    auto C_old = time::history(C);
    auto transportEquation = equ::createEquation(C);
    problem.validate();

    while (time.value() < time.end()) {
        time.advance();
        C_old.save(C, time.dt());
        transportEquation.reset();
        equ::ddt(transportEquation, storage, C_old);
        equ::div(transportEquation, phi);
        equ::laplacian(transportEquation, D, -1);
        equ::source(transportEquation, S);
        const auto result = equ::solve(transportEquation, C, linearOptions);

        reporter.record({{"time", time.value()}, {"residual", result.relative_residual}});
        if (!result.converged()) return SolverResult{result.status};
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return SolverResult::completed();
}
const SolverRegistration transport("transport", runTransport);
}
