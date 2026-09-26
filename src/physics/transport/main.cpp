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
    auto& phi = problem.createFaceScalarField("phi");
    const auto& physical = problem.physics();
    const double storage = physical.positive("storage");
    const double D = physical.nonnegative("diffusivity");
    const double S = physical.number("source");

    const int writeInterval = readWriteInterval(problem);
    auto time = time::start(problem);
    auto C_old = time::history(C);
    auto transportEquation = equ::createEquation(problem, "transport", C, {"convection", "diffusion"});
    phi = math::flux(U, transportEquation.options());
    problem.validate();

    while (time.value() < time.end()) {
        time.advance();
        C_old.save(C, time.dt());
        // BDF2 also requires a second-order time value for deferred convection
        // and diffusion corrections; the first step retains Euler startup.
        if (problem.methods().time == TimeMethod::BDF2 && C_old.levels() >= 2) {
            const double ratio = C_old.dt() / C_old.previousDt();
            C.assignScaled(1.0 + ratio, C_old.previous());
            C.addScaled(-ratio, C_old.older());
        }
        transportEquation.reset();
        equ::ddt(transportEquation, storage, C_old);
        equ::div(transportEquation, phi, 1.0, "convection");
        equ::laplacian(transportEquation, D, -1, "diffusion");
        equ::source(transportEquation, S);
        const auto result = equ::solve(transportEquation);

        reporter.record({{"time", time.value()}, {"residual", result.relative_residual}});
        if (!result.converged()) return SolverResult{result.status};
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return SolverResult::completed();
}
const SolverRegistration transport("transport", runTransport);
}
