#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/builtin_descriptions.h"
#include "babelsim/time_discretization.h"
#include "babelsim/solver.h"
#include "babelsim/monitor.h"

namespace babelsim {
SolverResult runTransport(Case& problem) {
    const monitor::Reporter reporter("transport");
    const auto& spec = builtin::transport();
    const auto config = problem.read.config(spec);
    const auto method = problem.read.methods().time;
    const int writeInterval = readWriteInterval(problem);

    auto& C = problem.load.scalarField("C");
    auto& U = problem.load.vectorField("U");
    auto& phi = problem.create.faceScalarField("phi");
    auto transportEquation = equ::createEquation(problem, spec.equation("transport"), C);
    phi = math::flux(U, transportEquation.options());
    time::HistoryGroup histories;
    auto& C_old = histories.addHistory(C);
    auto time = time::makeStepper(problem);
    problem.finalizeInitialization();

    while (time.value() < time.end()) {
        time.advance();
        histories.save(time);
        if (method == TimeMethod::BDF2 && C_old.levels() >= 2)
            time::predict(C, C_old, method);
        transportEquation.reset();
        equ::ddt(transportEquation, config.storage, C_old);
        equ::div(transportEquation, phi, 1.0, "convection");
        equ::laplacian(transportEquation, config.diffusivity, -1, "diffusion");
        equ::source(transportEquation, config.source);
        const auto result = equ::solve(transportEquation);

        reporter.record({{"time", time.value()}, {"residual", result.relative_residual}});
        if (!result.converged()) return SolverResult{result.status};
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return SolverResult::completed();
}
const SolverRegistration transport("transport", runTransport, &builtin::transport());
}
