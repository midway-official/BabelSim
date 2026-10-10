#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/time.h"

using namespace babelsim;

SolverResult solveHeat(Case& problem) {
    auto& T = problem.scalarField("T");
    const double rhoCp = problem.physics().positive("density") *
                         problem.physics().positive("heatCapacity");
    const double conductivity = problem.physics().nonnegative("conductivity");
    const double source = problem.physics().number("source");
    const auto writeInterval = readWriteInterval(problem);
    const auto control = readEquationControl(problem, "temperature", T, {"diffusion"});
    auto equation = equ::createEquation(T, control);
    auto stepper = time::start(problem);
    auto history = time::history(T);

    while (!stepper.finished()) {
        stepper.advance();
        history.save(T, stepper.dt());
        equation.reset();
        equ::ddt(equation, rhoCp, history);
        equ::laplacian(equation, conductivity, -1.0, "diffusion");
        equ::source(equation, source);
        const auto result = equ::solve(equation);
        if (!result.converged()) return SolverResult::notConverged();
        if (stepper.step() % writeInterval == 0 || stepper.finished())
            write(problem, stepper);
    }
    return SolverResult::completed();
}

const SolverRegistration registration("myHeat", solveHeat);

int main(int argc, char** argv) { return runApplication(argc, argv); }
