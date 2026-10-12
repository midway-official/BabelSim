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
struct CoupledConfig { double diffusivity, coupling, tolerance; int iterations; };
const babelsim::ConfigurationDescription<CoupledConfig>& description() {
    using namespace babelsim;
    static const auto spec = [] {
        ConfigurationDescription<CoupledConfig> result;
        result.name = "coupledScalar";
        result.explanation = "Two scalar diffusion equations with conservative exchange.";
        result.fields = {
            {"T", DescribedFieldType::Scalar, DescribedFieldLocation::Cell, FieldRole::Input, "1", "First scalar.", "1"},
            {"C", DescribedFieldType::Scalar, DescribedFieldLocation::Cell, FieldRole::Input, "1", "Second scalar.", "0"},
            {"total", DescribedFieldType::Scalar, DescribedFieldLocation::Cell, FieldRole::Output, "1", "Conserved sum.", "0", FieldInitialization::Computed, {}}};
        ParameterDescription diffusion;
        diffusion.name = "diffusivity"; diffusion.minimum = 0.0; diffusion.unit = "m^2/s";
        diffusion.explanation = "Common nonnegative diffusion coefficient.";
        result.bind(diffusion, &CoupledConfig::diffusivity);
        auto exchange = diffusion; exchange.name = "coupling"; exchange.unit = "1/s";
        exchange.explanation = "Nonnegative exchange coefficient.";
        result.bind(exchange, &CoupledConfig::coupling);
        ParameterDescription iterations;
        iterations.name = "couplingIterations"; iterations.dictionary = ParameterDictionary::Solution;
        iterations.required = false; iterations.default_value = 100; iterations.minimum = 1; iterations.maximum = 10000;
        iterations.explanation = "Maximum coupled iterations per time step.";
        result.bind(iterations, &CoupledConfig::iterations);
        auto tolerance = iterations; tolerance.name = "couplingTolerance"; tolerance.default_value = 1e-10;
        tolerance.minimum = 0.0; tolerance.exclusive_minimum = true; tolerance.maximum.reset();
        tolerance.explanation = "Relative scalar change acceptance tolerance.";
        result.bind(tolerance, &CoupledConfig::tolerance);
        result.equations = {{"temperature", "T", {"diffusion"}, "First exchange equation."},
            {"concentration", "C", {"diffusion"}, "Second exchange equation."}};
        result.time_methods = {TimeMethod::Euler};
        result.boundary_types = {"fixedValue", "fixedGradient", "zeroGradient", "symmetry", "inletOutlet"};
        return result;
    }();
    return spec;
}

babelsim::SolverResult runCoupled(babelsim::Case& problem) {
    using namespace babelsim;
    const auto& spec = description();
    const auto config = problem.read.config(spec);
    auto& T = problem.load.scalarField("T");
    auto& C = problem.load.scalarField("C");
    auto& total = problem.create.scalarField("total");
    problem.output(total);
    auto temperature = equ::createEquation(problem, spec.equation("temperature"), T);
    auto concentration = equ::createEquation(problem, spec.equation("concentration"), C);
    time::HistoryGroup histories;
    auto& oldT = histories.addHistory(T);
    auto& oldC = histories.addHistory(C);
    histories.trackForRollback(total);
    auto clock = time::makeStepper(problem);
    const int interval = readWriteInterval(problem);
    const monitor::Reporter reporter("coupledScalar");
    problem.finalizeInitialization();
    while (!clock.finished()) {
        auto step = time::beginStep(clock, histories);
        bool settled = false;
        for (int iteration = 0; iteration < config.iterations; ++iteration) {
            const ScalarField previousT = T;
            const ScalarField previousC = C;
            temperature.reset();
            equ::ddt(temperature, 1.0, oldT);
            equ::laplacian(temperature, config.diffusivity, -1, "diffusion");
            equ::reaction(temperature, config.coupling);
            equ::source(temperature, C, config.coupling);
            const auto solveT = equ::solve(temperature);
            if (!solveT.converged()) return SolverResult{solveT.status};
            concentration.reset();
            equ::ddt(concentration, 1.0, oldC);
            equ::laplacian(concentration, config.diffusivity, -1, "diffusion");
            equ::reaction(concentration, config.coupling);
            equ::source(concentration, T, config.coupling);
            const auto solveC = equ::solve(concentration);
            if (!solveC.converged()) return SolverResult{solveC.status};
            const double change = std::max(diagnostics::relativeChange(T, previousT),
                                           diagnostics::relativeChange(C, previousC));
            if (change <= config.tolerance) {
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
const babelsim::SolverRegistration registration("coupledScalar", runCoupled, &description());
} // namespace
