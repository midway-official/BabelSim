#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/builtin_descriptions.h"
#include "babelsim/geometry.h"
#include "babelsim/solver.h"
#include "../RANS/api.h"
#include "babelsim/monitor.h"

namespace babelsim {
SolverResult runTransientSimple(Case& problem) {
    const monitor::Reporter reporter("Transient SIMPLE");
    const auto& spec = builtin::transientSimple();
    const auto config = problem.read.config(spec);

    auto& U = problem.load.vectorField("U");
    auto& p = problem.load.scalarField("p");
    auto& phi = problem.create.faceScalarField("phi");
    auto& pPrime = problem.create.homogeneousCorrection(p, "pPrime");
    auto momentumEquation = equ::createEquation(problem, spec.equation("momentum"), U);
    auto pressureCorrectionEquation = equ::createEquation(problem, spec.equation("pressureCorrection"), pPrime);
    const auto pressureOptions = pressureCorrectionEquation.options("diffusion");
    const int pressureSolves = pressureOptions.diffusion == DiffusionMethod::Orthogonal
        ? 1 : config.coupling.nonOrthogonalCorrections + 1;

    phi = math::flux(U, momentumEquation.options());
    auto turbulence = rans::create(problem, U, phi, config.fluid, momentumEquation.options());
    const auto& muEff = turbulence.effectiveViscosity();
    const auto V = geometry::cellVolumes(problem.mesh());
    const auto Sf = geometry::faceAreaVectors(problem.mesh());
    const auto Af = geometry::faceAreas(problem.mesh());
    const int writeInterval = readWriteInterval(problem);
    time::HistoryGroup histories;
    auto& U_old = histories.addHistory(U);
    histories.trackForRollback(p, phi);
    turbulence.registerTimeState(histories);
    auto time = time::makeStepper(problem);
    problem.finalizeInitialization();

    while (time.value() < time.end()) {
        time.advance();
        histories.save(time);
        turbulence.beginStep();
        bool converged = false;
        for (int iter = 0; iter < config.coupling.maxIterations; ++iter) {
            const VectorField U_previous = U; // previous SIMPLE iterate, not physical time history

            // Assemble the momentum equation.
            momentumEquation.reset();
            equ::ddt(momentumEquation, config.fluid.density, U_old);
            equ::div(momentumEquation, phi, config.fluid.density, "convection");
            equ::laplacian(momentumEquation, muEff, -1, "diffusion");
            equ::source(momentumEquation, -math::grad(p, pressureOptions));
            if (turbulence)
            equ::source(momentumEquation, math::div(
                turbulence.deviatoricStressRemainder(U), momentumEquation.options()));

            // Measure the current outer iterate before relaxation or prediction.
            const double rU = diagnostics::relativeResidual(momentumEquation, U);
            if (!diagnostics::all(std::isfinite(rU))) return SolverResult::numericalFailure();
            equ::relax(momentumEquation, U_previous, config.coupling.velocityRelaxation);
            // Preserve SIMPLE row normalization; V/aP below uses these scaled rows.
            equ::scale(momentumEquation, config.coupling.velocityRelaxation);
            const auto aP = momentumEquation.diagonal();
            const auto rAU = V / aP;
            const auto velocitySolve = equ::solve(momentumEquation);
            if (!diagnostics::all(velocitySolve.healthy())) return SolverResult::numericalFailure();

            // Rhie-Chow interpolation; physical boundary fluxes remain unchanged.
            const auto gradP = math::grad(p, pressureOptions);
            auto phiH = math::flux(U, momentumEquation.options());
            const auto interpolatedGradientFlux =
                math::dot(math::interpolate(rAU * gradP, pressureOptions), Sf);
            const auto normalGradientFlux =
                math::interpolate(rAU, pressureOptions.coefficientOptions()) * math::normalGradient(p, gradP, pressureOptions) * Af;
            math::add(interpolatedGradientFlux, phiH, math::FaceRegion::Interior);
            math::subtract(normalGradientFlux, phiH, math::FaceRegion::Interior);
            const auto divPhiH = math::div(phiH);

            pPrime.fill(0.0);
            bool pressureConverged = true;
            double pressureLinearResidual = 0.0;
            for (int correction = 0; correction < pressureSolves; ++correction) {
                pressureCorrectionEquation.reset();
                equ::laplacian(pressureCorrectionEquation, rAU, -1, "diffusion");
                equ::source(pressureCorrectionEquation, -divPhiH);
                pressureCorrectionEquation.referenceIfUnanchored(0.0);
                const auto pressureSolve = equ::solve(pressureCorrectionEquation);
                if (!diagnostics::all(pressureSolve.healthy())) return SolverResult::numericalFailure();
                pressureConverged = pressureConverged && pressureSolve.converged();
                pressureLinearResidual = pressureSolve.relative_residual;
            }

            p += config.coupling.pressureRelaxation * pPrime;
            auto pressureFlux = equ::faceFlux(pressureCorrectionEquation, pPrime);
            // Preserve this transient algorithm's coupled pressure/velocity relaxation.
            pPrime *= config.coupling.pressureRelaxation;
            pressureFlux *= config.coupling.pressureRelaxation;
            U -= rAU * math::grad(pPrime, pressureOptions);
            phi = phiH + pressureFlux;

            bool turbulenceConverged = true;
            double dTurbulence = 0.0, rTurbulence = 0.0;
            if (turbulence) {
                const auto result = turbulence.solveTransport();
                dTurbulence = result.relativeChange();
                rTurbulence = result.initialResidual();
                if (!diagnostics::all(result.healthy() && std::isfinite(dTurbulence) && std::isfinite(rTurbulence))) return SolverResult::numericalFailure();
                turbulenceConverged = result.linearConverged()
                    && dTurbulence <= turbulence.tolerance()
                    && rTurbulence <= turbulence.tolerance();
            }

            const double dU = diagnostics::relativeChange(U, U_previous);
            const double dP = diagnostics::relativeMagnitude(pPrime, p);
            const auto mass = diagnostics::fluxBalance(phi);
            if (!diagnostics::all(std::isfinite(rU) && std::isfinite(dU)
                && std::isfinite(dP) && std::isfinite(mass.relative))) return SolverResult::numericalFailure();
            const bool linearConverged = velocitySolve.converged() && pressureConverged;
            converged = diagnostics::all(linearConverged && turbulenceConverged
                && rU <= config.convergence.momentum && dU <= config.convergence.velocity
                && dP <= config.convergence.pressureCorrection && mass.relative <= config.convergence.continuity);

            reporter.iteration(iter + 1, config.coupling.maxIterations, converged, {
                {"mass", mass.relative}, {"dU", dU}, {"rU", rU}, {"dP", dP},
                {"linP", pressureLinearResidual}, {"dTurb", dTurbulence}, {"rTurb", rTurbulence},
                {"linear", linearConverged ? "ok" : "inexact"}, {"converged", converged}});
            if (converged) break;
        }
        if (!converged) return SolverResult::notConverged();
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return SolverResult::completed();
}
const SolverRegistration transient_simple("transientSimple", runTransientSimple, &builtin::transientSimple());
}
