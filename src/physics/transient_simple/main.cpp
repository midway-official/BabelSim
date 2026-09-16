#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/geometry.h"
#include "babelsim/solver.h"
#include "../RANS/api.h"
#include "babelsim/monitor.h"
#include <limits>

namespace babelsim {
SolverResult runTransientSimple(Case& problem) {
    const monitor::Reporter reporter("Transient SIMPLE");
    auto& U = problem.vectorField("U");
    auto& p = problem.scalarField("p");
    auto& phi = problem.createFaceScalarField("phi");
    const auto& physical = problem.physics();
    const double rho = physical.positive("density");

    const auto& methods = problem.methods(); // loaded once by Case/runtime
    const auto V = geometry::cellVolumes(problem.mesh());
    const auto Sf = geometry::faceAreaVectors(problem.mesh());
    const auto Af = geometry::faceAreas(problem.mesh());
    const auto velocitySolver = readLinearControl(problem, U);
    const auto pressureSolver = readLinearControl(problem, p);
    const auto& settings = problem.solution();
    const int maxIterations = settings.integer("maxIterations", 1000, 1, std::numeric_limits<int>::max());
    const int nonOrthogonalCorrections = settings.integer("nonOrthogonalCorrections", 1, 0, 20);
    const double alphaU = settings.fraction("velocityRelaxation", 0.7);
    const double alphaP = settings.fraction("pressureRelaxation", 0.3);
    const double massTolerance = settings.positive("continuityTolerance", 1e-8);
    const double velocityTolerance = settings.positive("velocityTolerance", 1e-7);
    const double momentumTolerance = settings.positive("momentumTolerance", 1e-6);
    const double pressureTolerance = settings.positive("pressureCorrectionTolerance", 1e-6);

    auto turbulence = rans::load(problem, U, phi);
    const auto& muEff = turbulence.effectiveViscosity();
    auto pPrime = field::homogeneousLike(p, "pPrime");
    const int pressureSolves = methods.diffusionFor(pPrime.name()) == DiffusionMethod::Orthogonal
        ? 1 : nonOrthogonalCorrections + 1;
    auto momentumEquation = equ::createEquation(U);
    auto pressureCorrectionEquation = equ::createEquation(pPrime);
    phi = math::flux(U);
    if (methods.time == TimeMethod::Steady)
        throw std::invalid_argument("transient SIMPLE requires a transient time scheme");
    auto time = time::start(problem);
    const int writeInterval = readWriteInterval(problem);
    auto U_old = time::history(U);
    problem.validate();

    while (time.value() < time.end()) {
        time.advance();
        U_old.save(U, time.dt());
        if (turbulence) turbulence.saveOld(time.dt());
        bool converged = false;
        for (int iter = 0; iter < maxIterations; ++iter) {
            const VectorField U_previous = U; // previous SIMPLE iterate, not physical time history

            // Assemble the momentum equation.
            momentumEquation.reset();
            equ::ddt(momentumEquation, rho, U_old);
            equ::div(momentumEquation, phi, rho);
            equ::laplacian(momentumEquation, muEff, -1);
            equ::source(momentumEquation, -math::grad(p));
            if (turbulence)
                equ::source(momentumEquation, math::div(turbulence.deviatoricStressRemainder(U)));

            // Measure the current outer iterate before relaxation or prediction.
            const double rU = diagnostics::relativeResidual(momentumEquation, U);
            if (!diagnostics::all(std::isfinite(rU))) return SolverResult::numericalFailure();
            equ::relax(momentumEquation, U_previous, alphaU);
            // Preserve SIMPLE row normalization; V/aP below uses these scaled rows.
            equ::scale(momentumEquation, alphaU);
            const auto aP = momentumEquation.diagonal();
            const auto rAU = V / aP;
            const auto velocitySolve = equ::solve(momentumEquation, velocitySolver);
            if (!diagnostics::all(velocitySolve.healthy())) return SolverResult::numericalFailure();

            // Rhie-Chow interpolation; physical boundary fluxes remain unchanged.
            const auto gradP = math::grad(p);
            auto phiH = math::flux(U);
            const auto interpolatedGradientFlux =
                math::dot(math::interpolate(rAU * gradP), Sf);
            const auto normalGradientFlux =
                math::interpolate(rAU) * math::normalGradient(p, gradP) * Af;
            math::add(interpolatedGradientFlux, phiH, math::FaceRegion::Interior);
            math::subtract(normalGradientFlux, phiH, math::FaceRegion::Interior);
            const auto divPhiH = math::div(phiH);

            pPrime.fill(0.0);
            bool pressureConverged = true;
            double pressureLinearResidual = 0.0;
            for (int correction = 0; correction < pressureSolves; ++correction) {
                pressureCorrectionEquation.reset();
                equ::laplacian(pressureCorrectionEquation, rAU, -1);
                equ::source(pressureCorrectionEquation, -divPhiH);
                pressureCorrectionEquation.referenceIfUnanchored(0.0);
                const auto pressureSolve = equ::solve(pressureCorrectionEquation, pressureSolver);
                if (!diagnostics::all(pressureSolve.healthy())) return SolverResult::numericalFailure();
                pressureConverged = pressureConverged && pressureSolve.converged();
                pressureLinearResidual = pressureSolve.relative_residual;
            }

            p += alphaP * pPrime;
            auto pressureFlux = equ::faceFlux(pressureCorrectionEquation, pPrime);
            // Preserve this transient algorithm's coupled pressure/velocity relaxation.
            pPrime *= alphaP;
            pressureFlux *= alphaP;
            U -= rAU * math::grad(pPrime);
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
                && rU <= momentumTolerance && dU <= velocityTolerance
                && dP <= pressureTolerance && mass.relative <= massTolerance);

            reporter.iteration(iter + 1, maxIterations, converged, {
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
const SolverRegistration transient_simple("transientSimple", runTransientSimple);
}
