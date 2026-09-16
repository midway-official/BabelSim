#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/solver.h"
#include "../RANS/api.h"
#include <iostream>
#include <limits>

namespace babelsim {
int runSimple(Case& problem) {
    auto& U = problem.vectorField("U");
    auto& p = problem.scalarField("p");
    auto& phi = problem.faceField("phi");
    const auto& physical = problem.physics();
    const double rho = physical.positive("density");

    const auto& methods = loadMethods(problem);
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
    const auto& muEff = turbulence.viscosity();
    auto pPrime = math::createHomogeneousField(p);
    const int pressureSolves = methods.diffusionFor(pPrime.name()) == DiffusionMethod::Orthogonal
        ? 1 : nonOrthogonalCorrections + 1;
    auto momentumEquation = equ::createEquation(U);
    auto pressureCorrectionEquation = equ::createEquation(pPrime);
    phi = math::flux(U);
    if (methods.time != TimeMethod::Steady)
        throw std::invalid_argument("steady SIMPLE requires a steady time scheme");
    problem.validate();

    bool converged = false;
    for (int iter = 0; iter < maxIterations; ++iter) {
        const auto U_previous = math::copy(U);

        // Assemble the momentum equation.
        equ::reset(momentumEquation);
        equ::div(momentumEquation, phi, rho);
        equ::laplacian(momentumEquation, muEff, -1.0);
        equ::source(momentumEquation, -math::grad(p));
        if (turbulence) {
            const auto gradU = math::grad(U);
            const auto stress = muEff * (math::transpose(gradU)
                - (2.0 / 3.0) * math::isotropic(math::trace(gradU)));
            equ::source(momentumEquation, math::div(stress));
        }

        // Measure the current outer iterate before relaxation or prediction.
        const double rU = equ::relativeResidual(momentumEquation, U);
        if (!diagnostics::all(std::isfinite(rU))) return 2;
        equ::relax(momentumEquation, U_previous, alphaU);
        equ::scale(momentumEquation, alphaU); // retain the established SIMPLE row normalization
        const auto rAU = equ::response(momentumEquation);
        const auto velocitySolve = equ::solve(momentumEquation, U, velocitySolver);
        if (!diagnostics::all(velocitySolve.healthy())) return 2;

        // Rhie-Chow interpolation; physical boundary fluxes remain unchanged.
        const auto gradP = math::grad(p);
        auto phiH = math::flux(U);
        math::add(math::flux(math::interpolate(rAU * gradP)), phiH, math::FaceRegion::Interior);
        math::subtract(math::flux(math::interpolate(rAU), p, gradP), phiH, math::FaceRegion::Interior);
        const auto divPhiH = math::div(phiH);

        pPrime.fill(0.0);
        bool pressureConverged = true;
        double pressureLinearResidual = 0.0;
        for (int correction = 0; correction < pressureSolves; ++correction) {
            equ::reset(pressureCorrectionEquation);
            equ::laplacian(pressureCorrectionEquation, rAU, -1.0);
            equ::source(pressureCorrectionEquation, -divPhiH);
            equ::reference(pressureCorrectionEquation, 0.0);
            const auto pressureSolve = equ::solve(pressureCorrectionEquation, pPrime, pressureSolver);
            if (!diagnostics::all(pressureSolve.healthy())) return 2;
            pressureConverged = pressureConverged && pressureSolve.converged();
            pressureLinearResidual = pressureSolve.relative_residual;
        }

        p += alphaP * pPrime;
        auto pressureFlux = equ::faceFlux(pressureCorrectionEquation, pPrime);
        U -= rAU * math::grad(pPrime);
        phi = phiH + pressureFlux;

        bool turbulenceConverged = true;
        double dTurbulence = 0.0, rTurbulence = 0.0;
        if (turbulence) {
            const auto result = turbulence.correct();
            dTurbulence = turbulence.relativeChange();
            rTurbulence = turbulence.relativeResidual();
            if (!diagnostics::all(result.healthy() && std::isfinite(dTurbulence) && std::isfinite(rTurbulence))) return 2;
            turbulenceConverged = result.converged()
                && dTurbulence <= turbulence.tolerance()
                && rTurbulence <= turbulence.tolerance();
        }

        const double dU = diagnostics::relativeChange(U, U_previous);
        const double dP = diagnostics::relativeMagnitude(pPrime, p);
        const auto mass = diagnostics::fluxBalance(phi);
        if (!diagnostics::all(std::isfinite(rU) && std::isfinite(dU)
            && std::isfinite(dP) && std::isfinite(mass.relative))) return 2;
        const bool linearConverged = velocitySolve.converged() && pressureConverged;
        converged = diagnostics::all(linearConverged && turbulenceConverged
            && rU <= momentumTolerance && dU <= velocityTolerance
            && dP <= pressureTolerance && mass.relative <= massTolerance);

        if (primaryProcess() && (iter == 0 || (iter + 1) % 100 == 0 || converged || iter + 1 == maxIterations))
            std::cout << "SIMPLE " << iter + 1 << " mass=" << mass.relative
                      << " dU=" << dU << " rU=" << rU << " dP=" << dP
                      << " linP=" << pressureLinearResidual
                      << " dTurb=" << dTurbulence << " rTurb=" << rTurbulence
                      << " linear=" << (linearConverged ? "ok" : "inexact")
                      << " converged=" << (converged ? "true" : "false") << '\n';
        if (converged) break;
    }
    if (!converged) return 2;
    write(problem, readTimeControl(problem).start, 0);
    return 0;
}
const SolverRegistration simple("simple", runSimple);
}
