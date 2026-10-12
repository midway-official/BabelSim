#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/builtin_descriptions.h"
#include "babelsim/time_discretization.h"
#include "babelsim/geometry.h"
#include "babelsim/solver.h"
#include "../RANS/api.h"
#include "babelsim/monitor.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace babelsim {
namespace {

// Static mesh, volumetric phi. These are the same old-time RHS weights as
// equ::ddt, including Euler startup and a shortened/variable BDF2 step.
// See OpenFOAM-10 backwardDdtScheme::fvcDdtPhiCorr and ddtScheme::fvcDdtPhiCoeff.
ScalarField temporalFluxCorrection(
    const time::History<Vec3>& velocity,
    const time::History<double>& flux,
    TimeMethod method, const OperatorOptions& options)
{
    time::requireMatching(velocity, flux);
    if (method == TimeMethod::Steady)
        throw std::invalid_argument("PISO requires transient histories");
    const auto weights = time::weights(method, velocity);
    const double previousWeight = weights.previous, olderWeight = weights.older;
    const auto mismatch = flux.previous() - math::flux(velocity.previous(), options);
    auto correction = previousWeight * mismatch;
    if (olderWeight != 0.0)
        correction += olderWeight * (flux.older() - math::flux(velocity.older(), options));

    // Flux-normalised coupling: avoid amplifying a history discrepancy on a
    // stagnant/reversing face. No dimensional epsilon changes the flow scale.
    const auto coupling = fieldBinary(mismatch, flux.previous(), [](double e, double f) {
        if (e == 0.0) return 1.0;
        if (f == 0.0) return 0.0;
        return std::max(0.0, 1.0 - std::abs(e) / std::abs(f));
    });
    ScalarField result(flux.previous().mesh(), FieldLocation::Face, "piso.ddtFlux");
    result.fill(0.0);
    // Includes inter-partition faces; prescribed physical boundary fluxes are
    // handled by U/p constraints, never by old-time interpolation corrections.
    math::add(coupling * correction, result, math::FaceRegion::Interior);
    return result;
}

} // namespace

// 瞬态 PISO。每遍预测--修正外迭代组装并求解一次动量预测方程，随后用 nCorrectors
// 次校正重新计算非对角动量作用 H(U)，再投影通量；不能重复投影已经守恒的 phi。
// 修正步不做压力欠松弛。瞬态 SIMPLE 则相反：
// 步内反复做带欠松弛的动量/压力迭代直到收敛，没有独立修正结构。
// maxIterations 把整个预测--修正流程重复为 Picard 外层迭代，只用于非线性项与
// 湍流模型的耦合，默认 1 即标准 PISO。
// 动量预测可以欠松弛：修正步始终施加完整修正，因此修正后通量的守恒性与
// velocityRelaxation 无关。默认 1 为不欠松弛的 PISO；若指定小于 1，使用实际
// 松弛后矩阵的对角响应，不再通过 SIMPLE 风格的额外缩放取消该响应。
SolverResult runPiso(Case& problem) {
    const monitor::Reporter reporter("PISO");
    const auto& spec = builtin::piso();
    const auto config = problem.read.config(spec);
    const auto& methods = problem.read.methods();

    auto& U = problem.load.vectorField("U");
    auto& p = problem.load.scalarField("p");
    auto& phi = problem.create.faceScalarField("phi");
    auto& pPrime = problem.create.homogeneousCorrection(p, "pPrime");
    auto momentumEquation = equ::createEquation(problem, spec.equation("momentum"), U);
    const auto pressureControl = readEquationControl(problem, spec.equation("pressureCorrection"), p);
    auto pressureEquation = equ::createEquation(p, pressureControl);
    // The auxiliary initial projection reuses the pressure numerical profile.
    auto potentialEquation = equ::createEquation(pPrime, pressureControl);
    const auto pressureOptions = pressureEquation.options("diffusion");
    const int pressureSolves = pressureOptions.diffusion == DiffusionMethod::Orthogonal
        ? 1 : config.coupling.nonOrthogonalCorrections + 1;
    const int initialPressureSolves = pressureOptions.diffusion == DiffusionMethod::Orthogonal
        ? 1 : config.coupling.initialNonOrthogonalCorrections + 1;

    phi = math::flux(U, pressureOptions);
    if (config.coupling.initializePotentialFlow) {
        // Establish a conservative initial face flux before it enters convection
        // and the first physical time history. The potential is not pressure.
        const auto initialImbalance = math::div(phi);
        for (int correction = 0; correction < initialPressureSolves; ++correction) {
            potentialEquation.reset();
            equ::laplacian(potentialEquation, 1.0, -1, "diffusion");
            equ::source(potentialEquation, -initialImbalance);
            potentialEquation.referenceIfUnanchored(0.0);
            const auto initialSolve = equ::solve(potentialEquation);
            if (!diagnostics::all(initialSolve.converged())) return SolverResult::numericalFailure();
        }
        U -= math::grad(pPrime, pressureOptions);
        phi += equ::faceFlux(potentialEquation, pPrime);
        U.setBoundaryFlux(phi);
        pPrime.fill(0.0);
        reporter.record({{"initialMass", diagnostics::fluxBalance(phi).relative},
            {"initialUmax", math::max(math::map(U, [](Vec3 v) { return norm(v); }))},
            {"initialPressureSolves", initialPressureSolves}});
    }
    // Model inputs and histories use the final projected initial state.
    auto turbulence = rans::create(problem, U, phi, config.fluid, momentumEquation.options());
    const auto& muEff = turbulence.effectiveViscosity();
    const bool segregatedTurbulence = config.turbulence.mode == "segregated";
    if (segregatedTurbulence && config.coupling.maxIterations != 1)
        throw std::invalid_argument("segregated turbulence requires one PISO momentum pass");
    if (segregatedTurbulence && turbulence
        && turbulence.relaxation() != 1.0)
        throw std::invalid_argument("single-pass transient turbulence requires turbulenceRelaxation 1");
    // With one momentum pass, converge the model's relaxed transport equations
    // at the SAME time level. Otherwise relaxation changes physical time rates.
    // Group model Picard solves before rebuilding the flow equations. History
    // is still saved exactly once per physical step. The legacy multi-pass
    // default remains one model correction per flow pass.
    const int turbulenceCorrections = config.turbulence.correctors ? *config.turbulence.correctors
        : (turbulence && config.coupling.maxIterations == 1 && !segregatedTurbulence ? config.turbulence.maxIterations : 1);
    if (segregatedTurbulence && turbulenceCorrections != 1)
        throw std::invalid_argument("segregated turbulence requires one transport correction");
    const auto V = geometry::cellVolumes(problem.mesh());
    const int writeInterval = readWriteInterval(problem);
    time::HistoryGroup histories;
    auto& U_old = histories.addHistory(U);
    auto& phi_old = histories.addHistory(phi);
    histories.trackForRollback(p);
    turbulence.registerTimeState(histories);
    auto time = time::makeStepper(problem);
    problem.finalizeInitialization();

    while (time.value() < time.end()) {
        time.advance();
        histories.save(time);
        turbulence.beginStep();
        const auto temporalFlux = config.fluid.density * temporalFluxCorrection(
            U_old, phi_old, methods.time, pressureOptions);

        bool converged = false;
        for (int iter = 0; iter < config.coupling.maxIterations; ++iter) {
            const VectorField U_previous = U; // 上一个外层迭代值，不是物理时间历史

            // 动量预测：每个时间步求解一次，除非要求额外的外层迭代。
            momentumEquation.reset();
            equ::ddt(momentumEquation, config.fluid.density, U_old);
            equ::div(momentumEquation, phi, config.fluid.density, "convection");
            equ::laplacian(momentumEquation, muEff, -1, "diffusion");
            const auto predictorGradP = math::grad(p, pressureOptions);
            equ::source(momentumEquation, -predictorGradP);
            if (turbulence)
            equ::source(momentumEquation, math::div(
                turbulence.deviatoricStressRemainder(U), momentumEquation.options()));

            // 在松弛或预测之前测量当前外层迭代的方程残差。
            const double rU = diagnostics::relativeResidual(momentumEquation, U);
            if (!diagnostics::all(std::isfinite(rU))) return SolverResult::numericalFailure();
            equ::relax(momentumEquation, U_previous, config.coupling.velocityRelaxation);
            const auto aP = momentumEquation.diagonal();
            const auto rAU = V / aP;
            const auto inverseAP = 1.0 / aP;
            const auto rAUf = math::interpolate(rAU, pressureOptions.coefficientOptions());
            // Integrated non-pressure RHS of the frozen momentum equation.
            // Relaxation changes aP and b; the pressure source remains -V grad(p).
            const auto nonPressureRhs = momentumEquation.rhs() + V * predictorGradP;
            const auto velocitySolve = equ::solve(momentumEquation);
            if (!diagnostics::all(velocitySolve.healthy())) return SolverResult::numericalFailure();

            bool pressureConverged = true;
            double pressureLinearResidual = 0.0;
            double dP = 0.0, firstCouplingResidual = 0.0, couplingResidual = 0.0;
            int correctionsUsed = 0;
            for (int corrector = 0; corrector < config.coupling.maxCorrectors; ++corrector) {
                // A = diag(aP) + N. HbyA = (b0 - N U)/aP must use the
                // current corrected U on EVERY pass. equ::apply synchronises
                // the full matrix action without exposing backend storage.
                VectorField HbyA = U; // retain physical U boundary constraints
                HbyA = inverseAP * (nonPressureRhs - equ::apply(momentumEquation, U) + aP * U);
                HbyA.setBoundaryFlux(phi);
                const ScalarField previousPressure = p;
                auto phiH = math::flux(HbyA, pressureOptions);
                math::add(rAUf * temporalFlux, phiH, math::FaceRegion::Interior);
                const auto imbalance = math::div(phiH);
                // Solve for total pressure. A limited non-orthogonal flux is
                // nonlinear: L(p + dp) != L(p) + L(dp). Splitting those fluxes
                // leaves the corrected velocity and conservative phi inconsistent.
                for (int correction = 0; correction < pressureSolves; ++correction) {
                    pressureEquation.reset();
                    equ::laplacian(pressureEquation, rAUf, -1, "diffusion");
                    equ::source(pressureEquation, -imbalance);
                    pressureEquation.referenceIfUnanchored(0.0);
                    const auto pressureSolve = equ::solve(pressureEquation);
                    if (!diagnostics::all(pressureSolve.healthy())) return SolverResult::numericalFailure();
                    pressureConverged = pressureConverged && pressureSolve.converged();
                    pressureLinearResidual = std::max(pressureLinearResidual, pressureSolve.relative_residual);
                }
                // Use the same frozen explicit correction as the assembled
                // pressure equation, never a separately reconstructed flux.
                phi = phiH + equ::faceFlux(pressureEquation, p);
                const auto correctedGradP = math::grad(p, pressureOptions);
                U = HbyA;
                math::subtract(rAU, correctedGradP, U);
                U.setBoundaryFlux(phi);
                dP = std::max(dP, diagnostics::relativeMagnitude(p - previousPressure, p));
                const auto action = equ::apply(momentumEquation, U);
                const auto pressureSource = V * correctedGradP;
                const auto defect = nonPressureRhs - action - pressureSource;
                const auto norms = math::normsL2({nonPressureRhs, action, pressureSource, defect});
                const double scale = norms[0] + norms[1] + norms[2];
                couplingResidual = scale > 0.0 ? norms[3] / scale : norms[3];
                if (corrector == 0) firstCouplingResidual = couplingResidual;
                correctionsUsed = corrector + 1;
                const bool correctionDone = correctionsUsed >= config.coupling.correctors
                    && diagnostics::all(couplingResidual <= config.coupling.tolerance);
                if (corrector == 0 || correctionDone || correctionsUsed == config.coupling.maxCorrectors)
                    reporter.record({{"corrector", corrector + 1},
                        {"Umax", math::max(math::map(U, [](Vec3 v) { return norm(v); }))},
                        {"coupling", couplingResidual}});
                if (correctionDone) break;
            }

            bool turbulenceConverged = true;
            bool turbulenceLinearConverged = true;
            double dTurbulence = 0.0, rTurbulence = 0.0;
            if (turbulence) {
                for (int correction = 0; correction < turbulenceCorrections; ++correction) {
                    const auto result = turbulence.solveTransport();
                    dTurbulence = result.relativeChange();
                    rTurbulence = result.initialResidual();
                    if (!diagnostics::all(result.healthy() && std::isfinite(dTurbulence) && std::isfinite(rTurbulence))) return SolverResult::numericalFailure();
                    turbulenceLinearConverged = diagnostics::all(result.linearConverged());
                    turbulenceConverged = diagnostics::all(turbulenceLinearConverged
                        && dTurbulence <= turbulence.tolerance()
                        && rTurbulence <= turbulence.tolerance());
                    if (!turbulenceLinearConverged || turbulenceConverged) break;
                }
            }

            const double dU = diagnostics::relativeChange(U, U_previous);
            const auto mass = diagnostics::fluxBalance(phi);
            const auto courant = diagnostics::courantNumber(phi, time.dt());
            if (!diagnostics::all(std::isfinite(rU) && std::isfinite(dU)
                && std::isfinite(dP) && std::isfinite(couplingResidual)
                && std::isfinite(mass.relative))) return SolverResult::numericalFailure();
            const bool linearConverged = velocitySolve.converged() && pressureConverged
                && turbulenceLinearConverged;
            // A conservative flux does not excuse a failed linear solve.
            // Physical changes between time levels are not steady convergence.
            const bool stepAccepted = diagnostics::all(linearConverged && mass.relative <= config.convergence.continuity
                && couplingResidual <= config.coupling.tolerance);
            // 单次 PISO 不要求物理速度变化趋零；湍流的松弛输运仍须在本时间层收敛。
            const bool outerSettled = diagnostics::all(rU <= config.convergence.momentum
                && dU <= config.convergence.velocity && dP <= config.convergence.pressureCorrection && turbulenceConverged);
            converged = stepAccepted && (config.coupling.maxIterations == 1
                ? (segregatedTurbulence || turbulenceConverged) : outerSettled);

            reporter.iteration(iter + 1, config.coupling.maxIterations, converged, {
                {"time", time.value()}, {"correctors", correctionsUsed},
                {"CoMax", courant.maximum}, {"CoMean", courant.mean},
                {"accepted", stepAccepted}, {"mass", mass.relative}, {"dU", dU}, {"rU", rU}, {"dP", dP},
                {"rCouplingFirst", firstCouplingResidual}, {"rCoupling", couplingResidual},
                {"linU", velocitySolve.relative_residual}, {"linP", pressureLinearResidual},
                {"dTurb", dTurbulence}, {"rTurb", rTurbulence},
                {"turbulenceSettled", turbulenceConverged},
                {"linear", linearConverged ? "ok" : "inexact"}, {"converged", converged}});
            if (!stepAccepted) return SolverResult::notConverged();
            if (converged) break;
        }
        if (!converged) return SolverResult::notConverged();
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return SolverResult::completed();
}
const SolverRegistration piso("piso", runPiso, &builtin::piso());
}
