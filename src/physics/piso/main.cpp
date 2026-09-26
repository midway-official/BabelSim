#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/geometry.h"
#include "babelsim/solver.h"
#include "../RANS/api.h"
#include "babelsim/monitor.h"
#include <algorithm>
#include <cmath>
#include <limits>
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
    if (velocity.levels() == 0 || velocity.levels() != flux.levels()
        || velocity.dt() != flux.dt() || velocity.previousDt() != flux.previousDt()
        || method == TimeMethod::Steady)
        throw std::invalid_argument("PISO requires matching transient U/phi histories");
    const double dt = velocity.dt();
    double previousWeight = 1.0 / dt, olderWeight = 0.0;
    if (method == TimeMethod::BDF2 && velocity.levels() >= 2) {
        const double ratio = dt / velocity.previousDt();
        previousWeight = (1.0 + ratio) / dt;
        olderWeight = -ratio * ratio / ((1.0 + ratio) * dt);
    }
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

// 瞬态 PISO。每个物理时间步只组装并求解一次动量预测方程，随后用 nCorrectors
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
    auto& U = problem.vectorField("U");
    auto& p = problem.scalarField("p");
    auto& phi = problem.createFaceScalarField("phi");
    const auto& physical = problem.physics();
    const double rho = physical.positive("density");

    const auto& methods = problem.methods(); // 方法字典由 Case/runtime 解析一次
    const auto V = geometry::cellVolumes(problem.mesh());
    const auto Sf = geometry::faceAreaVectors(problem.mesh());
    const auto Af = geometry::faceAreas(problem.mesh());
    const auto& settings = problem.solution();
    // 每个时间步的预测--修正外层遍数；1 就是标准 PISO。
    const int maxIterations = settings.integer("maxIterations", 1, 1, std::numeric_limits<int>::max());
    // 外层每遍的动量/压力校正次数，每次用更新后的 U 重建 HbyA。
    const int correctors = settings.integer("nCorrectors", 2, 1, 100);
    const int nonOrthogonalCorrections = settings.integer("nonOrthogonalCorrections", 1, 0, 20);
    const double alphaU = settings.fraction("velocityRelaxation", 1.0);
    const double massTolerance = settings.positive("continuityTolerance", 1e-8);
    const double velocityTolerance = settings.positive("velocityTolerance", 1e-7);
    const double momentumTolerance = settings.positive("momentumTolerance", 1e-6);
    const double pressureTolerance = settings.positive("pressureCorrectionTolerance", 1e-6);

    auto turbulence = rans::load(problem, U, phi);
    // With one momentum pass, converge the model's relaxed transport equations
    // at the SAME time level. Otherwise relaxation changes physical time rates.
    const int turbulenceCorrections = turbulence && maxIterations == 1
        ? settings.integer("turbulenceMaxIterations", 1000, 1, std::numeric_limits<int>::max()) : 1;
    const auto& muEff = turbulence.effectiveViscosity();
    auto pPrime = field::homogeneousLike(p, "pPrime");
    auto momentumEquation = equ::createEquation(problem, "momentum", U, {"convection", "diffusion"});
    auto pressureEquation = equ::createEquation(problem, "pressureCorrection", pPrime, {"diffusion"});
    const auto pressureOptions = pressureEquation.options("diffusion");
    const int pressureSolves = pressureOptions.diffusion == DiffusionMethod::Orthogonal
        ? 1 : nonOrthogonalCorrections + 1;
    phi = math::flux(U, pressureOptions);
    if (methods.time == TimeMethod::Steady)
        throw std::invalid_argument("PISO requires a transient time scheme");
    auto time = time::start(problem);
    const int writeInterval = readWriteInterval(problem);
    auto U_old = time::history(U);
    auto phi_old = time::history(phi);
    problem.validate();

    while (time.value() < time.end()) {
        time.advance();
        U_old.save(U, time.dt());
        phi_old.save(phi, time.dt());
        if (turbulence) turbulence.saveOld(time.dt());
        const auto temporalFlux = rho * temporalFluxCorrection(
            U_old, phi_old, methods.time, pressureOptions);

        bool converged = false;
        for (int iter = 0; iter < maxIterations; ++iter) {
            const VectorField U_previous = U; // 上一个外层迭代值，不是物理时间历史

            // 动量预测：每个时间步求解一次，除非要求额外的外层迭代。
            momentumEquation.reset();
            equ::ddt(momentumEquation, rho, U_old);
            equ::div(momentumEquation, phi, rho, "convection");
            equ::laplacian(momentumEquation, muEff, -1, "diffusion");
            const auto predictorGradP = math::grad(p, pressureOptions);
            equ::source(momentumEquation, -predictorGradP);
            if (turbulence)
            equ::source(momentumEquation, math::div(
                turbulence.deviatoricStressRemainder(U), momentumEquation.options()));

            // 在松弛或预测之前测量当前外层迭代的方程残差。
            const double rU = diagnostics::relativeResidual(momentumEquation, U);
            if (!diagnostics::all(std::isfinite(rU))) return SolverResult::numericalFailure();
            equ::relax(momentumEquation, U_previous, alphaU);
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
            for (int corrector = 0; corrector < correctors; ++corrector) {
                // A = diag(aP) + N. HbyA = (b0 - N U)/aP must use the
                // current corrected U on EVERY pass. equ::apply synchronises
                // the full matrix action without exposing backend storage.
                VectorField HbyA = U; // retain physical U boundary constraints
                HbyA = inverseAP * (nonPressureRhs - equ::apply(momentumEquation, U) + aP * U);
                HbyA.setBoundaryFlux(phi);
                const auto gradP = math::grad(p, pressureOptions);
                // Use calculated traces for this known predictor so its
                // interpolated pressure term cancels also on physical faces.
                // HbyA itself retains the physical velocity constraints.
                const auto predictor = HbyA - rAU * gradP;
                auto phiH = math::flux(predictor, pressureOptions);
                const auto interpolatedGradientFlux =
                    math::dot(math::interpolate(rAU * gradP, pressureOptions), Sf);
                math::add(interpolatedGradientFlux, phiH);
                math::add(rAUf * temporalFlux, phiH, math::FaceRegion::Interior);
                // Include pressure Dirichlet boundary fluxes (e.g. outlet).
                // Zero-gradient/symmetry pressure contributes zero at walls.
                phi = phiH - rAUf * math::normalGradient(p, gradP, pressureOptions) * Af;
                U = predictor;
                pPrime.fill(0.0);
                const auto imbalance = math::div(phi);
                for (int correction = 0; correction < pressureSolves; ++correction) {
                    pressureEquation.reset();
                    equ::laplacian(pressureEquation, rAU, -1, "diffusion");
                    equ::source(pressureEquation, -imbalance);
                    pressureEquation.referenceIfUnanchored(0.0);
                    const auto pressureSolve = equ::solve(pressureEquation);
                    if (!diagnostics::all(pressureSolve.healthy())) return SolverResult::numericalFailure();
                    pressureConverged = pressureConverged && pressureSolve.converged();
                    pressureLinearResidual = std::max(pressureLinearResidual, pressureSolve.relative_residual);
                }
                const auto correctionFlux = equ::faceFlux(pressureEquation, pPrime);
                p += pPrime;
                U -= rAU * math::grad(pPrime, pressureOptions);
                phi += correctionFlux;
                U.setBoundaryFlux(phi);
                dP = std::max(dP, diagnostics::relativeMagnitude(pPrime, p));
                const auto action = equ::apply(momentumEquation, U);
                const auto pressureSource = V * math::grad(p, pressureOptions);
                const double scale = math::normL2(nonPressureRhs)
                    + math::normL2(action) + math::normL2(pressureSource);
                const double defect = math::normL2(nonPressureRhs - action - pressureSource);
                couplingResidual = scale > 0.0 ? defect / scale : defect;
                if (corrector == 0) firstCouplingResidual = couplingResidual;
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
            if (!diagnostics::all(std::isfinite(rU) && std::isfinite(dU)
                && std::isfinite(dP) && std::isfinite(couplingResidual)
                && std::isfinite(mass.relative))) return SolverResult::numericalFailure();
            const bool linearConverged = velocitySolve.converged() && pressureConverged
                && turbulenceLinearConverged;
            // A conservative flux does not excuse a failed linear solve.
            // Physical changes between time levels are not steady convergence.
            const bool stepAccepted = diagnostics::all(linearConverged && mass.relative <= massTolerance);
            // 单次 PISO 不要求物理速度变化趋零；湍流的松弛输运仍须在本时间层收敛。
            const bool outerSettled = diagnostics::all(rU <= momentumTolerance
                && dU <= velocityTolerance && dP <= pressureTolerance && turbulenceConverged);
            converged = stepAccepted && (maxIterations == 1 ? turbulenceConverged : outerSettled);

            reporter.iteration(iter + 1, maxIterations, converged, {
                {"mass", mass.relative}, {"dU", dU}, {"rU", rU}, {"dP", dP},
                {"rCouplingFirst", firstCouplingResidual}, {"rCoupling", couplingResidual},
                {"linU", velocitySolve.relative_residual}, {"linP", pressureLinearResidual},
                {"dTurb", dTurbulence}, {"rTurb", rTurbulence},
                {"linear", linearConverged ? "ok" : "inexact"}, {"converged", converged}});
            if (!stepAccepted) return SolverResult::notConverged();
            if (converged) break;
        }
        if (!converged) return SolverResult::notConverged();
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return SolverResult::completed();
}
const SolverRegistration piso("piso", runPiso);
}
