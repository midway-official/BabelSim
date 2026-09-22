#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/geometry.h"
#include "babelsim/solver.h"
#include "../RANS/api.h"
#include "babelsim/monitor.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace babelsim {

// 瞬态 PISO。每个物理时间步只组装并求解一次动量预测方程，随后用 nCorrectors
// 次压力修正依次消掉上一次修正后通量上仍然存在的质量不平衡；修正步不做压力
// 欠松弛，这是 PISO 的算子分裂（predictor + correctors）。瞬态 SIMPLE 则相反：
// 步内反复做带欠松弛的动量/压力迭代直到收敛，没有独立修正结构。
// maxIterations 把整个预测--修正流程重复为 Picard 外层迭代，只用于非线性项与
// 湍流模型的耦合，默认 1 即标准 PISO。
// 动量预测可以欠松弛：修正步始终施加完整修正，因此修正后通量的守恒性与
// velocityRelaxation 无关，该系数只决定耦合路径与时间步上限。默认取 0.7 与
// SIMPLE 家族一致；取 1.0 是不做松弛的教科书 PISO，Co 偏大时容易失稳。
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
    // 外层每遍的压力修正次数；源项是上一次修正后剩余的质量不平衡。
    const int correctors = settings.integer("nCorrectors", 2, 1, 100);
    const int nonOrthogonalCorrections = settings.integer("nonOrthogonalCorrections", 1, 0, 20);
    const double alphaU = settings.fraction("velocityRelaxation", 0.7);
    const double massTolerance = settings.positive("continuityTolerance", 1e-8);
    const double velocityTolerance = settings.positive("velocityTolerance", 1e-7);
    const double momentumTolerance = settings.positive("momentumTolerance", 1e-6);
    const double pressureTolerance = settings.positive("pressureCorrectionTolerance", 1e-6);

    auto turbulence = rans::load(problem, U, phi);
    const auto& muEff = turbulence.effectiveViscosity();
    auto pPrime = field::homogeneousLike(p, "pPrime");
    auto momentumEquation = equ::createEquation(problem, "momentum", U, {"convection", "diffusion"});
    auto pressureEquation = equ::createEquation(problem, "pressureCorrection", pPrime, {"diffusion"});
    const auto pressureOptions = pressureEquation.options("diffusion");
    const int pressureSolves = pressureOptions.diffusion == DiffusionMethod::Orthogonal
        ? 1 : nonOrthogonalCorrections + 1;
    phi = math::flux(U, momentumEquation.options());
    if (methods.time == TimeMethod::Steady)
        throw std::invalid_argument("PISO requires a transient time scheme");
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
            const VectorField U_previous = U; // 上一个外层迭代值，不是物理时间历史

            // 动量预测：每个时间步求解一次，除非要求额外的外层迭代。
            momentumEquation.reset();
            equ::ddt(momentumEquation, rho, U_old);
            equ::div(momentumEquation, phi, rho, "convection");
            equ::laplacian(momentumEquation, muEff, -1, "diffusion");
            equ::source(momentumEquation, -math::grad(p, pressureOptions));
            if (turbulence)
            equ::source(momentumEquation, math::div(
                turbulence.deviatoricStressRemainder(U), momentumEquation.options()));

            // 在松弛或预测之前测量当前外层迭代的方程残差。
            const double rU = diagnostics::relativeResidual(momentumEquation, U);
            if (!diagnostics::all(std::isfinite(rU))) return SolverResult::numericalFailure();
            equ::relax(momentumEquation, U_previous, alphaU);
            // 保持 SIMPLE 行归一化；下面的 V/aP 使用这些缩放后的行。
            equ::scale(momentumEquation, alphaU);
            const auto aP = momentumEquation.diagonal();
            const auto rAU = V / aP;
            const auto velocitySolve = equ::solve(momentumEquation);
            if (!diagnostics::all(velocitySolve.healthy())) return SolverResult::numericalFailure();

            // Rhie–Chow 插值；物理边界通量保持不变。
            const auto gradP = math::grad(p, pressureOptions);
            auto phiH = math::flux(U, momentumEquation.options());
            const auto interpolatedGradientFlux =
                math::dot(math::interpolate(rAU * gradP, pressureOptions), Sf);
            const auto normalGradientFlux =
                math::interpolate(rAU, pressureOptions.coefficientOptions()) * math::normalGradient(p, gradP, pressureOptions) * Af;
            math::add(interpolatedGradientFlux, phiH, math::FaceRegion::Interior);
            math::subtract(normalGradientFlux, phiH, math::FaceRegion::Interior);
            phi = phiH;

            // PISO 修正结构：第 k 次修正的压力方程以当前通量上剩余的质量不平衡为
            // 源项，并把完整修正量叠加到 p、U 和 phi 上。累加修正通量与用修正后的
            // 速度、压力重新构造 Rhie–Chow 一致通量在数学上等价，因此 RC 项每个
            // 外层只用构造一次；压力修正本身不做欠松弛。
            bool pressureConverged = true;
            double pressureLinearResidual = 0.0;
            for (int corrector = 0; corrector < correctors; ++corrector) {
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
                    pressureLinearResidual = pressureSolve.relative_residual;
                }
                const auto correctionFlux = equ::faceFlux(pressureEquation, pPrime);
                p += pPrime;
                U -= rAU * math::grad(pPrime, pressureOptions);
                phi += correctionFlux;
            }

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
            // 时间步的接受判据只有守恒：PISO 的时间步由预测--修正分裂定义，动量
            // 预测和湍流输运每个外层都只推进一步，它们的线性求解容差（linU/linP）
            // 与湍流稳态容差（dTurb/rTurb）只报告，不影响修正后通量的守恒性。
            const bool stepAccepted = diagnostics::all(mass.relative <= massTolerance);
            // 单次预测--修正已经是完整的 PISO 时间步；只有额外的外层迭代才要求
            // 动量、速度和压力修正的外层变化以及湍流输运在步内收敛。
            const bool outerSettled = diagnostics::all(rU <= momentumTolerance
                && dU <= velocityTolerance && dP <= pressureTolerance && turbulenceConverged);
            converged = stepAccepted && (maxIterations == 1 || outerSettled);

            reporter.iteration(iter + 1, maxIterations, converged, {
                {"mass", mass.relative}, {"dU", dU}, {"rU", rU}, {"dP", dP},
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
