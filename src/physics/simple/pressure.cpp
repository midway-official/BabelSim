#include "internal/simple_state.h"

namespace babelsim {
// 压力方程的第一次求解给出正交部分，后续求解用于显式非正交修正。调用点只表达
// “继续修正非正交项”，不再暴露 pass 编号或首末次分支。
class NonOrthogonalCorrections {
public:
    NonOrthogonalCorrections(DiffusionMethod method, int corrections)
        : m_remaining(method == DiffusionMethod::Orthogonal ? 1 : corrections + 1)
    {}

    bool correctNonOrthogonal() {
        if (m_remaining == 0) return false;
        --m_remaining;
        return true;
    }

private:
    int m_remaining;
};

void SimpleSolver::State::initializePressureCorrectionBoundaries() {
    m_has_fixed_pressure = setHomogeneousCorrectionBoundaries(
        m_algorithm.p_prime, m_p);
}

SimpleSolver::State::PressureEquationResult SimpleSolver::State::solvePressure() {
    // pEqn：div(phiHbyA) = div(rAU grad(pPrime))。非正交显式修正由通用
    // fvm::laplacian 处理；这里只保留 SIMPLE 所需的语义化修正循环。
    ScalarField& p = m_p;
    ScalarField& pPrime = m_algorithm.p_prime;
    ScalarField& rAU = m_algorithm.rAU;
    ScalarField& phiHbyA = m_algorithm.phiHbyA;
    const ScalarField& divPhiHbyA = m_workspace.div_phiHbyA;

    predictMomentumFlux();
    m_workspace.preparePressureEquation(phiHbyA);
    pPrime.fill(0.0);

    PressureEquationResult result;
    result.healthy = true;
    result.linear_converged = true;
    NonOrthogonalCorrections corrections(
        m_methods.diffusionFor(pPrime.name()), m_control.non_orthogonal_corrections);
    while (corrections.correctNonOrthogonal()) {
        result.linear = simple::solvePressureCorrectionEquation(
            -fvm::laplacian(rAU, pPrime) ==
                -fvm::source(divPhiHbyA),
            !m_has_fixed_pressure);
        result.healthy = result.healthy && result.linear.healthy();
        result.linear_converged = result.linear_converged && result.linear.converged();
    }

    p.addScaled(m_control.pressure_relaxation, pPrime);
    return result;
}

void SimpleSolver::State::correctVelocity() {
    VectorField& U = m_U;
    const ScalarField& pPrime = m_algorithm.p_prime;
    const ScalarField& rAU = m_algorithm.rAU;
    fvc::subtract(rAU, fvc::grad(pPrime), U);
}

void SimpleSolver::State::correctFlux() {
    ScalarField& phi = m_phi;
    const ScalarField& phiHbyA = m_algorithm.phiHbyA;
    const ScalarField& pPrime = m_algorithm.p_prime;
    const ScalarField& rAU = m_algorithm.rAU;

    phi.assign(phiHbyA);
    fvc::subtract(fvc::flux(rAU, pPrime), phi);
}

}  // babelsim 命名空间
