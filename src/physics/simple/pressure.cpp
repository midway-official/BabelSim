#include "state.h"

namespace babelsim {
SimpleSolver::State::PressureEquationResult SimpleSolver::State::solvePressure() {
    // pEqn：div(phiHbyA) = div(rAU grad(pPrime))。非正交显式修正由通用
    // eqn::laplacian 处理；这里只保留 SIMPLE 所需的语义化修正循环。
    ScalarField& p = m_p;
    ScalarField& pPrime = m_algorithm.p_prime;
    ScalarField& rAU = m_algorithm.rAU;
    ScalarField& phiHbyA = m_algorithm.phiHbyA;
    const ScalarField& divPhiHbyA = m_workspace.div_phiHbyA;

    predictMomentumFlux();
    math::evaluate(math::div(phiHbyA), m_workspace.div_phiHbyA);
    pPrime.fill(0.0);

    PressureEquationResult result;
    result.healthy = true;
    result.linear_converged = true;
    // 正交格式求解一次；非正交格式追加指定次数的显式修正，顺序与原算法一致。
    const int pressure_solves = m_methods.diffusionFor(pPrime.name()) == DiffusionMethod::Orthogonal
        ? 1 : m_control.non_orthogonal_corrections + 1;
    for (int correction = 0; correction < pressure_solves; ++correction) {
        result.linear = solve(
            -eqn::laplacian(rAU, pPrime) ==
                -eqn::source(divPhiHbyA),
            m_has_fixed_pressure ? EquationControl{} : referenceValue(0.0));
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
    math::subtract(rAU, math::grad(pPrime), U);
}

void SimpleSolver::State::correctFlux() {
    ScalarField& phi = m_phi;
    const ScalarField& phiHbyA = m_algorithm.phiHbyA;
    const ScalarField& pPrime = m_algorithm.p_prime;
    const ScalarField& rAU = m_algorithm.rAU;

    phi.assign(phiHbyA);
    math::subtract(math::flux(rAU, pPrime), phi);
}

}  // babelsim 命名空间
