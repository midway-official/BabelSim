#include "state.h"

namespace babelsim {
void SimpleSolver::solvePressure() {
    State& state = *m_state;
    state.requireStep(State::Step::Momentum);
    // pEqn：div(phiHbyA) = div(rAU grad(pPrime))。非正交显式修正由通用
    // eqn::laplacian 处理；这里只保留 SIMPLE 所需的语义化修正循环。
    ScalarField& p = state.m_p;
    ScalarField& pPrime = state.m_p_prime;
    ScalarField& rAU = state.m_rAU;
    ScalarField& phiHbyA = state.m_phiHbyA;
    ScalarField& divPhiHbyA = state.m_div_phiHbyA;

    state.predictMomentumFlux();
    math::evaluate(math::div(phiHbyA), divPhiHbyA);
    pPrime.fill(0.0);

    SolveResult& pressure = state.m_result.pressure;
    state.m_pressure_healthy = true;
    state.m_pressure_converged = true;
    // 正交格式求解一次；非正交格式追加指定次数的显式修正，顺序与原算法一致。
    const int pressure_solves = state.m_methods.diffusionFor(pPrime.name()) == DiffusionMethod::Orthogonal
        ? 1 : state.m_control.non_orthogonal_corrections + 1;
    for (int correction = 0; correction < pressure_solves; ++correction) {
        pressure = solve(
            -eqn::laplacian(rAU, pPrime) ==
                -eqn::source(divPhiHbyA),
            state.m_has_fixed_pressure ? EquationControl{} : referenceValue(0.0));
        state.m_pressure_healthy = state.m_pressure_healthy && pressure.healthy();
        state.m_pressure_converged = state.m_pressure_converged && pressure.converged();
    }

    p.addScaled(state.m_control.pressure_relaxation, pPrime);
    state.m_step = State::Step::Pressure;
}

void SimpleSolver::correctVelocity() {
    State& state = *m_state;
    state.requireStep(State::Step::Pressure);
    math::subtract(state.m_rAU, math::grad(state.m_p_prime), state.m_U);
    state.m_step = State::Step::Velocity;
}

void SimpleSolver::correctFlux() {
    State& state = *m_state;
    state.requireStep(State::Step::Velocity);
    state.m_phi.assign(state.m_phiHbyA);
    math::subtract(math::flux(state.m_rAU, state.m_p_prime), state.m_phi);
    state.m_step = State::Step::Flux;
}

}  // babelsim 命名空间
