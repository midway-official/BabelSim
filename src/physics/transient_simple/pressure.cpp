#include "state.h"

namespace babelsim {

void TransientSimpleAlgorithm::solvePressure() {
    State& state = *m_state;
    state.requireStep(State::Step::Momentum);

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
    const int pressure_solves =
        state.m_methods.diffusionFor(pPrime.name()) == DiffusionMethod::Orthogonal
            ? 1 : state.m_control.non_orthogonal_corrections + 1;
    for (int correction = 0; correction < pressure_solves; ++correction) {
        pressure = solve(
            -eqn::laplacian(rAU, pPrime) == -eqn::source(divPhiHbyA),
            state.m_has_fixed_pressure ? EquationControl{} : referenceValue(0.0));
        state.m_pressure_healthy = state.m_pressure_healthy && pressure.healthy();
        state.m_pressure_converged = state.m_pressure_converged && pressure.converged();
    }

    p.addScaled(state.m_control.pressure_relaxation, pPrime);
    state.m_step = State::Step::Pressure;
}

void TransientSimpleAlgorithm::correctVelocity() {
    State& state = *m_state;
    state.requireStep(State::Step::Pressure);
    math::subtract(state.m_rAU, math::grad(state.m_p_prime), state.m_U);
    state.m_step = State::Step::Velocity;
}

void TransientSimpleAlgorithm::correctFlux() {
    State& state = *m_state;
    state.requireStep(State::Step::Velocity);
    state.m_phi.assign(state.m_phiHbyA);
    math::subtract(math::flux(state.m_rAU, state.m_p_prime), state.m_phi);
    state.m_step = State::Step::Flux;
}

void TransientSimpleAlgorithm::correctTurbulence() {
    State& state = *m_state;
    state.requireStep(State::Step::Flux);
    SimpleIterationResult& result = state.m_result;
    if (state.m_turbulence) {
        result.turbulence = rans::correct(*state.m_turbulence);
        result.relative_turbulence_change = rans::relativeChange(*state.m_turbulence);
    } else {
        result.turbulence = {SolveStatus::Converged, 0, 0.0, 0.0, 0.0};
        result.relative_turbulence_change = 0.0;
    }
    state.m_step = State::Step::Turbulence;
}

}  // babelsim 命名空间
