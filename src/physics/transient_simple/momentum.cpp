#include "state.h"

#include <stdexcept>

namespace babelsim {

void TransientSimpleAlgorithm::solveMomentum() {
    State& state = *m_state;
    if (state.m_iteration < 0)
        throw std::logic_error("call beginTimeStep before solving transient momentum");
    if (state.m_step == State::Step::Complete) state.m_step = State::Step::Ready;
    state.requireStep(State::Step::Ready);
    state.m_previous_velocity.assign(state.m_U);

    state.m_result.velocity = solveWithResponse(state.momentumEquation(), state.m_rAU,
        relaxed(state.m_control.velocity_relaxation));
    state.m_step = State::Step::Momentum;
}

VectorEquationDefinition TransientSimpleAlgorithm::State::momentumEquation() {
    VectorExpression diffusion = m_turbulence
        ? eqn::laplacian(m_effective_viscosity, m_U)
        : eqn::laplacian(m_fluid.dynamic_viscosity, m_U);
    if (m_turbulence) {
        evaluateStressCorrection(m_U, m_phi, m_effective_viscosity,
            m_stress_gradient, m_stress_correction, m_stress_divergence);
        diffusion = diffusion + eqn::source(m_stress_divergence);
    }
    return eqn::ddt(m_fluid.density, m_U) +
        eqn::div(m_fluid.density, m_phi, m_U) == -math::grad(m_p) + diffusion;
}

void TransientSimpleAlgorithm::State::predictMomentumFlux() {
    // 沿用稳态 SIMPLE 的同位网格动量插值，但工作场由瞬态算法独立拥有。
    math::evaluate(math::grad(m_p), m_grad_p);
    math::evaluate(math::flux(m_U), m_phiHbyA);
    m_rAU_grad_p.assignProduct(m_rAU, m_grad_p);
    math::evaluate(math::interpolate(m_rAU_grad_p), m_rAU_grad_p_face);
    math::evaluate(math::interpolate(m_rAU), m_rAU_face);
    math::add(math::flux(m_rAU_grad_p_face), m_phiHbyA, math::FaceRegion::Interior);
    math::subtract(math::flux(m_rAU_face, math::reconstruct(m_p, m_grad_p)),
                   m_phiHbyA, math::FaceRegion::Interior);
}

}  // babelsim 命名空间
