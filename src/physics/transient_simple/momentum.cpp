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

    VectorField& U = state.m_U;
    ScalarField& p = state.m_p;
    ScalarField& phi = state.m_phi;
    ScalarField& rAU = state.m_rAU;
    const VectorExpression diffusion = state.m_turbulence
        ? eqn::laplacian(state.m_effective_viscosity, U)
        : eqn::laplacian(state.m_fluid.dynamic_viscosity, U);

    // 与稳态动量方程使用相同的顶层表达，瞬态版本只增加物理时间导数。
    state.m_result.velocity = solveWithResponse(
        eqn::ddt(state.m_fluid.density, U) +
            eqn::div(state.m_fluid.density, phi, U) ==
            -math::grad(p) + diffusion,
        rAU, relaxed(state.m_control.velocity_relaxation));
    state.m_step = State::Step::Momentum;
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
