#include "state.h"

namespace babelsim {

void SteadySimpleAlgorithm::solveMomentum() {
    State& state = *m_state;
    if (state.m_step == State::Step::Complete) state.m_step = State::Step::Ready;
    state.requireStep(State::Step::Ready);
    state.m_previous_velocity.assign(state.m_U);
    // UEqn：求预测速度及对角体积响应。rAU 是当前算法的数学场，
    // 其计算复用公开 solveWithResponse，不让算法接触离散矩阵。
    state.m_result.velocity = solveWithResponse(state.momentumEquation(), state.m_rAU,
        relaxed(state.m_control.velocity_relaxation));
    state.m_step = State::Step::Momentum;
}

VectorEquationDefinition SteadySimpleAlgorithm::State::momentumEquation() {
    VectorExpression diffusion = m_turbulence
        ? eqn::laplacian(m_effective_viscosity, m_U)
        : eqn::laplacian(m_fluid.dynamic_viscosity, m_U);
    if (m_turbulence) {
        evaluateStressCorrection(m_U, m_phi, m_effective_viscosity,
            m_stress_gradient, m_stress_correction, m_stress_divergence);
        diffusion = diffusion + eqn::source(m_stress_divergence);
    }
    return eqn::div(m_fluid.density, m_phi, m_U) == -math::grad(m_p) + diffusion;
}

void SteadySimpleAlgorithm::State::predictMomentumFlux() {
    // phiHbyA：同位网格中的动量插值/Rhie-Chow 重构。这里是 SIMPLE 的私有数值
    // 步骤；物理通量 phi 只在压力方程完成后更新。
    math::evaluate(math::grad(m_p), m_grad_p);
    math::evaluate(math::flux(m_U), m_phiHbyA);
    m_rAU_grad_p.assignProduct(m_rAU, m_grad_p);
    math::evaluate(math::interpolate(m_rAU_grad_p), m_rAU_grad_p_face);
    math::evaluate(math::interpolate(m_rAU), m_rAU_face);
    // phiHbyA += Sf·interpolate(rAU grad(p)) - rAU_f Sf·grad(p)。
    // 只修正内部面，保持入口/壁面等物理边界通量；面选择与同步属于通用 math。
    math::add(math::flux(m_rAU_grad_p_face), m_phiHbyA, math::FaceRegion::Interior);
    math::subtract(math::flux(m_rAU_face, math::reconstruct(m_p, m_grad_p)),
                  m_phiHbyA, math::FaceRegion::Interior);
}

}  // babelsim 命名空间
