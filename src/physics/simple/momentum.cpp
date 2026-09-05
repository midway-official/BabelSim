#include "state.h"

namespace babelsim {

void SteadySimpleAlgorithm::solveMomentum() {
    State& state = *m_state;
    if (state.m_step == State::Step::Complete) state.m_step = State::Step::Ready;
    state.requireStep(State::Step::Ready);
    state.m_previous_velocity.assign(state.m_U);
    // UEqn：求预测速度及对角体积响应。rAU 是当前算法的数学场，
    // 其计算复用公开 solveWithResponse，不让算法接触离散矩阵。
    VectorField& U = state.m_U;
    ScalarField& p = state.m_p;
    ScalarField& phi = state.m_phi;
    ScalarField& rAU = state.m_rAU;
    const VectorExpression diffusion = state.m_turbulence
        ? eqn::laplacian(state.m_effective_viscosity, U)
        : eqn::laplacian(state.m_fluid.dynamic_viscosity, U);

    state.m_result.velocity = solveWithResponse(
        eqn::div(state.m_fluid.density, phi, U) ==
            -math::grad(p) + diffusion,
        rAU, relaxed(state.m_control.velocity_relaxation));
    state.m_step = State::Step::Momentum;
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
