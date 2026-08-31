#include "state.h"

namespace babelsim {

SolveResult SimpleSolver::State::solveMomentum() {
    // UEqn：求预测速度及对角体积响应。rAU 是当前算法的数学场，
    // 其计算复用公开 solveWithResponse，不让算法接触离散矩阵。
    VectorField& U = m_U;
    ScalarField& p = m_p;
    ScalarField& phi = m_phi;
    ScalarField& rAU = m_algorithm.rAU;

    return solveWithResponse(
        eqn::div(m_fluid.density, phi, U) ==
            -math::grad(p) + eqn::laplacian(m_fluid.dynamic_viscosity, U),
        rAU, relaxed(m_control.velocity_relaxation));
}

void SimpleSolver::State::predictMomentumFlux() {
    // phiHbyA：同位网格中的动量插值/Rhie-Chow 重构。这里是 SIMPLE 的私有数值
    // 步骤；物理通量 phi 只在压力方程完成后更新。
    NumericalWorkspace& work = m_workspace;
    math::evaluate(math::grad(m_p), work.grad_p);
    math::evaluate(math::flux(m_U), m_algorithm.phiHbyA);
    work.rAU_grad_p.assignProduct(m_algorithm.rAU, work.grad_p);
    math::evaluate(math::interpolate(work.rAU_grad_p), work.rAU_grad_p_face);
    math::evaluate(math::interpolate(m_algorithm.rAU), work.rAU_face);
    // phiHbyA += Sf·interpolate(rAU grad(p)) - rAU_f Sf·grad(p)。
    // 只修正内部面，保持入口/壁面等物理边界通量；面选择与同步属于通用 math。
    math::add(math::flux(work.rAU_grad_p_face), m_algorithm.phiHbyA,
             math::FaceRegion::Interior);
    math::subtract(math::flux(work.rAU_face, math::reconstruct(m_p, work.grad_p)),
                  m_algorithm.phiHbyA, math::FaceRegion::Interior);
}

}  // babelsim 命名空间
