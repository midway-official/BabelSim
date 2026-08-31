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
        fvm::div(m_fluid.density, phi, U) ==
            -fvc::grad(p) + fvm::laplacian(m_fluid.dynamic_viscosity, U),
        rAU, relaxed(m_control.velocity_relaxation));
}

void SimpleSolver::State::predictMomentumFlux() {
    // phiHbyA：同位网格中的动量插值/Rhie-Chow 重构。这里是 SIMPLE 的私有数值
    // 步骤；物理通量 phi 只在压力方程完成后更新。
    NumericalWorkspace& work = m_workspace;
    fvc::evaluate(fvc::grad(m_p), work.grad_p);
    fvc::evaluate(fvc::flux(m_U), m_algorithm.phiHbyA);
    work.rAU_grad_p.assignProduct(m_algorithm.rAU, work.grad_p);
    fvc::evaluate(fvc::interpolate(work.rAU_grad_p), work.rAU_grad_p_face);
    fvc::evaluate(fvc::interpolate(m_algorithm.rAU), work.rAU_face);
    // phiHbyA += Sf·interpolate(rAU grad(p)) - rAU_f Sf·grad(p)。
    // 只修正内部面，保持入口/壁面等物理边界通量；面选择与同步属于通用 fvc。
    fvc::add(fvc::flux(work.rAU_grad_p_face), m_algorithm.phiHbyA,
             fvc::FaceRegion::Interior);
    fvc::subtract(fvc::flux(work.rAU_face, fvc::reconstruct(m_p, work.grad_p)),
                  m_algorithm.phiHbyA, fvc::FaceRegion::Interior);
}

}  // babelsim 命名空间
