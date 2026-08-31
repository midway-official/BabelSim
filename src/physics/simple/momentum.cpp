#include "state.h"

namespace babelsim {

std::array<SolveResult, 3> SimpleSolver::State::solveMomentum() {
    // UEqn：求预测速度，同时从离散主对角提取 rAU。rAU 仅属于 SIMPLE 的算法状态，
    // 不进入通用 Field 或 fvm API。
    VectorField& U = m_U;
    ScalarField& p = m_p;
    ScalarField& phi = m_phi;
    ScalarField& rAU = m_algorithm.rAU;

    return simple::solveMomentumEquation(
        fvm::div(m_fluid.density, phi, U) ==
            -fvc::grad(p) + fvm::laplacian(m_fluid.dynamic_viscosity, U),
        m_control.velocity_relaxation, rAU);
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
    detail::applyMomentumInterpolation(
        m_p, work.grad_p, work.rAU_face, work.rAU_grad_p_face,
        m_algorithm.phiHbyA, m_methods.diffusionFor(m_p.name()));
}

}  // babelsim 命名空间
