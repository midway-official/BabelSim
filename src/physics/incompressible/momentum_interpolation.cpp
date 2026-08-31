#include "babelsim/incompressible.h"

#include "internal/vector_equation_control.h"

namespace babelsim {

std::array<SolveResult, 3> SimpleSolver::solveMomentum() {
    // UEqn：求预测速度，同时从离散主对角提取 rAU。rAU 仅属于 SIMPLE 的算法状态，
    // 不进入通用 Field 或 fvm API。
    VectorField& U = m_fields.velocity;
    ScalarField& p = m_fields.pressure;
    ScalarField& phi = m_fields.face_flux;
    ScalarField& rAU = m_algorithm.rAU;

    VectorEquationControl control;
    control.relaxation = m_control.velocity_relaxation;
    control.mobility = &rAU;
    return detail::solve(
        fvm::div(m_fluid.density, phi, U) ==
            -fvc::grad(p) + fvm::laplacian(m_fluid.dynamic_viscosity, U),
        control);
}

void SimpleSolver::savePreviousState() {
    m_algorithm.previous_velocity.assign(m_fields.velocity);
}

void SimpleSolver::predictMomentumFlux() {
    // phiHbyA：同位网格中的动量插值/Rhie-Chow 重构。这里是 SIMPLE 的私有数值
    // 步骤；物理通量 phi 只在压力方程完成后更新。
    m_workspace.predictMomentumFlux(
        m_methods, m_fields.velocity, m_fields.pressure,
        m_algorithm.rAU, m_algorithm.phiHbyA);
}

}  // babelsim 命名空间
