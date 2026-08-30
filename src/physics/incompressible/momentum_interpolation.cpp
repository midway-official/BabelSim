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

    m_algorithm.previous_velocity.assign(U);
    VectorEquationControl control;
    control.relaxation = m_control.velocity_relaxation;
    control.mobility = &rAU;
    return detail::solve(
        fvm::div(m_fluid.density, phi, U) ==
            -fvc::grad(p) + fvm::laplacian(m_fluid.dynamic_viscosity, U),
        control);
}

void SimpleSolver::predictMomentumFlux() {
    // phiHbyA：同位网格中的动量插值/Rhie-Chow 重构。这里是 SIMPLE 的私有数值
    // 步骤；物理通量 phi 只在压力方程完成后更新。
    const VectorField& U = m_fields.velocity;
    const ScalarField& p = m_fields.pressure;
    const ScalarField& rAU = m_algorithm.rAU;
    ScalarField& phiHbyA = m_algorithm.phiHbyA;

    fvc::evaluate(fvc::grad(p), m_workspace.grad_p);
    fvc::evaluate(fvc::flux(U), phiHbyA);
    for (Index cell : m_mesh.owned_cells) {
        m_workspace.rAU_grad_p[cell] =
            rAU[cell] * m_workspace.grad_p[cell];
    }
    fvc::evaluate(
        fvc::interpolate(m_workspace.rAU_grad_p),
        m_workspace.rAU_grad_p_face);
    fvc::evaluate(
        fvc::interpolate(rAU), m_workspace.rAU_face);
    for (Index face : m_mesh.owned_faces) {
        const std::size_t index = static_cast<std::size_t>(face);
        if (m_mesh.face_neighbour[index] == invalid_index) continue;
        phiHbyA[face] +=
            dot(m_workspace.rAU_grad_p_face[face], m_mesh.face_area_vectors[index]) -
            m_workspace.rAU_face[face] * fvc::integratedNormalGradient(
                p, m_workspace.grad_p, face, m_methods.diffusion);
    }
}

}  // babelsim 命名空间
