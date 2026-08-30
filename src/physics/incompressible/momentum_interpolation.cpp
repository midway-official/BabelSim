#include "babelsim/incompressible.h"

#include "internal/vector_equation_control.h"

namespace babelsim {

std::array<SolveResult, 3> SimpleSolver::solveMomentum() {
    // UEqn：求预测速度，同时从离散主对角提取 rAU。rAU 仅属于 SIMPLE 的算法状态，
    // 不进入通用 Field 或 fvm API。
    m_workspace.previous_velocity.assign(m_fields.velocity);
    m_workspace.mass_flux.assignScaled(m_fluid.density, m_fields.face_flux);
    VectorEquationControl control;
    control.relaxation = m_control.velocity_relaxation;
    control.mobility = &m_workspace.mobility;
    const std::array<SolveResult, 3> result = RunTime::current().solve(
        fvm::div(m_workspace.mass_flux, m_fields.velocity) ==
            -fvc::grad(m_fields.pressure) +
            fvm::laplacian(m_fluid.dynamic_viscosity, m_fields.velocity),
        control);

    // phiHbyA：同位网格中的动量插值/Rhie-Chow 重构。这里是 SIMPLE 的私有数值
    // 步骤；SIMPLE 主循环只把它理解为“动量方程产生预测通量”。
    fvc::evaluate(fvc::grad(m_fields.pressure), m_workspace.pressure_gradient);
    fvc::evaluate(fvc::flux(m_fields.velocity), m_fields.face_flux);
    for (Index cell : m_mesh.owned_cells) {
        m_workspace.pressure_response[cell] =
            m_workspace.mobility[cell] * m_workspace.pressure_gradient[cell];
    }
    fvc::evaluate(
        fvc::interpolate(m_workspace.pressure_response),
        m_workspace.face_pressure_response);
    fvc::evaluate(
        fvc::interpolate(m_workspace.mobility), m_workspace.interpolation_mobility);
    for (Index face : m_mesh.owned_faces) {
        const std::size_t index = static_cast<std::size_t>(face);
        if (m_mesh.face_neighbour[index] == invalid_index) continue;
        m_fields.face_flux[face] +=
            dot(m_workspace.face_pressure_response[face], m_mesh.face_area_vectors[index]) -
            m_workspace.interpolation_mobility[face] * fvc::integratedNormalGradient(
                m_fields.pressure, m_workspace.pressure_gradient, face, m_methods.diffusion);
    }
    return result;
}

}  // babelsim 命名空间
