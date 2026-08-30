#include "babelsim/incompressible.h"

#include "internal/scalar_equation_control.h"

#include <cmath>

namespace babelsim {
namespace {

bool healthy(const SolveResult& result) {
    return result.status != SolveStatus::NumericalFailure &&
        std::isfinite(result.final_residual) && std::isfinite(result.relative_residual);
}

}  // 匿名命名空间

void SimpleSolver::initializePressureCorrectionBoundaries() {
    for (Index patch = 0; patch < static_cast<Index>(m_mesh.patches.size()); ++patch) {
        const BoundaryType pressure_type = m_fields.pressure.boundary(patch).type;
        if (pressure_type == BoundaryType::FixedValue) {
            m_has_fixed_pressure = true;
            m_workspace.pressure_correction.setBoundary(patch, fixedValue(0.0));
        } else if (pressure_type == BoundaryType::Symmetry) {
            m_workspace.pressure_correction.setBoundary(patch, symmetry());
        } else {
            m_workspace.pressure_correction.setBoundary(patch, zeroGradient());
        }
    }
}

SimpleSolver::PressureStep SimpleSolver::solvePressure() {
    // pEqn：div(phiHbyA) = div(rAU grad(pPrime))。非正交显式修正由通用
    // fvm::laplacian 处理；这里只保留 SIMPLE 所需的 correctNonOrthogonal 循环。
    fvc::evaluate(fvc::interpolate(m_workspace.mobility), m_workspace.face_mobility);
    fvc::evaluate(fvc::div(m_fields.face_flux), m_workspace.divergence);
    m_workspace.pressure_correction.fill(0.0);

    const int passes = m_methods.diffusion == DiffusionMethod::Orthogonal
        ? 1
        : m_control.non_orthogonal_corrections + 1;
    PressureStep result;
    result.healthy = true;
    result.linear_converged = true;
    ScalarEquationControl equation_control;
    equation_control.fix_reference = !m_has_fixed_pressure;
    for (int pass = 0; pass < passes; ++pass) {
        result.linear = RunTime::current().solve(
            -fvm::laplacian(m_workspace.face_mobility, m_workspace.pressure_correction) ==
                -fvm::source(m_workspace.divergence),
            equation_control);
        result.healthy = result.healthy && healthy(result.linear);
        result.linear_converged = result.linear_converged && result.linear.converged();
    }
    return result;
}

void SimpleSolver::correctVelocity() {
    fvc::evaluate(
        fvc::grad(m_workspace.pressure_correction), m_workspace.correction_gradient);
    for (Index cell : m_mesh.owned_cells) {
        m_fields.pressure[cell] +=
            m_control.pressure_relaxation * m_workspace.pressure_correction[cell];
        m_fields.velocity[cell] -=
            m_workspace.mobility[cell] * m_workspace.correction_gradient[cell];
    }
}

void SimpleSolver::correctFlux() {
    for (Index face : m_mesh.owned_faces) {
        const std::size_t index = static_cast<std::size_t>(face);
        const Index owner = m_mesh.face_owner[index];
        const Index neighbour = m_mesh.face_neighbour[index];
        if (neighbour != invalid_index) {
            const double owner_weight = m_mesh.face_owner_weights[index];
            const double face_mobility = owner_weight * m_workspace.mobility[owner] +
                (1.0 - owner_weight) * m_workspace.mobility[neighbour];
            m_fields.face_flux[face] -= face_mobility * fvc::integratedNormalGradient(
                m_workspace.pressure_correction, m_workspace.correction_gradient,
                face, m_methods.diffusion);
        } else if (m_fields.pressure.boundary(m_mesh.face_patch[index]).type ==
                   BoundaryType::FixedValue) {
            m_fields.face_flux[face] -= m_workspace.mobility[owner] *
                fvc::integratedNormalGradient(
                    m_workspace.pressure_correction, m_workspace.correction_gradient,
                    face, m_methods.diffusion);
        }
    }
}

}  // babelsim 命名空间
