#include "babelsim/incompressible.h"

#include "internal/scalar_equation_control.h"

#include <cmath>

namespace babelsim {
namespace {

bool healthy(const SolveResult& result) {
    return result.status != SolveStatus::NumericalFailure &&
        std::isfinite(result.final_residual) && std::isfinite(result.relative_residual);
}

// 压力方程的第一次求解给出正交部分，后续求解用于显式非正交修正。调用点只表达
// “继续修正非正交项”，不再暴露 pass 编号或首末次分支。
class NonOrthogonalCorrections {
public:
    NonOrthogonalCorrections(DiffusionMethod method, int corrections)
        : m_remaining(method == DiffusionMethod::Orthogonal ? 1 : corrections + 1)
    {}

    bool correctNonOrthogonal() {
        if (m_remaining == 0) return false;
        --m_remaining;
        return true;
    }

private:
    int m_remaining;
};

}  // 匿名命名空间

void SimpleSolver::initializePressureCorrectionBoundaries() {
    for (Index patch = 0; patch < static_cast<Index>(m_mesh.patches.size()); ++patch) {
        const BoundaryType pressure_type = m_fields.pressure.boundary(patch).type;
        if (pressure_type == BoundaryType::FixedValue) {
            m_has_fixed_pressure = true;
            m_algorithm.p_prime.setBoundary(patch, fixedValue(0.0));
        } else if (pressure_type == BoundaryType::Symmetry) {
            m_algorithm.p_prime.setBoundary(patch, symmetry());
        } else {
            m_algorithm.p_prime.setBoundary(patch, zeroGradient());
        }
    }
}

SimpleSolver::PressureEquationResult SimpleSolver::solvePressure() {
    // pEqn：div(phiHbyA) = div(rAU grad(pPrime))。非正交显式修正由通用
    // fvm::laplacian 处理；这里只保留 SIMPLE 所需的语义化修正循环。
    ScalarField& p = m_fields.pressure;
    ScalarField& pPrime = m_algorithm.p_prime;
    ScalarField& phiHbyA = m_algorithm.phiHbyA;
    ScalarField& rAUFace = m_workspace.rAU_face;
    ScalarField& divPhiHbyA = m_workspace.div_phiHbyA;

    predictMomentumFlux();
    fvc::evaluate(fvc::div(phiHbyA), divPhiHbyA);
    pPrime.fill(0.0);

    PressureEquationResult result;
    result.healthy = true;
    result.linear_converged = true;
    ScalarEquationControl equation_control;
    equation_control.fix_reference = !m_has_fixed_pressure;
    NonOrthogonalCorrections corrections(
        m_methods.diffusion, m_control.non_orthogonal_corrections);
    while (corrections.correctNonOrthogonal()) {
        result.linear = detail::solve(
            -fvm::laplacian(rAUFace, pPrime) == -fvm::source(divPhiHbyA),
            equation_control);
        result.healthy = result.healthy && healthy(result.linear);
        result.linear_converged = result.linear_converged && result.linear.converged();
    }

    for (Index cell : m_mesh.owned_cells) {
        p[cell] += m_control.pressure_relaxation * pPrime[cell];
    }
    return result;
}

void SimpleSolver::correctVelocity() {
    VectorField& U = m_fields.velocity;
    const ScalarField& pPrime = m_algorithm.p_prime;
    const ScalarField& rAU = m_algorithm.rAU;
    VectorField& gradPPrime = m_workspace.grad_p_prime;

    fvc::evaluate(fvc::grad(pPrime), gradPPrime);
    for (Index cell : m_mesh.owned_cells) {
        U[cell] -= rAU[cell] * gradPPrime[cell];
    }
}

void SimpleSolver::correctFlux() {
    ScalarField& phi = m_fields.face_flux;
    const ScalarField& phiHbyA = m_algorithm.phiHbyA;
    const ScalarField& pPrime = m_algorithm.p_prime;
    const ScalarField& rAUFace = m_workspace.rAU_face;
    const VectorField& gradPPrime = m_workspace.grad_p_prime;

    phi.assign(phiHbyA);
    for (Index face : m_mesh.owned_faces) {
        const std::size_t index = static_cast<std::size_t>(face);
        const Index neighbour = m_mesh.face_neighbour[index];
        const bool fixed_pressure = neighbour == invalid_index &&
            m_fields.pressure.boundary(m_mesh.face_patch[index]).type ==
                BoundaryType::FixedValue;
        if (neighbour != invalid_index || fixed_pressure) {
            phi[face] -= rAUFace[face] * fvc::integratedNormalGradient(
                pPrime, gradPPrime, face, m_methods.diffusion);
        }
    }
}

}  // babelsim 命名空间
