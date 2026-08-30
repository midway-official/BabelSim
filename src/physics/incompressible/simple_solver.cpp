#include "babelsim/incompressible.h"

#include <cmath>
#include <stdexcept>

namespace babelsim {
namespace {

bool finitePositive(double value) {
    return value > 0.0 && std::isfinite(value);
}

bool healthy(const SolveResult& result) {
    return result.status != SolveStatus::NumericalFailure &&
        std::isfinite(result.final_residual) &&
        std::isfinite(result.relative_residual);
}

}  // 匿名命名空间

void FluidProperties::validate() const {
    if (!finitePositive(density) || !finitePositive(dynamic_viscosity)) {
        throw std::invalid_argument("density and dynamic viscosity must be positive");
    }
}

void SimpleControl::validate() const {
    if (max_iterations <= 0 || non_orthogonal_corrections < 0 ||
        non_orthogonal_corrections > 20 ||
        !(velocity_relaxation > 0.0 && velocity_relaxation <= 1.0) ||
        !(pressure_relaxation > 0.0 && pressure_relaxation <= 1.0) ||
        !finitePositive(continuity_tolerance) || !finitePositive(velocity_tolerance)) {
        throw std::invalid_argument("SIMPLE controls are invalid");
    }
    velocity_solver.validate();
    pressure_solver.validate();
}

SimpleSolver::SimpleSolver(
    RunTime& run_time,
    IncompressibleFields& fields,
    FluidProperties fluid,
    SimpleControl control)
    : m_fields(fields),
      m_mesh(fields.velocity.mesh()),
      m_fluid(fluid),
      m_control(control),
      m_methods(run_time.methods()),
      m_algorithm(m_mesh),
      m_workspace(m_mesh)
{
    if (&run_time.mesh() != &m_mesh || &m_fields.pressure.mesh() != &m_mesh ||
        &m_fields.face_flux.mesh() != &m_mesh ||
        m_fields.velocity.location() != FieldLocation::Cell ||
        m_fields.pressure.location() != FieldLocation::Cell ||
        m_fields.face_flux.location() != FieldLocation::Face) {
        throw std::invalid_argument("incompressible fields do not match the run time mesh");
    }
    if (m_methods.time != TimeMethod::Steady) {
        throw std::invalid_argument(
            "SimpleSolver is steady; use a transient pressure-velocity algorithm for non-steady time methods");
    }
    m_fluid.validate();
    m_control.validate();
    initializePressureCorrectionBoundaries();
    fvc::evaluate(fvc::flux(m_fields.velocity), m_fields.face_flux);
}

SimpleIterationResult SimpleSolver::iterate() {
    SimpleIterationResult result;

    // 该顺序就是 SIMPLE：UEqn 生成 rAU/预测通量，pEqn 在非正交循环中求 pPrime，
    // 随后分别修正 U 和 phi，最后用全局连续性与速度变化判定外迭代。
    result.velocity = solveMomentum();
    const PressureEquationResult pressure = solvePressure();
    correctVelocity();
    correctFlux();
    checkContinuityAndConvergence(result, pressure);
    return result;
}

void SimpleSolver::checkContinuityAndConvergence(
    SimpleIterationResult& result,
    const PressureEquationResult& pressure) const
{
    result.pressure = pressure.linear;
    result.relative_velocity_change = diagnostics::relativeChange(
        m_fields.velocity, m_algorithm.previous_velocity);
    result.relative_pressure_correction = diagnostics::relativeMagnitude(
        m_algorithm.p_prime, m_fields.pressure);
    result.continuity = diagnostics::fluxBalance(m_fields.face_flux);

    bool local_healthy = pressure.healthy &&
        std::isfinite(result.relative_velocity_change) &&
        std::isfinite(result.continuity.relative);
    bool local_linear_converged = pressure.linear_converged;
    for (const SolveResult& velocity_result : result.velocity) {
        local_healthy = local_healthy && healthy(velocity_result);
        local_linear_converged =
            local_linear_converged && velocity_result.converged();
    }
    result.healthy = diagnostics::all(local_healthy);
    result.linear_converged = diagnostics::all(local_linear_converged);
    result.converged = result.healthy && result.linear_converged && diagnostics::all(
        result.continuity.relative <= m_control.continuity_tolerance &&
        result.relative_velocity_change <= m_control.velocity_tolerance);
}

}  // babelsim 命名空间
