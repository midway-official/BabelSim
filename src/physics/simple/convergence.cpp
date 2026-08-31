#include "state.h"
#include "babelsim/runtime.h"

#include <cmath>
#include <iomanip>
#include <iostream>

namespace babelsim {
void SimpleSolver::State::checkContinuityAndConvergence(
    SimpleIterationResult& result,
    const PressureEquationResult& pressure) const
{
    result.pressure = pressure.linear;
    result.relative_velocity_change = diagnostics::relativeChange(
        m_U, m_algorithm.previous_velocity);
    result.relative_pressure_correction = diagnostics::relativeMagnitude(
        m_algorithm.p_prime, m_p);
    result.continuity = diagnostics::fluxBalance(m_phi);

    bool local_healthy = pressure.healthy &&
        std::isfinite(result.relative_velocity_change) &&
        std::isfinite(result.continuity.relative);
    bool local_linear_converged = pressure.linear_converged;
    local_healthy = local_healthy && result.velocity.healthy();
    local_linear_converged = local_linear_converged && result.velocity.converged();
    result.healthy = diagnostics::all(local_healthy);
    result.linear_converged = diagnostics::all(local_linear_converged);
    result.converged = result.healthy && result.linear_converged && diagnostics::all(
        result.continuity.relative <= m_control.continuity_tolerance &&
        result.relative_velocity_change <= m_control.velocity_tolerance);
}

void SimpleSolver::State::report() const {
    if (!m_log || !RunTime::current().primary()) return;
    if (m_iteration != 1 && m_iteration % 100 != 0 && !m_result.converged &&
        m_result.healthy && m_iteration != m_control.max_iterations) return;
    std::cout << "SIMPLE " << m_iteration << std::scientific << std::setprecision(6)
              << " mass=" << m_result.continuity.relative
              << " dU=" << m_result.relative_velocity_change
              << " linP=" << m_result.pressure.relative_residual
              << " linear=" << (m_result.linear_converged ? "ok" : "inexact")
              << " converged=" << (m_result.converged ? "true" : "false") << '\n';
}

}  // babelsim 命名空间
