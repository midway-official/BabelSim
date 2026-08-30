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
    : m_run_time(run_time),
      m_fields(fields),
      m_mesh(fields.velocity.mesh()),
      m_fluid(fluid),
      m_control(control),
      m_pressure_correction(m_mesh, FieldLocation::Cell, "pPrime"),
      m_pressure_gradient(m_mesh, FieldLocation::Cell, "gradP"),
      m_correction_gradient(m_mesh, FieldLocation::Cell, "gradPPrime"),
      m_mass_flux(m_mesh, FieldLocation::Face, "rhoPhi"),
      m_mobility(m_mesh, FieldLocation::Cell, "rAU"),
      m_face_mobility(m_mesh, FieldLocation::Face, "rAUFace"),
      m_divergence(m_mesh, FieldLocation::Cell, "divPhi"),
      m_previous_velocity(m_mesh, FieldLocation::Cell, "UPrevious"),
      m_interpolation_workspace(m_mesh)
{
    if (&m_run_time.mesh() != &m_mesh || &m_fields.pressure.mesh() != &m_mesh ||
        &m_fields.face_flux.mesh() != &m_mesh ||
        m_fields.velocity.location() != FieldLocation::Cell ||
        m_fields.pressure.location() != FieldLocation::Cell ||
        m_fields.face_flux.location() != FieldLocation::Face) {
        throw std::invalid_argument("incompressible fields do not match the run time mesh");
    }
    if (m_run_time.methods().time != TimeMethod::Steady) {
        throw std::invalid_argument(
            "SimpleSolver is steady; use a transient pressure-velocity algorithm for non-steady time methods");
    }
    m_fluid.validate();
    m_control.validate();
    for (Index patch = 0; patch < static_cast<Index>(m_mesh.patches.size()); ++patch) {
        const BoundaryType pressure_type = m_fields.pressure.boundary(patch).type;
        if (pressure_type == BoundaryType::FixedValue) {
            m_has_fixed_pressure = true;
            m_pressure_correction.setBoundary(
                patch, BoundaryCondition<double>::fixedValue(0.0));
        } else if (pressure_type == BoundaryType::Symmetry) {
            m_pressure_correction.setBoundary(
                patch, BoundaryCondition<double>::symmetry());
        } else {
            m_pressure_correction.setBoundary(
                patch, BoundaryCondition<double>::zeroGradient());
        }
    }
    m_run_time.evaluate(fvc::flux(m_fields.velocity), m_fields.face_flux);
}

SimpleIterationResult SimpleSolver::iterate() {
    SimpleIterationResult result;

    // UEqn: div(rho*phi,U) = -grad(p) + laplacian(mu,U)。RunTime 在此处才
    // 执行同步、离散装配与线性求解；SIMPLE 只保留算法语义。
    m_run_time.copy(m_fields.velocity, m_previous_velocity);
    m_run_time.scale(m_fluid.density, m_fields.face_flux, m_mass_flux);
    VectorEquationControl momentum_control;
    momentum_control.relaxation = m_control.velocity_relaxation;
    momentum_control.mobility = &m_mobility;
    result.velocity = solve(
        m_run_time,
        fvm::div(m_mass_flux, m_fields.velocity) ==
            -fvc::grad(m_fields.pressure) +
            fvm::laplacian(m_fluid.dynamic_viscosity, m_fields.velocity),
        momentum_control);

    m_run_time.evaluate(fvc::grad(m_fields.pressure), m_pressure_gradient);
    MomentumInterpolation::apply(
        m_run_time, m_mesh, m_fields.velocity, m_fields.pressure,
        m_mobility, m_pressure_gradient, m_fields.face_flux,
        m_interpolation_workspace, m_run_time.methods().interpolation,
        m_run_time.methods().gradient, m_run_time.methods().diffusion);

    // pEqn: rAU 的面插值与 div(phi) 形成压力修正方程。修正扩散时重复该方程，
    // 以当前 pPrime 作为显式非正交项，和 OpenFOAM 的 correctNonOrthogonal 循环一致。
    const int pressure_passes =
        m_run_time.methods().diffusion == DiffusionMethod::Orthogonal
        ? 1 : m_control.non_orthogonal_corrections + 1;
    bool pressure_healthy = true;
    bool pressure_linear_converged = true;
    for (int pass = 0; pass < pressure_passes; ++pass) {
        result.pressure = PressureCorrection::solve(
            m_run_time, m_pressure_correction, m_mesh,
            m_fields.face_flux, m_mobility, m_fields.pressure,
            m_face_mobility, m_divergence, m_has_fixed_pressure,
            pass == 0);
        pressure_healthy = pressure_healthy && healthy(result.pressure);
        pressure_linear_converged =
            pressure_linear_converged && result.pressure.converged();
    }

    PressureCorrection::apply(
        m_run_time, m_mesh, m_control.pressure_relaxation,
        m_fields.pressure, m_fields.velocity, m_fields.face_flux,
        m_pressure_correction, m_mobility, m_correction_gradient,
        m_run_time.methods().diffusion);

    result.relative_velocity_change = m_run_time.relativeChange(
        m_fields.velocity, m_previous_velocity);
    result.relative_pressure_correction = m_run_time.relativeChange(
        m_pressure_correction, m_fields.pressure);
    result.continuity = m_run_time.fluxBalance(m_fields.face_flux);

    bool local_healthy = pressure_healthy &&
        std::isfinite(result.relative_velocity_change) &&
        std::isfinite(result.continuity.relative);
    bool local_linear_converged = pressure_linear_converged;
    for (const SolveResult& velocity_result : result.velocity) {
        local_healthy = local_healthy && healthy(velocity_result);
        local_linear_converged =
            local_linear_converged && velocity_result.converged();
    }
    result.healthy = m_run_time.all(local_healthy);
    result.linear_converged = m_run_time.all(local_linear_converged);
    result.converged = result.healthy && m_run_time.all(
        result.continuity.relative <= m_control.continuity_tolerance &&
        result.relative_velocity_change <= m_control.velocity_tolerance);
    return result;
}

}  // babelsim 命名空间
