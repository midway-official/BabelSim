#include "state.h"
#include "babelsim/case.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace babelsim {

TransientSimpleAlgorithm::State::State(VectorField& U, ScalarField& p, ScalarField& phi,
                                       FluidProperties fluid, SimpleControl control)
    : m_U(U), m_p(p), m_phi(phi), m_mesh(U.mesh()),
      m_fluid(fluid), m_control(control), m_methods(numericalMethods())
{
    if (&m_p.mesh() != &m_mesh || &m_phi.mesh() != &m_mesh ||
        m_U.location() != FieldLocation::Cell ||
        m_p.location() != FieldLocation::Cell || m_phi.location() != FieldLocation::Face) {
        throw std::invalid_argument("incompressible fields do not match the run time mesh");
    }
    if (m_methods.time == TimeMethod::Steady) {
        throw std::invalid_argument("transient SIMPLE requires Euler or BDF2 time discretization");
    }
    m_fluid.validate();
    m_control.validate();
    m_has_fixed_pressure = setHomogeneousCorrectionBoundaries(m_p_prime, m_p);
    math::evaluate(math::flux(m_U), m_phi);
}

TransientSimpleAlgorithm::TransientSimpleAlgorithm(
    IncompressibleFields& fields, FluidProperties fluid, SimpleControl control)
    : m_state(std::make_unique<State>(fields.velocity, fields.pressure,
                                     fields.face_flux, fluid, control))
{}

TransientSimpleAlgorithm::TransientSimpleAlgorithm(Case& problem)
    : m_state(std::make_unique<State>(
          problem.vectorField("U"), problem.scalarField("p"), problem.faceField("phi"),
          FluidProperties{problem.physics().positive("density"),
                          problem.physics().positive("dynamicViscosity")},
          readSimpleControl(problem.solution())))
{
    m_state->m_log = true;
}

TransientSimpleAlgorithm::~TransientSimpleAlgorithm() = default;

void TransientSimpleAlgorithm::State::requireStep(Step expected) const {
    if (m_step != expected) throw std::logic_error("transient SIMPLE steps called out of order");
}

void TransientSimpleAlgorithm::beginTimeStep() {
    State& state = *m_state;
    if (state.m_step != State::Step::Ready && state.m_step != State::Step::Complete)
        throw std::logic_error("cannot begin a time step during an incomplete SIMPLE correction");
    state.m_result = {};
    state.m_pressure_healthy = false;
    state.m_pressure_converged = false;
    state.m_step = State::Step::Ready;
    state.m_iteration = 0;
}

bool TransientSimpleAlgorithm::loop() const {
    const State& state = *m_state;
    if (state.m_iteration < 0)
        throw std::logic_error("call beginTimeStep after Case::loop before transient SIMPLE corrections");
    if (state.m_step != State::Step::Ready && state.m_step != State::Step::Complete)
        throw std::logic_error("transient SIMPLE correction has not completed");
    return state.m_iteration < state.m_control.max_iterations &&
           (state.m_iteration == 0 || (state.m_result.healthy && !state.m_result.converged));
}

void TransientSimpleAlgorithm::checkContinuity() {
    State& state = *m_state;
    state.requireStep(State::Step::Flux);
    SimpleIterationResult& result = state.m_result;
    result.relative_velocity_change = diagnostics::relativeChange(state.m_U, state.m_previous_velocity);
    result.relative_pressure_correction = diagnostics::relativeMagnitude(state.m_p_prime, state.m_p);
    result.continuity = diagnostics::fluxBalance(state.m_phi);

    result.healthy = diagnostics::all(state.m_pressure_healthy &&
        std::isfinite(result.relative_velocity_change) &&
        std::isfinite(result.continuity.relative) && result.velocity.healthy());
    result.linear_converged = diagnostics::all(
        state.m_pressure_converged && result.velocity.converged());
    result.converged = result.healthy && result.linear_converged && diagnostics::all(
        result.continuity.relative <= state.m_control.continuity_tolerance &&
        result.relative_velocity_change <= state.m_control.velocity_tolerance);
    state.m_step = State::Step::Complete;
    ++state.m_iteration;
    state.report();
}

bool TransientSimpleAlgorithm::converged() const { return m_state->m_result.converged; }

SimpleIterationResult TransientSimpleAlgorithm::iterate() {
    solveMomentum();
    solvePressure();
    correctVelocity();
    correctFlux();
    checkContinuity();
    return m_state->m_result;
}

void TransientSimpleAlgorithm::State::report() const {
    if (!m_log) return;
    if (m_iteration != 1 && m_iteration % 100 != 0 && !m_result.converged &&
        m_result.healthy && m_iteration != m_control.max_iterations) return;
    std::ostringstream message;
    message << "Transient SIMPLE " << m_iteration
            << std::scientific << std::setprecision(6)
            << " mass=" << m_result.continuity.relative
            << " dU=" << m_result.relative_velocity_change
            << " linP=" << m_result.pressure.relative_residual
            << " linear=" << (m_result.linear_converged ? "ok" : "inexact")
            << " converged=" << (m_result.converged ? "true" : "false");
    diagnostics::report(message.str());
}

}  // babelsim 命名空间
