#include "state.h"

#include <stdexcept>

namespace babelsim {

void SimpleSolver::State::requireStep(Step expected) const {
    if (m_step != expected) throw std::logic_error("SIMPLE steps called out of order");
}

bool SimpleSolver::loop() const {
    const State& state = *m_state;
    if (state.m_step != State::Step::Ready && state.m_step != State::Step::Complete)
        throw std::logic_error("SIMPLE iteration has not completed");
    return state.m_iteration < state.m_control.max_iterations &&
           (state.m_iteration == 0 || (state.m_result.healthy && !state.m_result.converged));
}

void SimpleSolver::solveMomentum() {
    State& state = *m_state;
    if (state.m_step == State::Step::Complete) state.m_step = State::Step::Ready;
    state.requireStep(State::Step::Ready);
    state.m_algorithm.previous_velocity.assign(state.m_U);
    state.m_result.velocity = state.solveMomentum();
    state.m_step = State::Step::Momentum;
}

void SimpleSolver::solvePressure() {
    m_state->requireStep(State::Step::Momentum);
    m_state->m_pressure_result = m_state->solvePressure();
    m_state->m_step = State::Step::Pressure;
}

void SimpleSolver::correctVelocity() {
    m_state->requireStep(State::Step::Pressure);
    m_state->correctVelocity();
    m_state->m_step = State::Step::Velocity;
}

void SimpleSolver::correctFlux() {
    m_state->requireStep(State::Step::Velocity);
    m_state->correctFlux();
    m_state->m_step = State::Step::Flux;
}

void SimpleSolver::checkContinuity() {
    State& state = *m_state;
    state.requireStep(State::Step::Flux);
    state.checkContinuityAndConvergence(state.m_result, state.m_pressure_result);
    state.m_step = State::Step::Complete;
    ++state.m_iteration;
    state.report();
}

bool SimpleSolver::converged() const { return m_state->m_result.converged; }

SimpleIterationResult SimpleSolver::iterate() {
    solveMomentum();
    solvePressure();
    correctVelocity();
    correctFlux();
    checkContinuity();
    return m_state->m_result;
}

}  // babelsim 命名空间
