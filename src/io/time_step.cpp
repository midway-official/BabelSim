#include "babelsim/time_step.h"
#include "babelsim/case.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>

namespace babelsim {

TimeStepper::TimeStepper(Case& problem)
    : case_(&problem), options_(problem.read.timeControl()), value_(options_.start_time),
      dt_(options_.delta_t) {
    options_.validate();
}
TimeStepper time::makeStepper(Case& problem) {
    return TimeStepper(problem);
}
const Mesh& TimeStepper::mesh() const {
    return case_->mesh();
}
void TimeStepper::advance() {
    advanceImpl(options_.delta_t, false);
}
void TimeStepper::advance(double requestedDt) {
    advanceImpl(requestedDt, true);
}
void TimeStepper::advanceImpl(double requestedDt, bool explicitStep) {
    auto& time = *this;
    if (time.attempt_active_)
        throw std::logic_error("resolve the pending step before advancing the clock");
    if (!(requestedDt > 0) || !std::isfinite(requestedDt))
        throw std::invalid_argument("requested step requires positive finite dt");
    if (time.finished())
        throw std::logic_error("cannot advance beyond endTime");
    const double tolerance = 32 * std::numeric_limits<double>::epsilon() *
                             std::max({std::abs(time.options_.start_time),
                                       std::abs(time.options_.end_time), time.options_.delta_t});
    const bool variable = time.variable_steps_ || explicitStep;
    double next =
        variable
            ? static_cast<double>(static_cast<long double>(time.value_) + requestedDt)
            : static_cast<double>(static_cast<long double>(time.options_.start_time) +
                                  static_cast<long double>(time.step_ + 1) * time.options_.delta_t);
    if (next >= time.options_.end_time - tolerance)
        next = time.options_.end_time;
    const double dt = next - time.value_;
    if (!(dt > 0))
        throw std::runtime_error("time step is below representable time precision");
    time.case_->setTime(next, time.step_ + 1, dt);
    time.value_ = next;
    time.dt_ = dt;
    time.variable_steps_ = variable;
    ++time.step_;
}
void write(Case& problem, const TimeStepper& time) {
    if (&problem != &time.owner())
        throw std::invalid_argument("time service belongs to another case");
    problem.setTime(time.value(), time.step(), time.dt());
    problem.write();
}

TimeOptions readTimeControl(const Case& problem) {
    const auto& c = problem.read.timeControl();
    return {c.start_time, c.end_time, c.delta_t};
}
} // namespace babelsim

namespace babelsim::time {
Step::Step(TimeStepper& clock, HistoryGroup& histories, double requestedDt)
    : clock_(&clock), histories_(&histories), value_(clock.value_), dt_(clock.dt_),
      step_(clock.step_), variable_(clock.variable_steps_) {
    if (clock.attempt_active_ || histories.active_)
        throw std::logic_error("cannot nest pending physical steps");
    snapshot_ = histories.checkpoint();
    try {
        if (requestedDt == 0)
            clock.advance();
        else
            clock.advance(requestedDt);
        histories.save(clock);
    } catch (...) {
        histories.restore(snapshot_);
        if (clock.step_ != step_)
            restoreClock();
        throw;
    }
    clock.attempt_active_ = true;
    histories.active_ = true;
    pending_ = true;
}
Step::Step(Step&& other) noexcept
    : clock_(other.clock_), histories_(other.histories_), snapshot_(std::move(other.snapshot_)),
      value_(other.value_), dt_(other.dt_), step_(other.step_), variable_(other.variable_),
      pending_(other.pending_) {
    other.pending_ = false;
}
Step::~Step() noexcept {
    if (pending_) {
        try {
            reject();
        } catch (...) {
            std::terminate();
        }
    }
}
void Step::restoreClock() {
    clock_->case_->setTime(value_, step_, dt_);
    clock_->value_ = value_;
    clock_->dt_ = dt_;
    clock_->step_ = step_;
    clock_->variable_steps_ = variable_;
}
void Step::accept() {
    if (!pending_)
        throw std::logic_error("physical step is already resolved");
    histories_->validate();
    if (clock_->step_ != step_ + 1 || histories_->lastStep_ != clock_->step_ ||
        clock_->case_->step() != clock_->step_ || clock_->case_->time() != clock_->value_)
        throw std::logic_error("clock or history checkpoint changed during pending step");
    clock_->attempt_active_ = false;
    histories_->active_ = false;
    pending_ = false;
}
void Step::reject() {
    if (!pending_)
        throw std::logic_error("physical step is already resolved");
    histories_->restore(snapshot_);
    restoreClock();
    clock_->attempt_active_ = false;
    histories_->active_ = false;
    pending_ = false;
}
}
