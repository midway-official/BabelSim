#include "babelsim/time_step.h"
#include "babelsim/case.h"
#include <exception>

namespace babelsim::time {
Step::Step(TimeStepper& clock,HistoryGroup& histories,double requestedDt)
    : clock_(&clock),histories_(&histories),value_(clock.value_),dt_(clock.dt_),
      step_(clock.step_),variable_(clock.variable_steps_) {
    if (clock.attempt_active_ || histories.active_)
        throw std::logic_error("cannot nest pending physical steps");
    snapshot_=histories.checkpoint();
    try {
        if (requestedDt==0) clock.advance(); else clock.advance(requestedDt);
        histories.save(clock);
    } catch (...) {
        histories.restore(snapshot_);
        if (clock.step_!=step_) restoreClock();
        throw;
    }
    clock.attempt_active_=true; histories.active_=true; pending_=true;
}
Step::Step(Step&& other) noexcept
    : clock_(other.clock_),histories_(other.histories_),snapshot_(std::move(other.snapshot_)),
      value_(other.value_),dt_(other.dt_),step_(other.step_),variable_(other.variable_),pending_(other.pending_) {
    other.pending_=false;
}
Step::~Step() noexcept {
    if (pending_) { try { reject(); } catch (...) { std::terminate(); } }
}
void Step::restoreClock() {
    clock_->case_->setTime(value_,step_,dt_);
    clock_->value_=value_; clock_->dt_=dt_; clock_->step_=step_; clock_->variable_steps_=variable_;
}
void Step::accept() {
    if (!pending_) throw std::logic_error("physical step is already resolved");
    histories_->validate();
    if (clock_->step_!=step_+1 || histories_->lastStep_!=clock_->step_
        || clock_->case_->step()!=clock_->step_ || clock_->case_->time()!=clock_->value_)
        throw std::logic_error("clock or history checkpoint changed during pending step");
    clock_->attempt_active_=false; histories_->active_=false; pending_=false;
}
void Step::reject() {
    if (!pending_) throw std::logic_error("physical step is already resolved");
    histories_->restore(snapshot_); restoreClock();
    clock_->attempt_active_=false; histories_->active_=false; pending_=false;
}
}
