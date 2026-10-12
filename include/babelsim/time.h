#pragma once

#include <cmath>
#include <stdexcept>

namespace babelsim {

// Physical time configuration; histories are explicit time::History values.
struct TimeControl {
    double start_time = 0.0;
    double end_time = 1.0;
    double delta_t = 1e-3;

    void validate() const {
        if (!std::isfinite(start_time) || !std::isfinite(end_time) ||
            !std::isfinite(delta_t) || end_time < start_time || delta_t <= 0.0) {
            throw std::invalid_argument("time control is invalid");
        }
    }
};

class Case;
class Mesh;
namespace time { class Step; }
// Explicit time service. Queries have no side effects; only advance() changes time.
class TimeStepper {
public:
    explicit TimeStepper(Case&);
    TimeStepper(const TimeStepper&) = delete;
    TimeStepper& operator=(const TimeStepper&) = delete;
    void advance();
    // An explicitly requested step is clipped at endTime. It does not change
    // the configured default dt used by subsequent advance() calls.
    void advance(double requestedDt);
    double value() const { return value_; }
    double end() const { return options_.end_time; }
    double dt() const { return dt_; }
    int step() const { return step_; }
    bool finished() const { return value_ >= options_.end_time; }
    const Case& owner() const {return *case_;}
    const Mesh& mesh() const;
private:
    void advanceImpl(double requestedDt, bool explicitStep);
    Case* case_;
    TimeControl options_;
    double value_, dt_;
    int step_=0;
    bool variable_steps_=false;
    bool attempt_active_=false;
    friend class time::Step;
};
namespace time {
TimeStepper makeStepper(Case&);
}

}  // babelsim 命名空间
