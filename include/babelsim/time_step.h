#pragma once
#include "babelsim/history.h"

namespace babelsim::time {
// Optional rollback transaction. Declare/watch ALL state before beginStep;
// output and other external side effects must occur only after accept().
class Step {
public:
    Step(TimeStepper&,HistoryGroup&,double requestedDt=0);
    ~Step() noexcept;
    Step(const Step&)=delete;
    Step& operator=(const Step&)=delete;
    Step(Step&&) noexcept;
    Step& operator=(Step&&)=delete;
    void accept();
    void reject();
    bool pending() const { return pending_; }
private:
    void restoreClock();
    TimeStepper* clock_;
    HistoryGroup* histories_;
    HistoryGroup::Snapshot snapshot_;
    double value_,dt_;
    int step_;
    bool variable_,pending_=false;
};
inline Step beginStep(TimeStepper& clock,HistoryGroup& histories) { return Step(clock,histories); }
inline Step beginStep(TimeStepper& clock,HistoryGroup& histories,double requestedDt) {
    if (!(requestedDt>0) || !std::isfinite(requestedDt))
        throw std::invalid_argument("requested step requires positive finite dt");
    return Step(clock,histories,requestedDt);
}
}
