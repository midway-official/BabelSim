#pragma once
#include "babelsim/history.h"
#include "babelsim/methods.h"

namespace babelsim::time {
// For capacity*d(x)/dt: current*x_new = previous*x_old + older*x_older
// in the time contribution to Ax=b. All coefficients have units 1/time.
// This describes the existing discretization; assembly still uses equ::ddt.
struct Weights {
    double current=0, previous=0, older=0;
    TimeMethod method=TimeMethod::Steady;
};
inline Weights weights(TimeMethod method, double dt, double previousDt=0, int levels=1) {
    if (method==TimeMethod::Steady) return {};
    if (method!=TimeMethod::Euler && method!=TimeMethod::BDF2)
        throw std::invalid_argument("unsupported time method");
    if (!(dt>0) || !std::isfinite(dt)) throw std::invalid_argument("time weights require positive finite dt");
    if (levels<1) throw std::logic_error("save history before evaluating transient time weights");
    if (method==TimeMethod::Euler || levels<2) {
        if (!std::isfinite(1/dt)) throw std::invalid_argument("time weights exceed representable precision");
        return {1/dt,1/dt,0,TimeMethod::Euler};
    }
    if (!(previousDt>0) || !std::isfinite(previousDt))
        throw std::invalid_argument("BDF2 weights require positive finite previous dt");
    const double ratio=dt/previousDt;
    const Weights result{(1+2*ratio)/((1+ratio)*dt),(1+ratio)/dt,
        -ratio*ratio/((1+ratio)*dt),TimeMethod::BDF2};
    if (!std::isfinite(result.current) || !std::isfinite(result.previous) || !std::isfinite(result.older))
        throw std::invalid_argument("time weights exceed representable precision");
    return result;
}
template<class T> Weights weights(TimeMethod method,const History<T>& history) {
    return weights(method,history.dt(),history.previousDt(),history.levels());
}
template<class A,class B> void requireMatching(const History<A>& a,const History<B>& b) {
    if (a.levels()!=b.levels() || a.dt()!=b.dt() || a.previousDt()!=b.previousDt()
        || &a.previous().mesh()!=&b.previous().mesh())
        throw std::invalid_argument("associated time histories are not synchronized");
}
// Predict only when the algorithm requests it. The target's physical boundaries
// stay intact. Euler/startup retains the saved old state; BDF2 uses linear time
// extrapolation for deferred spatial corrections, including variable steps.
template<class T> void predict(Field<T>& target,const History<T>& history,TimeMethod method) {
    const auto w=weights(method,history);
    if (&target.mesh()!=&history.previous().mesh() || target.location()!=history.previous().location())
        throw std::invalid_argument("time prediction requires matching field/history layouts");
    if (w.method==TimeMethod::Steady) return;
    if (w.method==TimeMethod::Euler) { target.assign(history.previous()); return; }
    const double ratio=history.dt()/history.previousDt();
    target.assignScaled(1+ratio,history.previous());
    target.addScaled(-ratio,history.older());
}
}
