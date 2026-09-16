#pragma once
#include "babelsim/field.h"

namespace babelsim::math {
// A numerical data container, not a time controller. Only saveOld() advances it.
template<class T> class History {
public:
    explicit History(const Field<T>& field): previous_(field), older_(field) {}
    void save(const Field<T>& field,double dt) {
        if(!(dt>0) || !std::isfinite(dt)) throw std::invalid_argument("history requires positive finite dt");
        older_=previous_; previous_=field;
        previous_dt_=dt_; dt_=dt; ++saved_;
    }
    const Field<T>& previous() const {return previous_;}
    const Field<T>& older() const {return older_;}
    double dt() const {return dt_;}
    double previousDt() const {return previous_dt_;}
    int levels() const {return saved_;}
private:
    Field<T> previous_, older_;
    double dt_=0, previous_dt_=0;
    int saved_=0;
};
template<class T> History<T> history(const Field<T>& field) {return History<T>(field);}
template<class T> void saveOld(History<T>& old,const Field<T>& field,double dt) {old.save(field,dt);}
}
