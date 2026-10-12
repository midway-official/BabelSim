#pragma once
#include "babelsim/field.h"
#include "babelsim/time.h"
#include <functional>
#include <memory>
#include <optional>

namespace babelsim::time {
// Explicit time levels. Only save() advances history, once per physical timestep.
template<class T> class History {
public:
    explicit History(const Field<T>& field): previous_(field), older_(field) {}
    void save(const Field<T>& field,double dt) {
        if(!(dt>0) || !std::isfinite(dt)) throw std::invalid_argument("history requires positive finite dt");
        if (&field.mesh() != &previous_.mesh() || field.location() != previous_.location())
            throw std::invalid_argument("history requires its original mesh and location");
        older_=previous_; previous_=field;
        previous_dt_=dt_; dt_=dt; ++saved_;
    }
    const Field<T>& previous() const {return previous_;}
    const Field<T>& older() const {return older_;}
    double dt() const {return dt_;}
    double previousDt() const {return previous_dt_;}
    int levels() const {return saved_;}
    void restoreSnapshot(const History& snapshot) {
        previous_.restoreSnapshot(snapshot.previous_);
        older_.restoreSnapshot(snapshot.older_);
        dt_=snapshot.dt_; previous_dt_=snapshot.previous_dt_; saved_=snapshot.saved_;
    }
private:
    Field<T> previous_, older_;
    double dt_=0, previous_dt_=0;
    int saved_=0;
    friend class HistoryGroup;
};
template<class T> History<T> history(const Field<T>& field) {return History<T>(field);}

// Heterogeneous history group: fields are borrowed, addHistory histories are
// owned; add(field, history) borrows an external history. Algorithms choose when to save.
class HistoryGroup {
    struct Metadata {
        int levels; double dt, previousDt;
        bool operator==(const Metadata& other) const {
            return levels==other.levels && dt==other.dt && previousDt==other.previousDt;
        }
    };
    struct Entry {
        void* field;
        void* history;
        const int* levels;
        const double* dt;
        const double* previousDt;
        Metadata metadata() const { return {*levels,*dt,*previousDt}; }
        void (*save)(Entry&,double);
        std::function<std::function<void()>()> checkpoint;
    };
public:
    HistoryGroup() = default;
    HistoryGroup(const HistoryGroup&) = delete;
    HistoryGroup& operator=(const HistoryGroup&) = delete;
    // The group owns these histories at stable addresses. Construction captures
    // the final initialized values; rollback snapshots remain explicitly opt-in.
    template<class T> History<T>& addHistory(Field<T>& field) {
        requireDeclaration();
        auto history = std::make_shared<History<T>>(field);
        owned_.push_back(history);
        try { add(field, *history); }
        catch (...) { owned_.pop_back(); throw; }
        return *history;
    }
    template<class T> void add(Field<T>& field, History<T>& history) {
        requireDeclaration();
        field.validateStorage();
        if (&field.mesh()!=&history.previous().mesh() || field.location()!=history.previous().location())
            throw std::invalid_argument("history group requires matching field/history layouts");
        requireMesh(field.mesh());
        for (const auto& entry:entries_)
            if (entry.field==&field || entry.history==&history)
                throw std::invalid_argument("duplicate history group field or history");
        for (const void* fieldId:tracked_)
            if (fieldId==&field) throw std::invalid_argument("field is already tracked; add its history first");
        const Metadata metadata{history.levels(),history.dt(),history.previousDt()};
        if (expected_ && !(metadata==*expected_))
            throw std::invalid_argument("cannot group histories at different time levels");
        Entry entry;
        entry.field=&field; entry.history=&history;
        // Metadata is read directly in the per-step checks. Type erasure needs
        // only one callback for the actual heterogeneous history save.
        entry.levels=&history.saved_; entry.dt=&history.dt_; entry.previousDt=&history.previous_dt_;
        entry.save=[](Entry& item,double dt) {
            static_cast<History<T>*>(item.history)->save(*static_cast<Field<T>*>(item.field),dt);
        };
        entry.checkpoint=[&field,&history] {
            auto fieldCopy=std::make_shared<Field<T>>(field);
            auto historyCopy=std::make_shared<History<T>>(history);
            return std::function<void()>([&field,&history,fieldCopy,historyCopy] {
                field.restoreSnapshot(*fieldCopy); history.restoreSnapshot(*historyCopy);
            });
        };
        entries_.push_back(std::move(entry)); expected_=metadata;
    }
    // Include pressure, derived coefficients or other fields in optional rollback
    // without making them a time-discretized unknown.
    template<class T> void trackForRollback(Field<T>& field) {
        requireDeclaration(); requireMesh(field.mesh()); field.validateStorage();
        for (const auto& entry:entries_)
            if (entry.field==&field) throw std::invalid_argument("history field is already tracked");
        for (const void* fieldId:tracked_)
            if (fieldId==&field) throw std::invalid_argument("duplicate tracked field");
        tracked_.push_back(&field);
        checkpoints_.push_back([&field] {
            auto copy=std::make_shared<Field<T>>(field);
            return std::function<void()>([&field,copy] { field.restoreSnapshot(*copy); });
        });
    }
    template<class T, class U, class... Rest>
    void trackForRollback(Field<T>& first, Field<U>& second, Rest&... rest) {
        trackForRollback(first);
        trackForRollback(second);
        (trackForRollback(rest), ...);
    }
    // Copyable model state (e.g. counters or a closure state struct). Registration
    // is explicit; opaque pointers or external side effects cannot be rolled back.
    template<class T> void watch(T& value) {
        requireDeclaration();
        for (const void* id:watched_)
            if (id==&value) throw std::invalid_argument("duplicate watched model state");
        watched_.push_back(&value);
        checkpoints_.push_back([&value] {
            auto copy=std::make_shared<T>(value);
            return std::function<void()>([&value,copy] { value=*copy; });
        });
    }
    template<class T> void watch(Field<T>&) = delete; // use trackForRollback(field)
    template<class T> void watch(History<T>&) = delete; // use add(field, history)
    void validate() const {
        for (const auto& entry:entries_)
            if (!expected_ || !(entry.metadata()==*expected_))
                throw std::logic_error("history group is out of sync: save all members through the group");
    }
    // Expert path when no TimeStepper exists; save(clock) also guards duplicate
    // saves within one physical step. Neither belongs in a nonlinear inner loop.
    void save(double dt) {
        if (clock_) throw std::logic_error("use save(clock) for a clock-bound history group");
        saveImpl(dt);
    }
private:
    void saveImpl(double dt) {
        if (active_) throw std::logic_error("cannot save histories during a pending step");
        validate();
        if (!(dt>0) || !std::isfinite(dt)) throw std::invalid_argument("history group requires positive finite dt");
        if (entries_.empty()) throw std::logic_error("history group has no histories");
        sealed_=true;
        for (auto& entry:entries_) entry.save(entry,dt);
        expected_=entries_.front().metadata();
    }
public:
    void save(const TimeStepper& clock) {
        if (clock.step()<=0) throw std::logic_error("advance the clock before saving histories");
        if (mesh_!=&clock.mesh())
            throw std::invalid_argument("history group and clock require the same mesh");
        if (clock_ && clock_!=&clock) throw std::invalid_argument("history group is bound to another clock");
        if (clock.step()<=lastStep_) throw std::logic_error("histories already saved for this physical step");
        if (clock.step()!=lastStep_+1) throw std::logic_error("history group skipped a physical step");
        saveImpl(clock.dt()); clock_=&clock; lastStep_=clock.step();
    }
    struct Snapshot {
    private:
        const HistoryGroup* owner=nullptr;
        std::vector<std::function<void()>> restore;
        std::optional<Metadata> expected;
        const TimeStepper* clock=nullptr;
        int lastStep=0;
        bool sealed=false;
        std::size_t registrations=0;
        friend class HistoryGroup;
    };
    Snapshot checkpoint() const {
        validate(); Snapshot snapshot; snapshot.owner=this;
        snapshot.expected=expected_; snapshot.clock=clock_; snapshot.lastStep=lastStep_; snapshot.sealed=sealed_;
        snapshot.registrations=entries_.size()+checkpoints_.size();
        for (const auto& entry:entries_) snapshot.restore.push_back(entry.checkpoint());
        for (const auto& make:checkpoints_) snapshot.restore.push_back(make());
        return snapshot;
    }
    void restore(const Snapshot& snapshot) {
        if (snapshot.owner!=this) throw std::invalid_argument("history snapshot belongs to another group");
        if (snapshot.registrations!=entries_.size()+checkpoints_.size())
            throw std::logic_error("history group registrations changed after checkpoint");
        for (const auto& restore:snapshot.restore) restore();
        expected_=snapshot.expected; clock_=snapshot.clock; lastStep_=snapshot.lastStep; sealed_=snapshot.sealed;
    }
private:
    void requireDeclaration() const {
        if (sealed_ || active_) throw std::logic_error("register all histories and state before time stepping");
    }
    void requireMesh(const Mesh& mesh) {
        if (mesh_ && mesh_!=&mesh) throw std::invalid_argument("history group fields require the same mesh");
        mesh_=&mesh;
    }
    const Mesh* mesh_=nullptr;
    std::vector<std::shared_ptr<void>> owned_;
    std::vector<Entry> entries_;
    std::vector<const void*> tracked_, watched_;
    std::vector<std::function<std::function<void()>()>> checkpoints_;
    std::optional<Metadata> expected_;
    const TimeStepper* clock_=nullptr;
    int lastStep_=0;
    bool sealed_=false, active_=false;
    friend class Step;
};
}
