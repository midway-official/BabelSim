#pragma once

#include "babelsim/solver.h"
#include "babelsim/methods.h"
#include "babelsim/time.h"

#include <memory>

namespace babelsim {
struct RuntimeControl {
    Methods methods;
    TimeControl time;
    LinearSolverConfig scalar_solver{};
    LinearSolverConfig vector_solver{};

    void validate() const;
};


namespace detail {
class FvmExecution;
FvmExecution& execution();
}

class RunTime {
public:
    // 若进程已经进入 MPI，构造时自动绑定当前通信器；否则为不需要 MPI_Init 的串行运行。
    static RunTime forMesh(const Mesh& mesh, RuntimeControl control = {});

    ~RunTime();
    RunTime(RunTime&&) = delete;
    RunTime& operator=(RunTime&&) = delete;
    RunTime(const RunTime&) = delete;
    RunTime& operator=(const RunTime&) = delete;

    const Mesh& mesh() const;
    const Methods& methods() const;
    bool loop();
    const TimeControl& timeControl() const;
    const LinearSolverConfig& linearControl(bool vector) const;
    // Set evaluation metadata only: no history advancement or output.
    void setTime(double value, int step, double dt);
    double time() const;
    double deltaT() const;
    int step() const;
    bool primary() const;

    // 仅供 math、equ、诊断和内部算法桥接使用。每个线程同时只能有一个活动运行域，
    // 使显式场运算与方程装配共享同一个执行域，同时避免 Field/Mesh 反向依赖 Runtime。
    static RunTime& current();

private:
    // Case records result-writing as an I/O phase; it is not part of the
    // numerical solver counters or convergence decision.
    void recordOutput(double seconds);
    // 性能快照属于 Case/Application 生命周期观测，不是 Physics Solver API。
    friend class Case;
    PerformanceCounters performance() const;
    explicit RunTime(const Mesh& mesh, RuntimeControl control);
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
    friend detail::FvmExecution& detail::execution();
};
}  // babelsim 命名空间
