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
    double time() const;
    double deltaT() const;
    int step() const;
    bool primary() const;

    // 仅供 eqn/math、诊断和内部算法桥接使用。每个线程同时只能有一个活动运行域，
    // 使 solve(equation) 的含义明确，同时避免 Field/Mesh 反向依赖 Runtime。
    static RunTime& current();

private:
    explicit RunTime(const Mesh& mesh, RuntimeControl control);
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
    friend detail::FvmExecution& detail::execution();
};
}  // babelsim 命名空间
