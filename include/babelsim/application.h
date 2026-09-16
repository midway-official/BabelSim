#pragma once

#include "babelsim/solver_control.h"

namespace babelsim {
class Case;
using ApplicationErrorHandler = void (*)(const char* message);

// Physics reports a numerical outcome; only the application maps it to an exit
// code. The same statuses are used for linear solves and whole solver programs.
struct SolverResult {
    SolveStatus status = SolveStatus::Converged;
    static SolverResult completed() { return {SolveStatus::Converged}; }
    static SolverResult notConverged() { return {SolveStatus::MaxIterations}; }
    static SolverResult numericalFailure() { return {SolveStatus::NumericalFailure}; }
};

// 每个 Solver 在自己的源文件命名空间作用域注册一次，名称使用字符串字面量。
// 注册只连接描述项：无堆分配、异常或 MPI 调用，校验与分派在 runApplication 内完成。
// 对象及名称必须在应用运行期间保持有效；不支持运行中的并发注册或动态插件卸载。
class SolverRegistration {
public:
    SolverRegistration(const char* name, SolverResult (*run)(Case&)) noexcept;
    ~SolverRegistration() noexcept;
    SolverRegistration(const SolverRegistration&) = delete;
    SolverRegistration& operator=(const SolverRegistration&) = delete;

private:
    const char* m_name;
    SolverResult (*m_run)(Case&);
    mutable const SolverRegistration* m_next;
    static const SolverRegistration*& first() noexcept;
    friend int runApplication(int argc, char* argv[], ApplicationErrorHandler);
};

// Dispatch and MPI lifetime only. No output or automatic result writes.
// The application may supply its own error handler; default operation is silent.
int runApplication(int argc, char* argv[], ApplicationErrorHandler onError = nullptr);

}  // babelsim 命名空间
