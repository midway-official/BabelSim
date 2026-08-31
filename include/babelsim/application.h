#pragma once

#include <cstddef>

namespace babelsim {
class Case;

// 显式选择表，不使用注册宏、全局 Registry 或动态工厂。外部 Solver 与内置 Solver
// 使用相同入口；初始化、错误传播、并行生命周期和最终输出由实现负责。
struct SolverEntry {
    const char* name;
    int (*run)(Case&);
};

int runApplication(int argc, char* argv[], const SolverEntry* solvers, std::size_t count);
inline int runApplication(int argc, char* argv[], SolverEntry solver) {
    return runApplication(argc, argv, &solver, 1);
}
template <std::size_t N>
int runApplication(int argc, char* argv[], const SolverEntry (&solvers)[N]) {
    return runApplication(argc, argv, solvers, N);
}

}  // babelsim 命名空间
