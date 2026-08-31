#pragma once

#include <memory>

namespace babelsim {
class Case;
class RunTime;
struct IncompressibleFields;
struct FluidProperties;
struct SimpleControl;
struct SimpleIterationResult;

// 稳态层流 SIMPLE：调用者只看到算法步骤，不需要管理 pPrime/rAU 或数值工作区。
// Case 或传入的 RunTime/物理场必须活得比算法更久。
class SimpleSolver {
public:
    explicit SimpleSolver(Case& problem);
    SimpleSolver(RunTime& run_time, IncompressibleFields& fields,
                 FluidProperties fluid, SimpleControl control);
    ~SimpleSolver();
    SimpleSolver(const SimpleSolver&) = delete;
    SimpleSolver& operator=(const SimpleSolver&) = delete;

    bool loop() const;
    void solveMomentum();
    void solvePressure();
    void correctVelocity();
    void correctFlux();
    void checkContinuity();
    bool converged() const;

    // 数值测试及嵌入式调用的单次完整外迭代，与上述五步共用同一条执行路径。
    SimpleIterationResult iterate();

private:
    struct State;
    std::unique_ptr<State> m_state;
};

}  // babelsim 命名空间
