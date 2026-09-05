#pragma once

#include "../simple_common.h"

#include <memory>

namespace babelsim {
class Case;

// 瞬态 SIMPLE 的模块私有算法对象。它服务本目录的运行入口与维护测试，
// 不是普通 Solver 作者需要依赖的 BabelSim 公共 API。
class TransientSimpleAlgorithm {
public:
    explicit TransientSimpleAlgorithm(Case& problem);
    TransientSimpleAlgorithm(IncompressibleFields& fields,
                             FluidProperties fluid, SimpleControl control);
    ~TransientSimpleAlgorithm();
    TransientSimpleAlgorithm(const TransientSimpleAlgorithm&) = delete;
    TransientSimpleAlgorithm& operator=(const TransientSimpleAlgorithm&) = delete;

    // Case::loop() 进入新物理时间步后调用一次；时间历史仍由框架维护。
    void beginTimeStep();
    bool loop() const;
    void solveMomentum();
    void solvePressure();
    void correctVelocity();
    void correctFlux();
    void correctTurbulence();
    void checkContinuity();
    bool converged() const;
    SimpleIterationResult iterate();

private:
    struct State;
    std::unique_ptr<State> m_state;
};

}  // babelsim 命名空间
