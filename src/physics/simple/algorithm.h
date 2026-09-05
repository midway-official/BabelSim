#pragma once

#include "../simple_common.h"

#include <memory>

namespace babelsim {
class Case;

// 稳态 SIMPLE 的模块私有算法对象，只供本目录入口和维护测试使用。
class SteadySimpleAlgorithm {
public:
    explicit SteadySimpleAlgorithm(Case& problem);
    SteadySimpleAlgorithm(IncompressibleFields& fields,
                          FluidProperties fluid, SimpleControl control);
    ~SteadySimpleAlgorithm();
    SteadySimpleAlgorithm(const SteadySimpleAlgorithm&) = delete;
    SteadySimpleAlgorithm& operator=(const SteadySimpleAlgorithm&) = delete;

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
