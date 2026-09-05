#pragma once

#include "algorithm.h"

namespace babelsim {

// 瞬态 SIMPLE 的物理场引用、跨校正算法量和可复用数学中间场。
// 所有临时场都属于本算法模块，不向 Solver 调用者或框架底层泄漏。
struct TransientSimpleAlgorithm::State {
    State(VectorField& U, ScalarField& p, ScalarField& phi,
          FluidProperties fluid, SimpleControl control, Case* problem = nullptr);

    // Case 或调用者拥有物理场和网格，必须活得比算法对象更久。
    VectorField& m_U;
    ScalarField& m_p;
    ScalarField& m_phi;
    const Mesh& m_mesh;
    FluidProperties m_fluid;
    SimpleControl m_control;
    Methods m_methods;
    ScalarField m_effective_viscosity{m_mesh, FieldLocation::Cell, "muEffective"};
    std::unique_ptr<rans::Model, void (*)(rans::Model*)> m_turbulence{nullptr, rans::destroy};

    ScalarField m_p_prime{m_mesh, FieldLocation::Cell, "pPrime"};
    ScalarField m_rAU{m_mesh, FieldLocation::Cell, "rAU"};
    ScalarField m_phiHbyA{m_mesh, FieldLocation::Face, "phiHbyA"};
    VectorField m_previous_velocity{m_mesh, FieldLocation::Cell, "UPrevious"};

    VectorField m_grad_p{m_mesh, FieldLocation::Cell, "gradP"};
    VectorField m_rAU_grad_p{m_mesh, FieldLocation::Cell, "rAUGradP"};
    VectorField m_rAU_grad_p_face{m_mesh, FieldLocation::Face, "rAUGradPFace"};
    ScalarField m_rAU_face{m_mesh, FieldLocation::Face, "rAUFace"};
    ScalarField m_div_phiHbyA{m_mesh, FieldLocation::Cell, "divPhiHbyA"};

    bool m_has_fixed_pressure = false;
    SimpleIterationResult m_result;
    bool m_pressure_healthy = false;
    bool m_pressure_converged = false;
    enum class Step { Ready, Momentum, Pressure, Velocity, Flux, Turbulence, Complete };
    Step m_step = Step::Ready;
    // -1 表示尚未由物理时间循环开启时间步。
    int m_iteration = -1;
    bool m_log = false;

    void predictMomentumFlux();
    void report() const;
    void requireStep(Step expected) const;
};

}  // babelsim 命名空间
