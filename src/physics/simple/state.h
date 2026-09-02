#pragma once

#include "babelsim/simple.h"
#include "babelsim/simple_control.h"

namespace babelsim {

// 只在 SIMPLE 模块内可见。用字段分组区分职责，不再为每组增加嵌套类型。
// 所有工作场在构造时分配一次；外迭代只更新数值，不创建整场临时数组。
struct SimpleSolver::State {
    State(VectorField& U, ScalarField& p, ScalarField& phi,
          FluidProperties fluid, SimpleControl control);

    // 借用物理场；Case/调用者拥有它们及网格，必须活得比算法更久。
    VectorField& m_U;
    ScalarField& m_p;
    ScalarField& m_phi;
    const Mesh& m_mesh;
    FluidProperties m_fluid;
    SimpleControl m_control;
    Methods m_methods;

    // 跨步骤保留的算法量。m_mesh 必须先于这些 Field 初始化；场名用于格式选择。
    ScalarField m_p_prime{m_mesh, FieldLocation::Cell, "pPrime"};
    ScalarField m_rAU{m_mesh, FieldLocation::Cell, "rAU"};
    ScalarField m_phiHbyA{m_mesh, FieldLocation::Face, "phiHbyA"};
    VectorField m_previous_velocity{m_mesh, FieldLocation::Cell, "UPrevious"};

    // 可复用的数学中间场；不是底层离散/通信工作区。
    VectorField m_grad_p{m_mesh, FieldLocation::Cell, "gradP"};
    VectorField m_rAU_grad_p{m_mesh, FieldLocation::Cell, "rAUGradP"};
    VectorField m_rAU_grad_p_face{m_mesh, FieldLocation::Face, "rAUGradPFace"};
    ScalarField m_rAU_face{m_mesh, FieldLocation::Face, "rAUFace"};
    ScalarField m_div_phiHbyA{m_mesh, FieldLocation::Cell, "divPhiHbyA"};

    bool m_has_fixed_pressure = false;
    SimpleIterationResult m_result;
    // 最后一次压力结果存在 m_result.pressure；两项累计所有非正交子求解的状态。
    bool m_pressure_healthy = false;
    bool m_pressure_converged = false;
    enum class Step { Ready, Momentum, Pressure, Velocity, Flux, Complete };
    Step m_step = Step::Ready;
    int m_iteration = 0;
    bool m_log = false;
    void predictMomentumFlux();
    void report() const;
    void requireStep(Step expected) const;
};
}  // babelsim 命名空间
