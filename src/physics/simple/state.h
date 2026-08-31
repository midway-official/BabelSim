#pragma once

#include "babelsim/simple.h"
#include "babelsim/simple_control.h"
#include "internal/simple_discretization.h"

namespace babelsim {

// 只在框架实现中可见。U/p/phi 属于物理场；算法量与预分配工作量分开存放。
// 本对象只在初始化时分配一次，外迭代不创建新的整场工作数组。
struct SimpleSolver::State {
    State(RunTime& run_time, VectorField& U, ScalarField& p, ScalarField& phi,
          FluidProperties fluid, SimpleControl control);

    VectorField& m_U;
    ScalarField& m_p;
    ScalarField& m_phi;
    const Mesh& m_mesh;
    FluidProperties m_fluid;
    SimpleControl m_control;
    Methods m_methods;
    // 算法状态跨一次 SIMPLE 外迭代的多个步骤存在，并直接对应 p'、rAU、phiHbyA。
    // 它与物理状态 U/p/phi 分离，也不包含可随时重算的梯度工作量。
    struct AlgorithmState {
        explicit AlgorithmState(const Mesh& mesh)
            : p_prime(mesh, FieldLocation::Cell, "pPrime"),
              rAU(mesh, FieldLocation::Cell, "rAU"),
              phiHbyA(mesh, FieldLocation::Face, "phiHbyA"),
              previous_velocity(mesh, FieldLocation::Cell, "UPrevious")
        {}

        ScalarField p_prime;
        ScalarField rAU;
        ScalarField phiHbyA;
        VectorField previous_velocity;
    };

    // 数值工作区仅为避免每次外迭代重新分配完整 Field；它没有独立物理生命周期。
    struct NumericalWorkspace {
        explicit NumericalWorkspace(const Mesh& mesh)
            : grad_p(mesh, FieldLocation::Cell, "gradP"),
              rAU_grad_p(mesh, FieldLocation::Cell, "rAUGradP"),
              rAU_grad_p_face(mesh, FieldLocation::Face, "rAUGradPFace"),
              rAU_face(mesh, FieldLocation::Face, "rAUFace"),
              div_phiHbyA(mesh, FieldLocation::Cell, "divPhiHbyA")
        {}

        VectorField grad_p;
        VectorField rAU_grad_p;
        VectorField rAU_grad_p_face;
        ScalarField rAU_face;
        ScalarField div_phiHbyA;
    };

    struct PressureEquationResult {
        SolveResult linear;
        bool healthy = false;
        bool linear_converged = false;
    };

    std::array<SolveResult, 3> solveMomentum();
    PressureEquationResult solvePressure();
    void predictMomentumFlux();
    void correctVelocity();
    void correctFlux();
    void checkContinuityAndConvergence(
        SimpleIterationResult& result,
        const PressureEquationResult& pressure) const;

    AlgorithmState m_algorithm;
    NumericalWorkspace m_workspace;
    bool m_has_fixed_pressure = false;
    SimpleIterationResult m_result;
    PressureEquationResult m_pressure_result;
    enum class Step { Ready, Momentum, Pressure, Velocity, Flux, Complete };
    Step m_step = Step::Ready;
    int m_iteration = 0;
    bool m_log = false;
    void report() const;
    void requireStep(Step expected) const;
};
}  // babelsim 命名空间
