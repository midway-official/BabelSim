#pragma once

#include "babelsim/runtime.h"

#include <array>

namespace babelsim {

// 不可压缩牛顿流体的常物性模型。若未来引入变黏度模型，仍由该层提供等价物性场，
// 而不是让 SIMPLE 接触底层存储或 Case 解析。
struct FluidProperties {
    double density = 1.0;
    double dynamic_viscosity = 1e-3;

    void validate() const;
};

// SIMPLE 只拥有算法控制；空间方法、时间方法和线性后端属于 RunTime 配置。
struct SimpleControl {
    int max_iterations = 1000;
    int non_orthogonal_corrections = 1;
    double velocity_relaxation = 0.7;
    double pressure_relaxation = 0.3;
    double continuity_tolerance = 1e-8;
    double velocity_tolerance = 1e-7;
    LinearSolverConfig velocity_solver{};
    LinearSolverConfig pressure_solver{
        LinearSolverType::ConjugateGradient,
        PreconditionerType::IncompleteCholesky,
        1e-12,
        1e-8,
        1000,
        false,
    };

    void validate() const;
};

inline RuntimeControl simpleRunTimeControl(
    const Methods& methods,
    const SimpleControl& control)
{
    RuntimeControl result;
    result.methods = methods;
    result.scalar_solver = control.pressure_solver;
    result.vector_solver = control.velocity_solver;
    return result;
}

struct IncompressibleFields {
    explicit IncompressibleFields(const Mesh& mesh)
        : velocity(mesh, FieldLocation::Cell, "U"),
          pressure(mesh, FieldLocation::Cell, "p"),
          face_flux(mesh, FieldLocation::Face, "phi")
    {}

    VectorField velocity;
    ScalarField pressure;
    ScalarField face_flux;
};

// 通用 RunTime 计算的面通量平衡在不可压缩流中就是连续性度量。保留熟悉的名称，
// 但不把全局归约和数据分区暴露给 SIMPLE。
using ContinuityMetrics = FluxBalance;

struct SimpleIterationResult {
    std::array<SolveResult, 3> velocity;
    SolveResult pressure;
    ContinuityMetrics continuity;
    double relative_velocity_change = 0.0;
    double relative_pressure_correction = 0.0;
    bool healthy = false;
    bool linear_converged = false;
    bool converged = false;
};

// 稳态层流不可压缩 SIMPLE 算法。它是 Algorithm-driven Solver：内部组织多个
// Equation-driven 方程与修正步骤，RunTime 隐藏 FVM 装配、线性代数和并行后端。
class SimpleSolver {
public:
    SimpleSolver(
        RunTime& run_time,
        IncompressibleFields& fields,
        FluidProperties fluid,
        SimpleControl control);
    ~SimpleSolver() = default;
    SimpleSolver(const SimpleSolver&) = delete;
    SimpleSolver& operator=(const SimpleSolver&) = delete;

    // 单次调用正好对应 SIMPLE 的一个外迭代：动量、Rhie-Chow、压力修正、速度/通量
    // 修正和连续性检查。外层停止判据是 RunTime 内部全局归约后的值。
    SimpleIterationResult iterate();

private:
    IncompressibleFields& m_fields;
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

        void predictMomentumFlux(
            const Methods& methods,
            const VectorField& U,
            const ScalarField& p,
            const ScalarField& rAU,
            ScalarField& phiHbyA);
        void preparePressureEquation(const ScalarField& phiHbyA);

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

    void savePreviousState();
    std::array<SolveResult, 3> solveMomentum();
    PressureEquationResult solvePressure();
    void predictMomentumFlux();
    void correctVelocity();
    void correctFlux();
    void checkContinuityAndConvergence(
        SimpleIterationResult& result,
        const PressureEquationResult& pressure) const;
    void initializePressureCorrectionBoundaries();

    AlgorithmState m_algorithm;
    NumericalWorkspace m_workspace;
    bool m_has_fixed_pressure = false;
};

}  // babelsim 命名空间
