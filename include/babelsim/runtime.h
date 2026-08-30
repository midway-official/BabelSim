#pragma once

#include "babelsim/fvc.h"
#include "babelsim/fvm.h"
#include "babelsim/methods.h"
#include "babelsim/solver_control.h"
#include "babelsim/time.h"

#include <array>
#include <memory>

namespace babelsim {

struct MomentumInterpolation;
struct PressureCorrection;

// Solver 只提供方法、时间和线性收敛控制；MPI 通信器、halo、LDU 和后端库都由
// RunTime 的私有实现管理。不同 Field 可在一次运行中复用同一个 RunTime。
struct RuntimeControl {
    Methods methods;
    TimeControl time;
    LinearSolverConfig scalar_solver{};
    LinearSolverConfig vector_solver{};

    void validate() const;
};

// 方程级控制仍保持数学含义：欠松弛和由离散主对角导出的单元 mobility。不会暴露
// 稀疏矩阵、行号或存储格式；SIMPLE 使用 mobility 形成 rAU。
struct VectorEquationControl {
    double relaxation = 1.0;
    ScalarField* mobility = nullptr;
};

struct ScalarEquationControl {
    bool fix_reference = false;
};

// 面通量在每个控制体上的守恒误差。它是通用有限体积诊断量；不可压缩 SIMPLE 将
// 它命名为连续性残差。归约由 RunTime 完成，调用者不会接触 rank 或通信器。
struct FluxBalance {
    double l1 = 0.0;
    double l2 = 0.0;
    double maximum = 0.0;
    double relative = 0.0;
};

class RunTime {
public:
    // 若进程已经进入 MPI，构造时自动绑定当前通信器；否则为不需要 MPI_Init 的串行运行。
    static RunTime forMesh(const Mesh& mesh, RuntimeControl control = {});

    ~RunTime();
    RunTime(RunTime&&) noexcept;
    RunTime& operator=(RunTime&&) noexcept;
    RunTime(const RunTime&) = delete;
    RunTime& operator=(const RunTime&) = delete;

    const Mesh& mesh() const;
    const Methods& methods() const;
    bool loop();
    double time() const;
    double deltaT() const;
    int step() const;
    bool primary() const;

    // 高层场与收敛工具：它们只表达“复制”“相对变化”“守恒”和“全部成立”，
    // 不暴露 owned/ghost、全局归约或 MPI 生命周期。
    void copy(const VectorField& source, VectorField& destination);
    void scale(double factor, const ScalarField& source, ScalarField& destination);
    double relativeChange(const VectorField& current, const VectorField& previous) const;
    double relativeChange(const ScalarField& current, const ScalarField& previous) const;
    FluxBalance fluxBalance(const ScalarField& face_flux) const;
    bool all(bool local_condition) const;

    // 唯一的隐式方程求解入口。Equation 是轻量数学描述；这里才会同步输入、
    // 组装离散方程、选择串行/分布式线性后端并写回解场。
    SolveResult solve(const ScalarEquationDefinition& equation);
    SolveResult solve(
        const ScalarEquationDefinition& equation,
        ScalarEquationControl control);
    std::array<SolveResult, 3> solve(const VectorEquationDefinition& equation);
    std::array<SolveResult, 3> solve(
        const VectorEquationDefinition& equation,
        VectorEquationControl control);

    // fvc 显式求值入口。调用者提供结果场，避免隐式大 Field 临时对象和重复分配。
    void evaluate(fvc::ScalarGradient operation, VectorField& result);
    void evaluate(fvc::VectorGradient operation, TensorField& result);
    void evaluate(fvc::FaceFlux operation, ScalarField& result);
    void evaluate(fvc::FaceDivergence operation, ScalarField& result);
    void evaluate(fvc::VectorDivergence operation, ScalarField& result);
    void evaluate(fvc::ScalarConvection operation, ScalarField& result);
    void evaluate(fvc::VectorConvection operation, VectorField& result);
    void evaluate(fvc::ScalarInterpolation operation, ScalarField& result);
    void evaluate(fvc::VectorInterpolation operation, VectorField& result);
    void evaluate(fvc::ScalarReconstruction operation, ScalarField& result);
    void evaluate(fvc::VectorReconstruction operation, VectorField& result);
    void evaluate(fvc::ScalarLaplacian operation, ScalarField& result);

private:
    explicit RunTime(const Mesh& mesh, RuntimeControl control);
    void synchronize(ScalarField& field);
    void synchronize(VectorField& field);
    void synchronize(TensorField& field);

    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;

    friend struct MomentumInterpolation;
    friend struct PressureCorrection;
};

inline SolveResult solve(RunTime& run_time, const ScalarEquationDefinition& equation) {
    return run_time.solve(equation);
}

inline std::array<SolveResult, 3> solve(
    RunTime& run_time, const VectorEquationDefinition& equation)
{
    return run_time.solve(equation);
}

inline std::array<SolveResult, 3> solve(
    RunTime& run_time,
    const VectorEquationDefinition& equation,
    VectorEquationControl control)
{
    return run_time.solve(equation, control);
}

}  // babelsim 命名空间
