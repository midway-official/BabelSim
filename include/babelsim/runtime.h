#pragma once

#include "babelsim/fvc.h"
#include "babelsim/fvm.h"
#include "babelsim/methods.h"
#include "babelsim/solver_control.h"
#include "babelsim/time.h"

#include <array>
#include <memory>

namespace babelsim {

class RunTime;

// Solver 只提供方法、时间和线性收敛控制；MPI 通信器、halo、LDU 和后端库都由
// RunTime 的私有实现管理。不同 Field 可在一次运行中复用同一个 RunTime。
struct RuntimeControl {
    Methods methods;
    TimeControl time;
    LinearSolverConfig scalar_solver{};
    LinearSolverConfig vector_solver{};

    void validate() const;
};

struct VectorEquationControl;
class SimpleSolver;

struct ScalarEquationControl;

// 面通量在每个控制体上的守恒误差。它是通用有限体积诊断量；不可压缩 SIMPLE 将
// 它命名为连续性残差。归约由 RunTime 完成，调用者不会接触 rank 或通信器。
struct FluxBalance {
    double l1 = 0.0;
    double l2 = 0.0;
    double maximum = 0.0;
    double relative = 0.0;
};

// Runtime 之外的 Solver API：显式量属于 fvc，收敛与守恒量属于 diagnostics，
// 隐式方程由 solve() 处理。它们自动使用当前线程唯一活动的 RunTime，因此 Solver
// 不需要在每个数学操作中传递执行对象。
SolveResult solve(const ScalarEquationDefinition& equation);
std::array<SolveResult, 3> solve(const VectorEquationDefinition& equation);

namespace fvc {
void evaluate(ScalarGradient operation, VectorField& result);
void evaluate(VectorGradient operation, TensorField& result);
void evaluate(FaceFlux operation, ScalarField& result);
void evaluate(FaceDivergence operation, ScalarField& result);
void evaluate(VectorDivergence operation, ScalarField& result);
void evaluate(ScalarConvection operation, ScalarField& result);
void evaluate(VectorConvection operation, VectorField& result);
void evaluate(ScalarInterpolation operation, ScalarField& result);
void evaluate(VectorInterpolation operation, VectorField& result);
void evaluate(ScalarReconstruction operation, ScalarField& result);
void evaluate(VectorReconstruction operation, VectorField& result);
void evaluate(ScalarLaplacian operation, ScalarField& result);
}  // fvc 命名空间

namespace diagnostics {
double relativeChange(const VectorField& current, const VectorField& previous);
double relativeChange(const ScalarField& current, const ScalarField& previous);
double relativeMagnitude(const ScalarField& value, const ScalarField& reference);
FluxBalance fluxBalance(const ScalarField& face_flux);
bool all(bool local_condition);
}  // diagnostics 命名空间

class RunTime {
public:
    // 若进程已经进入 MPI，构造时自动绑定当前通信器；否则为不需要 MPI_Init 的串行运行。
    static RunTime forMesh(const Mesh& mesh, RuntimeControl control = {});

    ~RunTime();
    RunTime(RunTime&&) = delete;
    RunTime& operator=(RunTime&&) = delete;
    RunTime(const RunTime&) = delete;
    RunTime& operator=(const RunTime&) = delete;

    const Mesh& mesh() const;
    const Methods& methods() const;
    bool loop();
    double time() const;
    double deltaT() const;
    int step() const;
    bool primary() const;

    // 仅供 fvm/fvc、诊断和内部算法桥接使用。每个线程同时只能有一个活动运行域，
    // 使 solve(equation) 的含义明确，同时避免 Field/Mesh 反向依赖 Runtime。
    static RunTime& current();

private:
    explicit RunTime(const Mesh& mesh, RuntimeControl control);
    void synchronize(ScalarField& field);
    void synchronize(VectorField& field);
    void synchronize(TensorField& field);
    SolveResult solve(const ScalarEquationDefinition& equation);
    SolveResult solve(const ScalarEquationDefinition& equation, ScalarEquationControl control);
    std::array<SolveResult, 3> solve(const VectorEquationDefinition& equation);
    std::array<SolveResult, 3> solve(
        const VectorEquationDefinition& equation,
        VectorEquationControl control);
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
    double relativeChange(const VectorField& current, const VectorField& previous) const;
    double relativeChange(const ScalarField& current, const ScalarField& previous) const;
    double relativeMagnitude(const ScalarField& value, const ScalarField& reference) const;
    FluxBalance fluxBalance(const ScalarField& face_flux) const;
    bool all(bool local_condition) const;

    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;

    friend SolveResult solve(const ScalarEquationDefinition& equation);
    friend std::array<SolveResult, 3> solve(const VectorEquationDefinition& equation);
    friend class SimpleSolver;
    friend void fvc::evaluate(fvc::ScalarGradient, VectorField&);
    friend void fvc::evaluate(fvc::VectorGradient, TensorField&);
    friend void fvc::evaluate(fvc::FaceFlux, ScalarField&);
    friend void fvc::evaluate(fvc::FaceDivergence, ScalarField&);
    friend void fvc::evaluate(fvc::VectorDivergence, ScalarField&);
    friend void fvc::evaluate(fvc::ScalarConvection, ScalarField&);
    friend void fvc::evaluate(fvc::VectorConvection, VectorField&);
    friend void fvc::evaluate(fvc::ScalarInterpolation, ScalarField&);
    friend void fvc::evaluate(fvc::VectorInterpolation, VectorField&);
    friend void fvc::evaluate(fvc::ScalarReconstruction, ScalarField&);
    friend void fvc::evaluate(fvc::VectorReconstruction, VectorField&);
    friend void fvc::evaluate(fvc::ScalarLaplacian, ScalarField&);
    friend double diagnostics::relativeChange(const VectorField&, const VectorField&);
    friend double diagnostics::relativeChange(const ScalarField&, const ScalarField&);
    friend double diagnostics::relativeMagnitude(const ScalarField&, const ScalarField&);
    friend FluxBalance diagnostics::fluxBalance(const ScalarField&);
    friend bool diagnostics::all(bool);
};

}  // babelsim 命名空间
