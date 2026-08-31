#pragma once

#include "babelsim/solver.h"
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
struct ScalarEquationControl;

namespace detail {
SolveResult solve(
    const ScalarEquationDefinition& equation,
    ScalarEquationControl control);
std::array<SolveResult, 3> solve(
    const VectorEquationDefinition& equation,
    VectorEquationControl control);
}  // detail 命名空间


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
    SolveResult solve(const ScalarEquationDefinition& equation, ScalarEquationControl control);
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
    void subtract(
        const ScalarField& coefficient,
        fvc::ScalarGradient operation,
        VectorField& target);
    void subtract(fvc::ScalarDiffusionFlux operation, ScalarField& target);
    double relativeChange(const VectorField& current, const VectorField& previous) const;
    double relativeChange(const ScalarField& current, const ScalarField& previous) const;
    double relativeMagnitude(const ScalarField& value, const ScalarField& reference) const;
    FluxBalance fluxBalance(const ScalarField& face_flux) const;
    bool all(bool local_condition) const;

    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;

    friend SolveResult solve(const ScalarEquationDefinition& equation);
    friend SolveResult solve(const VectorEquationDefinition& equation);
    friend SolveResult detail::solve(
        const ScalarEquationDefinition&, ScalarEquationControl);
    friend std::array<SolveResult, 3> detail::solve(
        const VectorEquationDefinition&, VectorEquationControl);
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
    friend void fvc::subtract(
        const ScalarField&, fvc::ScalarGradient, VectorField&);
    friend void fvc::subtract(fvc::ScalarDiffusionFlux, ScalarField&);
    friend double diagnostics::relativeChange(const VectorField&, const VectorField&);
    friend double diagnostics::relativeChange(const ScalarField&, const ScalarField&);
    friend double diagnostics::relativeMagnitude(const ScalarField&, const ScalarField&);
    friend FluxBalance diagnostics::fluxBalance(const ScalarField&);
    friend bool diagnostics::all(bool);
};

}  // babelsim 命名空间
