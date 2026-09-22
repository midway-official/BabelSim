#include "internal/compute_backend.h"

#include "babelsim/equ.h"
#include "babelsim/parallel.h"
#include "babelsim/runtime.h"

#include "test_util.h"

#include <algorithm>
#include <iostream>
#include <memory>

using namespace babelsim;

namespace {

class RecordingBackend final : public detail::ComputeBackend {
public:
    void synchronize(ScalarField&) override { ++synchronizations; }
    void synchronize(VectorField&) override { ++synchronizations; }
    void synchronize(TensorField&) override { ++synchronizations; }

    void sum(const double* local, double* global, int count) const override {
        std::copy_n(local, count, global);
    }

    void maximum(const double* local, double* global, int count) const override {
        std::copy_n(local, count, global);
    }

    bool all(bool local_condition) const override { return local_condition; }

    SolveResult solve(
        const ScalarDiscreteEquation& equation, ScalarField&)
    {
        equation.validateStorage();
        ++scalar_solves;
        return {SolveStatus::Converged, 0, 0.0, 0.0, 0.0};
    }

    std::array<SolveResult, 3> solve(
        const VectorDiscreteEquation& equation, VectorField&)
    {
        equation.validateStorage();
        ++vector_solves;
        const SolveResult result{SolveStatus::Converged, 0, 0.0, 0.0, 0.0};
        return {result, result, result};
    }

    SolveResult solve(const ScalarDiscreteEquation& equation, ScalarField& field,
                      const LinearSolverConfig& config) override {
        selected = config; ++configured_solves; return solve(equation, field);
    }
    std::array<SolveResult,3> solve(const VectorDiscreteEquation& equation, VectorField& field,
                                  const LinearSolverConfig& config) override {
        selected = config; ++configured_solves; return solve(equation, field);
    }
    LinearSolverConfig selected;
    int configured_solves = 0;
    int scalar_solves = 0;
    int vector_solves = 0;
    int synchronizations = 0;
};

RecordingBackend* recording = nullptr;

}  // 匿名命名空间

// 计算后端只有一个替换点：实现 makeComputeBackend() 即可整体换掉默认的 Eigen/MPI
// 装配与求解实现。本测试安装记录型后端，验证 equ:: 与 math:: 只通过 ComputeBackend
// 契约访问后端，而不依赖任何具体矩阵或通信实现。
namespace babelsim::detail {
std::unique_ptr<ComputeBackend> makeComputeBackend(
    const Mesh&, ParallelContext)
{
    auto backend = std::make_unique<RecordingBackend>();
    recording = backend.get();
    return backend;
}
}  // babelsim::detail 命名空间

int main() {
    const Mesh mesh = makeHexBox({2, 1, 1}, {0, 0, 0}, {2, 1, 1});
    RuntimeControl control;
    control.methods.time = TimeMethod::Euler;
    RunTime run_time = RunTime::forMesh(mesh, control);
    require(recording != nullptr, "the run did not install the replaceable backend");

    ScalarField temperature(mesh, FieldLocation::Cell, "T");
    ScalarField old_temperature(mesh, FieldLocation::Cell, "Told");
    auto scalar = testEquation(temperature);
    equ::ddt(scalar, 1.0, old_temperature, 0.1);
    equ::source(scalar, 1.0);
    require(
        equ::solve(scalar).converged(),
        "replaceable backend did not solve a scalar equation");

    VectorField velocity(mesh, FieldLocation::Cell, "U");
    VectorField old_velocity(mesh, FieldLocation::Cell, "Uold");
    auto vector = testEquation(velocity);
    equ::ddt(vector, 1.0, old_velocity, 0.1);
    require(
        equ::solve(vector).converged(),
        "replaceable backend did not solve a vector equation");

    require(recording->scalar_solves == 1, "scalar equation bypassed compute backend");
    require(recording->vector_solves == 1, "vector equation bypassed compute backend");

    const int synchronizations = recording->synchronizations;
    VectorField gradient(mesh, FieldLocation::Cell, "gradT");
    const auto gradientOptions = testEquationControl("gradient").spatial;
    math::evaluate(math::ScalarGradient{temperature, gradientOptions}, gradient);
    require(recording->synchronizations == synchronizations + 2,
            "explicit operator bypassed backend synchronization");
    require(diagnostics::all(true), "global logical reduction bypassed compute backend");
    const int configured_before = recording->configured_solves;
    EquationControl first;
    first.name = "temperature";
    first.spatial.interpolation = InterpolationMethod::Corrected;
    first.spatial.gradient = GradientMethod::LeastSquares;
    first.spatial.convection = ConvectionMethod::Upwind;
    first.spatial.diffusion = DiffusionMethod::Corrected;
    first.linear.max_iterations = 17;
    auto named = equ::createEquation(temperature, first);
    equ::reaction(named, 1.0);
    equ::solve(named);
    require(recording->configured_solves == configured_before + 1 && recording->selected.max_iterations == 17,
            "default solve ignored bound equation configuration");
    auto copied = named.copy(); named.reset();
    equ::solve(copied);
    require(recording->selected.max_iterations == 17, "copy lost linear configuration");
    auto explicitConfig = first.linear; explicitConfig.max_iterations = 29;
    equ::solve(copied, explicitConfig);
    require(recording->selected.max_iterations == 29, "explicit solve did not override binding");
    equ::solve(copied);
    require(recording->selected.max_iterations == 17, "explicit solve mutated bound configuration");
    first.name = "momentum"; first.linear.max_iterations = 41;
    auto namedVector = equ::createEquation(velocity, first);
    equ::reaction(namedVector, 1.0); equ::solve(namedVector);
    require(recording->selected.max_iterations == 41, "vector equation ignored its configuration");
    std::cout << "backend_interface_test: replaceable coarse-grained backend passed\n";
}
