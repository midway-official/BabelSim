#include "internal/compute_backend.h"
#include "internal/fvm_execution.h"

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
        const ScalarDiscreteEquation& equation, ScalarField&) override
    {
        equation.validateStorage();
        ++scalar_solves;
        return {SolveStatus::Converged, 0, 0.0, 0.0, 0.0};
    }

    std::array<SolveResult, 3> solve(
        const VectorDiscreteEquation& equation, VectorField&) override
    {
        equation.validateStorage();
        ++vector_solves;
        const SolveResult result{SolveStatus::Converged, 0, 0.0, 0.0, 0.0};
        return {result, result, result};
    }

    int scalar_solves = 0;
    int vector_solves = 0;
    int synchronizations = 0;
};

}  // 匿名命名空间

int main() {
    const Mesh mesh = Mesh::cartesian({2, 1, 1}, {0, 0, 0}, {2, 1, 1});
    auto backend = std::make_unique<RecordingBackend>();
    RecordingBackend* recording = backend.get();
    Methods methods;
    methods.time = TimeMethod::Euler;
    detail::FvmExecution fvm(mesh, methods, std::move(backend), 0.1);
    fvm.beginStep(0.1);

    ScalarField temperature(mesh, FieldLocation::Cell, "T");
    VectorField velocity(mesh, FieldLocation::Cell, "U");
    VectorField gradient(mesh, FieldLocation::Cell, "gradT");

    require(
        fvm.solve(eqn::ddt(temperature) == eqn::source(1.0), {}).converged(),
        "replaceable backend did not solve a scalar equation");
    require(
        fvm.solve(eqn::ddt(velocity) == Vec3{}, {}).front().converged(),
        "replaceable backend did not solve a vector equation");
    fvm.evaluate(math::grad(temperature), gradient);

    require(recording->scalar_solves == 1, "scalar equation bypassed compute backend");
    require(recording->vector_solves == 1, "vector equation bypassed compute backend");
    require(recording->synchronizations == 2, "explicit operator bypassed backend synchronization");
    require(fvm.all(true), "global logical reduction bypassed compute backend");
    std::cout << "backend_interface_test: replaceable coarse-grained backend passed\n";
}
