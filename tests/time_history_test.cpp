#include "babelsim/runtime.h"
#include "test_util.h"

#include <iostream>

using namespace babelsim;

void checkHistory(TimeMethod method) {
    const Mesh mesh = Mesh::cartesian({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
    ScalarField T(mesh, FieldLocation::Cell, "T", 0.0);
    VectorField U(mesh, FieldLocation::Cell, "U", {});
    RuntimeControl control;
    control.methods.time = method;
    control.time = {0.0, 0.3, 0.1};
    RunTime time = RunTime::forMesh(mesh, control);
    while (time.loop()) {
        // 内迭代不是时间推进：重复求解同一方程，结果仍在同一个物理时间层。
        for (int correction = 0; correction < 3; ++correction) {
            require(solve(fvm::ddt(T) == fvm::source(2.0)).converged(), "scalar solve failed");
            require(solve(fvm::ddt(U) == fvm::source(Vec3{1, 2, 3})).converged(),
                    "vector solve failed");
            require(near(T[0], 2.0 * time.time(), 1e-12), "inner solve advanced scalar history");
            require(near(U[0], Vec3{time.time(), 2*time.time(), 3*time.time()}, 1e-12),
                    "inner solve advanced vector history");
        }
    }
    require(time.step() == 3 && !time.loop(), "time loop ran extra steps");
}

int main() {
    checkHistory(TimeMethod::Euler);
    checkHistory(TimeMethod::BDF2);
    {
        const Mesh mesh = Mesh::cartesian({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
        ScalarField T(mesh, FieldLocation::Cell, "T", 0.0);
        RuntimeControl control;
        control.methods.time = TimeMethod::Euler;
        control.time = {0.0, 0.25, 0.1};
        RunTime time = RunTime::forMesh(mesh, control);
        while (time.loop()) {
            require(solve(fvm::ddt(T) == fvm::source(1.0)).converged(), "short-step solve failed");
            require(near(T[0], time.time(), 1e-12), "short final step used wrong deltaT");
        }
        require(time.step() == 3 && near(time.time(), 0.25), "endTime was not respected");
        require(near(time.deltaT(), 0.05), "incorrect last deltaT");
    }
    bool rejected = false;
    try {
        RuntimeControl control;
        control.methods.time = TimeMethod::BDF2;
        control.time = {0, 0.25, 0.1};
        control.validate();
    } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "variable-step BDF2 must not silently use uniform-step coefficients");
    {
        const Mesh mesh = Mesh::cartesian({6, 5, 1}, {0, 0, 0}, {1, 1, 1});
        VectorField U(mesh, FieldLocation::Cell, "U");
        for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch)
            U.boundary(patch) = fixedValue(Vec3{});
        RuntimeControl control;
        control.vector_solver.max_iterations = 1;
        control.vector_solver.absolute_tolerance = 1e-30;
        control.vector_solver.relative_tolerance = 1e-25;
        RunTime time = RunTime::forMesh(mesh, control);
        const SolveResult result = solve(
            -fvm::laplacian(0.7, U) == fvm::source(Vec3{0, 1, 0}));
        require(!result.converged(), "zero x/z components hid the unconverged y equation");
    }
    std::cout << "time_history_test: repeated scalar/vector solves, Euler/BDF2, endTime passed\n";
}
