#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/equ.h"
#include "babelsim/runtime.h"
#include "test_util.h"

#include <iostream>

using namespace babelsim;

void checkHistory(TimeMethod method) {
    const Mesh mesh = makeHexBox({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
    ScalarField T(mesh, FieldLocation::Cell, "T", 0.0);
    VectorField U(mesh, FieldLocation::Cell, "U", {});
    RuntimeControl control;
    control.methods.time = method;
    control.time = {0.0, 0.3, 0.1};
    RunTime time = RunTime::forMesh(mesh, control);
    time::History<double> scalar_history = time::history(T);
    time::History<Vec3> vector_history = time::history(U);
    while (time.loop()) {
        // 时间层只由 History::save 推进，每个物理步恰好一次。
        scalar_history.save(T, time.deltaT());
        vector_history.save(U, time.deltaT());
        // 内迭代不是时间推进：重复求解同一方程，结果仍在同一个物理时间层。
        for (int correction = 0; correction < 3; ++correction) {
            auto scalar = testEquation(T);
            equ::ddt(scalar, 1.0, scalar_history);
            equ::source(scalar, 2.0);
            require(equ::solve(scalar).converged(), "scalar solve failed");
            auto vector = testEquation(U);
            equ::ddt(vector, 1.0, vector_history);
            equ::source(vector, Vec3{1, 2, 3});
            require(equ::solve(vector).converged(), "vector solve failed");
            require(near(detail::fieldData(T)[0], 2.0 * time.time(), 1e-12), "inner solve advanced scalar history");
            require(near(detail::fieldData(U)[0], Vec3{time.time(), 2*time.time(), 3*time.time()}, 1e-12),
                    "inner solve advanced vector history");
        }
    }
    require(time.step() == 3 && !time.loop(), "time loop ran extra steps");
}

int main() {
    checkHistory(TimeMethod::Euler);
    checkHistory(TimeMethod::BDF2);
    {
        const Mesh mesh = makeHexBox({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
        ScalarField T(mesh, FieldLocation::Cell, "T", 0.0);
        RuntimeControl control;
        control.methods.time = TimeMethod::Euler;
        control.time = {0.0, 0.25, 0.1};
        RunTime time = RunTime::forMesh(mesh, control);
        time::History<double> history = time::history(T);
        while (time.loop()) {
            history.save(T, time.deltaT());
            auto equation = testEquation(T);
            equ::ddt(equation, 1.0, history);
            equ::source(equation, 1.0);
            require(equ::solve(equation).converged(), "short-step solve failed");
            require(near(detail::fieldData(T)[0], time.time(), 1e-12), "short final step used wrong deltaT");
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
        // RunTime::loop 只做等步长时间推进；变步长 BDF2 由显式 time::History 承担，
        // 因此这条路必须整体拒绝，而不是静默按等步长系数推进。
        const Mesh mesh = makeHexBox({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
        RunTime time = RunTime::forMesh(mesh, control);
        time.loop();
    } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "legacy loop must not silently use uniform-step BDF2 coefficients");
    {
        const Mesh mesh = makeHexBox({6, 5, 1}, {0, 0, 0}, {1, 1, 1});
        VectorField U(mesh, FieldLocation::Cell, "U");
        for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch)
            U.boundary(patch) = fixedValue(Vec3{});
        RuntimeControl control;
        RunTime time = RunTime::forMesh(mesh, control);
        auto equationControl = testEquationControl("velocity");
        equationControl.linear.max_iterations = 1;
        equationControl.linear.absolute_tolerance = 1e-30;
        equationControl.linear.relative_tolerance = 1e-25;
        auto equation = equ::createEquation(U, equationControl);
        equ::laplacian(equation, 0.7, -1.0);  // 左端 -div(0.7 grad U)
        equ::source(equation, Vec3{0, 1, 0});
        const SolveResult result = equ::solve(equation);
        require(!result.converged(), "zero x/z components hid the unconverged y equation");
    }
    std::cout << "time_history_test: repeated scalar/vector solves, Euler/BDF2, endTime passed\n";
}
