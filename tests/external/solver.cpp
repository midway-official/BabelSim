#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/solver.h"

using namespace babelsim;

// 外部方程驱动 Solver：不需要框架源码、专用 Case reader 或内置分派表。
int transport(Case& problem) {
    ScalarField& C = problem.scalarField("C");
    VectorField& U = problem.vectorField("U", Vec3{});
    ScalarField& phi = problem.faceFlux("phi", U);
    const double D = problem.physics().nonnegative("diffusivity");
    const double Q = problem.physics().number("source");
    while (problem.loop()) {
        if (!solve(eqn::ddt(C) + eqn::div(phi, C) ==
                   eqn::laplacian(D, C) + Q).converged()) return 2;
    }
    return 0;
}

const SolverRegistration transport_registration("transport_extension", transport);

// 外部算法驱动 Solver：两个 PDE 的定点耦合，收敛判断由公开全局诊断提供。
int coupled(Case& problem) {
    ScalarField& T = problem.scalarField("T", 1.0);
    ScalarField& C = problem.scalarField("C", 0.0);
    ScalarField& previousT = problem.scalarField("previousT", 0.0);
    ScalarField& previousC = problem.scalarField("previousC", 0.0);
    const double D = problem.physics().nonnegative("diffusivity");
    const double a = problem.physics().number("coupling");
    problem.output(T);
    problem.output(C);
    while (problem.loop()) {
        bool converged = false;
        for (int correction = 0; correction < 100; ++correction) {
            previousT.assign(T);
            previousC.assign(C);
            if (!solve(eqn::ddt(T) == eqn::laplacian(D, T) + eqn::source(a, C)).converged()) return 2;
            if (!solve(eqn::ddt(C) == eqn::laplacian(D, C) + eqn::source(a, T)).converged()) return 2;
            if (diagnostics::relativeChange(T, previousT) < 1e-12 &&
                diagnostics::relativeChange(C, previousC) < 1e-12) {
                converged = true;
                break;
            }
        }
        if (!converged) return 2;
    }
    return 0;
}

const SolverRegistration coupled_registration("coupled_extension", coupled);

// 验收非均匀矢量场源、动量响应、压力规范和派生场输出；不是新的流动物理模型。
int vectorResponse(Case& problem) {
    VectorField& U = problem.vectorField("U", Vec3{});
    VectorField& force = problem.vectorField("force", Vec3{});
    ScalarField& p = problem.scalarField("p", 0.0);
    ScalarField& rAU = problem.scalarField("rAU", 0.0);
    ScalarField& energy = problem.scalarField("energy", 0.0);
    const double strength = problem.physics().number("strength");
    force.evaluate([](Vec3 position) { return Vec3{1 + position.x, 2 + position.y, 3 + position.z}; });
    problem.validate();
    // 校验不应抢先关闭声明阶段；新的组合算法仍能声明自己的数学场。
    TensorField& stress = problem.tensorField("stress", Tensor3{});
    stress.evaluate([](Vec3) { return Tensor3{{Vec3{1, 2, 3}, Vec3{4, 5, 6}, Vec3{7, 8, 9}}}; });
    problem.faceVectorField("faceU");
    problem.faceTensorField("faceStress");
    problem.output(U);
    problem.output(p);
    problem.output(rAU);
    problem.output(energy);
    problem.output(stress);
    while (problem.loop()) {
        if (!solveWithResponse(eqn::ddt(U) == eqn::source(strength, force), rAU).converged()) return 2;
        if (!solve(-eqn::laplacian(1.0, p) == 0.0, referenceValue(3.0)).converged()) return 2;
        math::subtract(rAU, math::grad(p), U);
        energy.evaluate(U, [](Vec3 velocity) { return 0.5*squaredNorm(velocity); });
    }
    return 0;
}

const SolverRegistration vector_registration("vector_extension", vectorResponse);

int main(int argc, char* argv[]) {
    return runApplication(argc, argv);
}
