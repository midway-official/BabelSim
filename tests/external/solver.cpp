#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/solver.h"
#include "babelsim/equ.h"
#include <iostream>

using namespace babelsim;

// 外部方程驱动 Solver：不需要框架源码、专用 Case reader 或内置分派表。
SolverResult transport(Case& problem) {
    ScalarField& C = problem.scalarField("C");
    VectorField& U = problem.createVectorField("U", Vec3{});
    ScalarField& phi = problem.createFaceField("phi");
    phi = math::flux(U);
    const double D = problem.physics().nonnegative("diffusivity");
    const double Q = problem.physics().number("source");
    auto time=time::start(problem);
    auto old=time::history(C);
    auto A=equ::createEquation(C);
    while(time.value()<time.end()) {
        time.advance(); old.save(C, time.dt());
        A.reset(); equ::ddt(A,1.0,old); equ::div(A,phi);
        equ::laplacian(A,D, -1); equ::source(A,Q);
        if(!equ::solve(A,C).converged()) return SolverResult::notConverged();
        write(problem,time);
    }
    return SolverResult::completed();
}

const SolverRegistration transport_registration("transport_extension", transport);

// 外部算法驱动 Solver：两个 PDE 的定点耦合，收敛判断由公开全局诊断提供。
SolverResult coupled(Case& problem) {
    ScalarField& T = problem.createScalarField("T", 1.0);
    ScalarField& C = problem.createScalarField("C", 0.0);
    ScalarField& previousT = problem.createScalarField("previousT", 0.0);
    ScalarField& previousC = problem.createScalarField("previousC", 0.0);
    const double D = problem.physics().nonnegative("diffusivity");
    const double a = problem.physics().number("coupling");
    problem.output(T);
    problem.output(C);
    auto time=time::start(problem);
    auto oldT=time::history(T),oldC=time::history(C);
    auto A=equ::createEquation(T),B=equ::createEquation(C);
    while(time.value()<time.end()) {
        time.advance();
        oldT.save(T, time.dt()); oldC.save(C, time.dt());
        bool converged = false;
        for (int correction = 0; correction < 100; ++correction) {
            previousT.assign(T);
            previousC.assign(C);
            A.reset(); equ::ddt(A,1.0,oldT); equ::laplacian(A,D, -1); equ::source(A,a*C);
            if(!equ::solve(A,T).converged()) return SolverResult::notConverged();
            B.reset(); equ::ddt(B,1.0,oldC); equ::laplacian(B,D, -1); equ::source(B,a*T);
            if(!equ::solve(B,C).converged()) return SolverResult::notConverged();
            if (diagnostics::relativeChange(T, previousT) < 1e-12 &&
                diagnostics::relativeChange(C, previousC) < 1e-12) {
                converged = true;
                break;
            }
        }
        if (!converged) return SolverResult::notConverged();
        write(problem,time);
    }
    return SolverResult::completed();
}

const SolverRegistration coupled_registration("coupled_extension", coupled);

// 验收非均匀矢量场源、动量响应、压力规范和派生场输出；不是新的流动物理模型。
SolverResult vectorResponse(Case& problem) {
    VectorField& U = problem.createVectorField("U", Vec3{});
    VectorField& force = problem.createVectorField("force", Vec3{});
    ScalarField& p = problem.createScalarField("p", 0.0);
    ScalarField& rAU = problem.createScalarField("rAU", 0.0);
    ScalarField& energy = problem.createScalarField("energy", 0.0);
    const double strength = problem.physics().number("strength");
    force.evaluate([](Vec3 position) { return Vec3{1 + position.x, 2 + position.y, 3 + position.z}; });
    problem.validate();
    // 校验不应抢先关闭声明阶段；新的组合算法仍能声明自己的数学场。
    TensorField& stress = problem.createTensorField("stress", Tensor3{});
    stress.evaluate([](Vec3) { return Tensor3{{Vec3{1, 2, 3}, Vec3{4, 5, 6}, Vec3{7, 8, 9}}}; });
    problem.createFaceVectorField("faceU");
    problem.createFaceTensorField("faceStress");
    problem.output(U);
    problem.output(p);
    problem.output(rAU);
    problem.output(energy);
    problem.output(stress);
    auto time=time::start(problem);
    auto old=time::history(U);
    auto A=equ::createEquation(U);
    auto P=equ::createEquation(p);
    while(time.value()<time.end()) {
        time.advance(); old.save(U, time.dt());
        A.reset(); equ::ddt(A,1.0,old); equ::source(A,strength*force);
        rAU=A.volumeScaledInverseDiagonal();
        if(!equ::solve(A,U).converged()) return SolverResult::notConverged();
        P.reset(); equ::laplacian(P,1.0, -1); P.reference(0, 3.0);
        if(!equ::solve(P,p).converged()) return SolverResult::notConverged();
        math::subtract(rAU, math::grad(p), U);
        energy.evaluate(U, [](Vec3 velocity) { return 0.5*squaredNorm(velocity); });
        write(problem,time);
    }
    return SolverResult::completed();
}

const SolverRegistration vector_registration("vector_extension", vectorResponse);

int main(int argc, char* argv[]) {
    return runApplication(argc, argv, [](const char* message){ std::cerr << message << '\n'; });
}
