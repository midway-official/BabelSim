#include "babelsim/case.h"
#include "babelsim/solver.h"
#include "babelsim/equ.h"

namespace babelsim {

// 二次开发验收：单个普通函数描述双场耦合，没有新类、解析器、矩阵或通信代码。
// dT/dt = D laplacian(T) + a C；dC/dt = D laplacian(C) + a T。
int runCoupledScalar(Case& problem) {
    ScalarField& T = problem.scalarField("T");
    ScalarField& C = problem.scalarField("C");
    ScalarField& previous = problem.scalarField("previous", 0.0);
    ScalarField& previous_C = problem.scalarField("previousC", 0.0);
    const double D = problem.physics().nonnegative("diffusivity");
    const double a = problem.physics().number("coupling");
    const int corrections = problem.solution().integer("couplingIterations", 100);
    const double tolerance = problem.solution().number("couplingTolerance", 1e-12);

    loadMethods(problem);
    auto time=enableTime(problem);
    auto oldT=math::history(T), oldC=math::history(C);
    auto A=equ::matrix(T), B=equ::matrix(C);
    while(time.value()<time.end()) {
        advance(time);
        math::saveOld(oldT,T,time.dt());
        math::saveOld(oldC,C,time.dt());
        bool converged = false;
        for (int correction = 0; correction < corrections; ++correction) {
            previous.assign(T);
            previous_C.assign(C);
            equ::clear(A); equ::ddt(A,1.0,oldT); equ::laplacian(A,D,-1.0); equ::source(A,a*C);
            if(!equ::solve(A,T).converged()) return 2;
            equ::clear(B); equ::ddt(B,1.0,oldC); equ::laplacian(B,D,-1.0); equ::source(B,a*T);
            if(!equ::solve(B,C).converged()) return 2;
            if (diagnostics::relativeChange(T, previous) <= tolerance &&
                diagnostics::relativeChange(C, previous_C) <= tolerance) {
                converged = true;
                break;
            }
        }
        if (!converged) return 2;
        write(problem,time);
    }
    return 0;
}

}  // babelsim 命名空间
