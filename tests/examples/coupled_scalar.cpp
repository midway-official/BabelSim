#include "babelsim/case.h"
#include "babelsim/application.h"
#include "babelsim/solver.h"
#include "babelsim/equ.h"

namespace babelsim {

// 二次开发验收：单个普通函数描述双场耦合，没有新类、解析器、矩阵或通信代码。
// dT/dt = D laplacian(T) + a C；dC/dt = D laplacian(C) + a T。
SolverResult runCoupledScalar(Case& problem) {
    ScalarField& T = problem.load.scalarField("T");
    ScalarField& C = problem.load.scalarField("C");
    ScalarField& previous = problem.create.scalarField("previous", 0.0);
    ScalarField& previous_C = problem.create.scalarField("previousC", 0.0);
    const double D = problem.read.physics().nonnegative("diffusivity");
    const double a = problem.read.physics().number("coupling");
    const int corrections = problem.read.solution().integer("couplingIterations", 100);
    const double tolerance = problem.read.solution().number("couplingTolerance", 1e-12);

    auto time=time::makeStepper(problem);
    auto oldT=time::history(T), oldC=time::history(C);
    time::HistoryGroup histories;
    histories.add(T, oldT); histories.add(C, oldC);
    auto A=equ::createEquation(
        problem, "temperature", T, {"diffusion"});
    auto B=equ::createEquation(
        problem, "concentration", C, {"diffusion"});
    problem.finalizeInitialization();
    while(time.value()<time.end()) {
        time.advance();
        histories.save(time);
        bool converged = false;
        for (int correction = 0; correction < corrections; ++correction) {
            previous.assign(T);
            previous_C.assign(C);
            A.reset(); equ::ddt(A,1.0,oldT); equ::laplacian(A,D, -1); equ::source(A,a*C);
            if(!equ::solve(A,T).converged()) return SolverResult::notConverged();
            B.reset(); equ::ddt(B,1.0,oldC); equ::laplacian(B,D, -1); equ::source(B,a*T);
            if(!equ::solve(B,C).converged()) return SolverResult::notConverged();
            if (diagnostics::relativeChange(T, previous) <= tolerance &&
                diagnostics::relativeChange(C, previous_C) <= tolerance) {
                converged = true;
                break;
            }
        }
        if (!converged) return SolverResult::notConverged();
        write(problem,time);
    }
    return SolverResult::completed();
}

}  // babelsim 命名空间
