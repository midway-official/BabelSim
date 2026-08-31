#include "babelsim/case.h"
#include "babelsim/solver.h"

namespace babelsim {

// 二次开发验收：单个普通函数描述双场耦合，没有新类、解析器、矩阵或通信代码。
// dT/dt = D laplacian(T) + a C；dC/dt = D laplacian(C) + a T。
int runCoupledScalar(Case& problem) {
    ScalarField& T = problem.scalarField("T");
    ScalarField& C = problem.scalarField("C");
    ScalarField& previous = problem.scalarField("previous", 0.0);
    ScalarField& previous_C = problem.scalarField("previousC", 0.0);
    const double D = problem.properties().nonnegative("diffusivity");
    const double a = problem.properties().number("coupling");
    const int corrections = problem.solution().integer("couplingIterations", 100);
    const double tolerance = problem.solution().number("couplingTolerance", 1e-12);

    while (problem.loop()) {
        bool converged = false;
        for (int correction = 0; correction < corrections; ++correction) {
            previous.assign(T);
            previous_C.assign(C);
            if (!solve(fvm::ddt(T) == fvm::laplacian(D, T) + fvm::source(a, C)).converged()) return 2;
            if (!solve(fvm::ddt(C) == fvm::laplacian(D, C) + fvm::source(a, T)).converged()) return 2;
            if (diagnostics::relativeChange(T, previous) <= tolerance &&
                diagnostics::relativeChange(C, previous_C) <= tolerance) {
                converged = true;
                break;
            }
        }
        if (!converged) return 2;
    }
    return 0;
}

}  // babelsim 命名空间
