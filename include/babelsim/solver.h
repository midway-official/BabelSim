#pragma once

#include "babelsim/fvm.h"
#include "babelsim/solver_control.h"

namespace babelsim {

// 普通 Solver 的数学执行入口。此头文件不包含 RunTime、MPI 或离散矩阵实现。
// 面通量在每个控制体上的守恒误差。它是通用有限体积诊断量；不可压缩 SIMPLE 将
// 它命名为连续性残差。归约由 RunTime 完成，调用者不会接触 rank 或通信器。
struct FluxBalance {
    double l1 = 0.0;
    double l2 = 0.0;
    double maximum = 0.0;
    double relative = 0.0;
};

// Runtime 之外的 Solver API：显式量属于 fvc，收敛与守恒量属于 diagnostics，
// 隐式方程由 solve() 处理。它们自动使用当前线程唯一活动的 RunTime，因此 Solver
// 不需要在每个数学操作中传递执行对象。
SolveResult solve(const ScalarEquationDefinition& equation);
// 标量和矢量方程都返回一个方程级结果；分量细节仅留在内部数值诊断。
SolveResult solve(const VectorEquationDefinition& equation);

namespace fvc {
void evaluate(ScalarGradient operation, VectorField& result);
void evaluate(VectorGradient operation, TensorField& result);
void evaluate(FaceFlux operation, ScalarField& result);
void evaluate(FaceDivergence operation, ScalarField& result);
void evaluate(VectorDivergence operation, ScalarField& result);
void evaluate(ScalarConvection operation, ScalarField& result);
void evaluate(VectorConvection operation, VectorField& result);
void evaluate(ScalarInterpolation operation, ScalarField& result);
void evaluate(VectorInterpolation operation, VectorField& result);
void evaluate(ScalarReconstruction operation, ScalarField& result);
void evaluate(VectorReconstruction operation, VectorField& result);
void evaluate(ScalarLaplacian operation, ScalarField& result);
void subtract(
    const ScalarField& coefficient,
    ScalarGradient operation,
    VectorField& target);
void subtract(ScalarDiffusionFlux operation, ScalarField& target);
}  // fvc 命名空间

namespace diagnostics {
double relativeChange(const VectorField& current, const VectorField& previous);
double relativeChange(const ScalarField& current, const ScalarField& previous);
double relativeMagnitude(const ScalarField& value, const ScalarField& reference);
FluxBalance fluxBalance(const ScalarField& face_flux);
bool all(bool local_condition);
}  // diagnostics 命名空间

}  // babelsim 命名空间
