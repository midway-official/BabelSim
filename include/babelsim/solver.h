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
// 方程级数值控制，不含矩阵、工作区或通信参数。referenceValue 是零空间的定值规范：
// 在框架选定的固定参考单元设置值，适用于具有常数零空间且满足相容条件的标量方程。
struct EquationControl {
    double relaxation = 1.0;
    bool fix_reference = false;
    double reference_value = 0.0;
};

inline EquationControl relaxed(double factor) { return {factor, false, 0.0}; }
inline EquationControl referenceValue(double value) { return {1.0, true, value}; }

[[nodiscard]] SolveResult solve(
    const ScalarEquationDefinition& equation, EquationControl control = {});
// 标量和矢量方程都返回一个方程级结果；分量细节仅留在内部数值诊断。
[[nodiscard]] SolveResult solve(
    const VectorEquationDefinition& equation, EquationControl control = {});

// 求解矢量方程并返回对体源的对角响应 V/aP；用于 SIMPLE 等算法，而不是暴露 aP 存储。
// aP 为当前缩放行的对角，沿用 SIMPLE 约定；欠松弛时原始体源本身另乘 relaxation。
// response 必须是独立的同网格 cell 场，不得覆盖方程的输入系数或压力。
[[nodiscard]] SolveResult solveWithResponse(
    const VectorEquationDefinition& equation, ScalarField& response,
    EquationControl control = {});


namespace diagnostics {
double relativeChange(const VectorField& current, const VectorField& previous);
double relativeChange(const ScalarField& current, const ScalarField& previous);
double relativeMagnitude(const ScalarField& value, const ScalarField& reference);
FluxBalance fluxBalance(const ScalarField& face_flux);
bool all(bool local_condition);
}  // diagnostics 命名空间

}  // babelsim 命名空间
