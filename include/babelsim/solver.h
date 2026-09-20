#pragma once

#include "babelsim/field.h"
#include "babelsim/math.h"
#include "babelsim/methods.h"
#include "babelsim/solver_control.h"

namespace babelsim {

// 当前问题的只读数值方法；算法可以据此选择修正次数，不接触运行后端。
// 引用仅在当前 Case/运行域生存期内有效。
const Methods& numericalMethods();
// Query only; the solver owns all output.
bool primaryProcess();

// 普通 Solver 的数学执行入口。此头文件不包含 RunTime、MPI 或离散矩阵实现。
// 面通量在每个控制体上的守恒误差。它是通用有限体积诊断量；不可压缩 SIMPLE 将
// 它命名为连续性残差。归约由 RunTime 完成，调用者不会接触 rank 或通信器。
struct FluxBalance {
    double l1 = 0.0;
    double l2 = 0.0;
    double maximum = 0.0;
    double relative = 0.0;
};

// Runtime 之外的 Solver API：显式量属于 math，隐式方程由 equ:: 组装并用 equ::solve
// 求解，收敛与守恒量属于 diagnostics。它们自动使用当前线程唯一活动的 RunTime，
// 因此 Solver 不需要在每个数学操作中传递执行对象。
namespace diagnostics {
double relativeChange(const VectorField& current, const VectorField& previous);
double relativeChange(const ScalarField& current, const ScalarField& previous);
double relativeMagnitude(const ScalarField& value, const ScalarField& reference);
FluxBalance fluxBalance(const ScalarField& face_flux);
bool all(bool local_condition);
}  // diagnostics 命名空间

}  // babelsim 命名空间
