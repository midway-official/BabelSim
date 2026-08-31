#pragma once

#include "babelsim/field.h"

namespace babelsim {

// 执行层内部的方程控制，不包含任何 Solver 对象或算法生命周期。
// 数值算法可要求零参考约束、欠松弛或提取对角响应；普通 Solver 只使用 solve()。
struct VectorEquationControl {
    double relaxation = 1.0;
    ScalarField* mobility = nullptr;
};

}  // babelsim 命名空间
