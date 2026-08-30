#pragma once

#include "babelsim/field.h"

namespace babelsim {

// SIMPLE 从动量方程的离散主对角提取 rAU 所需的内部控制。它不属于 Public Solver API：
// 普通矢量 PDE 只写 fvm/fvc 方程并调用 solve()。
struct VectorEquationControl {
    double relaxation = 1.0;
    ScalarField* mobility = nullptr;
};

}  // babelsim 命名空间
