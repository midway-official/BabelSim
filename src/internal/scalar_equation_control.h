#pragma once

namespace babelsim {

// 压力方程在没有固定压力边界时需要固定一个全局参考值。该离散控制仅服务于
// SIMPLE 私有 pEqn，不属于普通标量 PDE 的 Public Solver API。
struct ScalarEquationControl {
    bool fix_reference = false;
};

}  // babelsim 命名空间
