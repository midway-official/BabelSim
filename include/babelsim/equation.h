#pragma once

// 兼容旧的框架级单元测试与内部后端。新 Physics Solver 不应包含本文件：
// `EquationDefinition` 位于 fvm.h，而真正的 LDU 载体只定义一次于 src/internal。
#include "internal/discrete_equation.h"

namespace babelsim {

template <typename T>
using Equation = detail::DiscreteEquation<T>;

using ScalarEquation = detail::ScalarDiscreteEquation;
using VectorEquation = detail::VectorDiscreteEquation;

}  // babelsim 命名空间
