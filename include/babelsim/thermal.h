#pragma once

#include "babelsim/runtime.h"

namespace babelsim {

// 热物性是独立的物理值对象；热传导算法不关心它来自常数、表格还是更高层材料模型。
struct ThermalProperties {
    double density = 1.0;
    double heat_capacity = 1.0;
    double conductivity = 1.0;

    double volumetricHeatCapacity() const { return density * heat_capacity; }
    void validate() const;
};

struct HeatResult {
    SolveResult linear;
    int steps = 0;
    bool converged = false;
};

// Equation-driven Solver 的黄金模板。其实现只写热方程和时间循环：历史场、离散
// 系统、线性求解、并行同步和全局归约均由 RunTime 隐藏。
HeatResult solveTransientHeat(
    RunTime& run_time,
    ScalarField& temperature,
    const ThermalProperties& material,
    double volumetric_source = 0.0);

}  // babelsim 命名空间
