#pragma once

#include "babelsim/runtime.h"

namespace babelsim {

// 标量输运是复用通用 FVM API 的最小非热学示例：储存、对流、扩散和体源均没有
// 专用矩阵或并行实现。
struct ScalarTransportResult {
    SolveResult linear;
    int steps = 0;
    bool converged = false;
};

ScalarTransportResult solveTransientScalarTransport(
    RunTime& run_time,
    ScalarField& scalar,
    const ScalarField& face_flux,
    double storage = 1.0,
    double diffusivity = 0.0,
    double source = 0.0);

}  // babelsim 命名空间
