#pragma once

#include <cmath>
#include <stdexcept>

namespace babelsim {

// 时间控制只描述案例需要的时间区间和步长；历史场和并行同步属于 Runtime 内部。
struct TimeControl {
    double start_time = 0.0;
    double end_time = 1.0;
    double delta_t = 1e-3;

    void validate() const {
        if (!std::isfinite(start_time) || !std::isfinite(end_time) ||
            !std::isfinite(delta_t) || end_time < start_time || delta_t <= 0.0) {
            throw std::invalid_argument("time control is invalid");
        }
    }
};

}  // babelsim 命名空间
