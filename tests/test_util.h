#pragma once

#include "babelsim/vector.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

inline void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

inline bool near(double actual, double expected, double tolerance = 1e-11) {
    return std::abs(actual - expected) <= tolerance *
        std::max({1.0, std::abs(actual), std::abs(expected)});
}

inline bool near(
    const babelsim::Vec3& actual,
    const babelsim::Vec3& expected,
    double tolerance = 1e-11)
{
    return near(actual.x, expected.x, tolerance) &&
        near(actual.y, expected.y, tolerance) &&
        near(actual.z, expected.z, tolerance);
}

