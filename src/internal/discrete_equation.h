#pragma once

#include "babelsim/mesh.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace babelsim::detail {

// 这是 FVM 后端的 LDU 存储，不是面向 Solver 的数学 Equation。它只在 Runtime、
// 离散和代数层之间传递，保留现有连续数组布局和分区面系数语义。
template <typename T>
struct DiscreteEquation {
    explicit DiscreteEquation(const Mesh& mesh_value)
        : mesh(&mesh_value),
          diagonal(static_cast<std::size_t>(mesh_value.cellCount()), 0.0),
          upper(static_cast<std::size_t>(mesh_value.faceCount()), 0.0),
          lower(static_cast<std::size_t>(mesh_value.faceCount()), 0.0),
          source(static_cast<std::size_t>(mesh_value.cellCount()), T{})
    {}

    DiscreteEquation(const Mesh&&) = delete;

    void validateStorage() const {
        if (mesh == nullptr ||
            diagonal.size() != static_cast<std::size_t>(mesh->cellCount()) ||
            source.size() != static_cast<std::size_t>(mesh->cellCount()) ||
            upper.size() != static_cast<std::size_t>(mesh->faceCount()) ||
            lower.size() != static_cast<std::size_t>(mesh->faceCount())) {
            throw std::logic_error("discrete equation storage invariant is violated");
        }
    }

    void reset() {
        std::fill(diagonal.begin(), diagonal.end(), 0.0);
        std::fill(upper.begin(), upper.end(), 0.0);
        std::fill(lower.begin(), lower.end(), 0.0);
        std::fill(source.begin(), source.end(), T{});
    }

    const Mesh* mesh;
    std::vector<double> diagonal;
    std::vector<double> upper;
    std::vector<double> lower;
    std::vector<T> source;
};

using ScalarDiscreteEquation = DiscreteEquation<double>;
using VectorDiscreteEquation = DiscreteEquation<Vec3>;

}  // babelsim::detail 命名空间
