#pragma once

#include "babelsim/mesh.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace babelsim {

// 框架维护接口：FVM 后端的 LDU 存储，不是 Solver 使用的数学 EquationDefinition。
// 离散、装配、执行和分布式代数共享此唯一表示；不再保留含混的 Equation 兼容别名。
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

}  // babelsim 命名空间
