#pragma once

#include "babelsim/mesh.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace babelsim {

// 系数为每个局部单元/面存储，使算子可用同一索引方案访问自有与幽灵单元。
// SparseAssembly 仅生成 owned 行；对内部 owner-neighbour 面 f，upper[f] 是
// owner 行到 neighbour 行的耦合，lower[f] 是 neighbour 行到 owner 行的耦合。
// DistributedLinearSolver 在 halo 矩阵向量乘中施加到幽灵 neighbour 的耦合。
template <typename T>
struct Equation {
    explicit Equation(const Mesh& mesh_value)
        : mesh(&mesh_value),
          diagonal(static_cast<std::size_t>(mesh_value.cellCount()), 0.0),
          upper(static_cast<std::size_t>(mesh_value.faceCount()), 0.0),
          lower(static_cast<std::size_t>(mesh_value.faceCount()), 0.0),
          source(static_cast<std::size_t>(mesh_value.cellCount()), T{})
    {}

    Equation(const Mesh&&) = delete;

    void validateStorage() const {
        if (mesh == nullptr ||
            diagonal.size() != static_cast<std::size_t>(mesh->cellCount()) ||
            source.size() != static_cast<std::size_t>(mesh->cellCount()) ||
            upper.size() != static_cast<std::size_t>(mesh->faceCount()) ||
            lower.size() != static_cast<std::size_t>(mesh->faceCount())) {
            throw std::logic_error("equation storage invariant is violated");
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

using ScalarEquation = Equation<double>;
using VectorEquation = Equation<Vec3>;

}  // babelsim 命名空间
