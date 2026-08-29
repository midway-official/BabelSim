#pragma once

#include "babelsim/mesh.h"

#include <algorithm>
#include <vector>

namespace babelsim {

// Coefficients are stored for every local cell/face so operators can use one
// indexing scheme for owned and ghost cells. SparseAssembly emits only owned
// rows; for an internal owner-neighbour face f, upper[f] is the owner-row
// coupling to the neighbour and lower[f] is the neighbour-row coupling to the
// owner. DistributedLinearSolver applies couplings to ghost neighbours in its
// halo matvec.
template <typename T>
struct Equation {
    explicit Equation(const Mesh& mesh_value)
        : mesh(&mesh_value),
          diagonal(static_cast<std::size_t>(mesh_value.cellCount()), 0.0),
          upper(static_cast<std::size_t>(mesh_value.faceCount()), 0.0),
          lower(static_cast<std::size_t>(mesh_value.faceCount()), 0.0),
          source(static_cast<std::size_t>(mesh_value.cellCount()), T{})
    {}

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

}  // namespace babelsim
