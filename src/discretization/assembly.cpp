#include "babelsim/assembly.h"

#include <stdexcept>
#include <vector>

namespace babelsim {
namespace {

const Mesh& equationMesh(const Mesh* mesh) {
    if (mesh == nullptr) {
        throw std::invalid_argument("equation has no mesh");
    }
    return *mesh;
}

Eigen::Index coefficientPosition(
    const Eigen::SparseMatrix<double>& matrix,
    Index row,
    Index column)
{
    const auto* outer = matrix.outerIndexPtr();
    const auto* inner = matrix.innerIndexPtr();
    for (Eigen::Index position = outer[column];
         position < outer[column + 1]; ++position) {
        if (inner[position] == row) {
            return position;
        }
    }
    throw std::logic_error("sparse assembly pattern is incomplete");
}

template <typename T>
Eigen::SparseMatrix<double> assembleMatrixImpl(const Equation<T>& equation) {
    SparseAssembly assembly(equationMesh(equation.mesh));
    assembly.update(equation);
    return assembly.matrix();
}

}  // 匿名命名空间

SparseAssembly::SparseAssembly(const Mesh& mesh)
    : mesh_(&mesh),
      matrix_(mesh.ownedCellCount(), mesh.ownedCellCount()),
      diagonal_positions_(
          static_cast<std::size_t>(mesh.cellCount()), Eigen::Index{-1}),
      upper_positions_(
          static_cast<std::size_t>(mesh.faceCount()), Eigen::Index{-1}),
      lower_positions_(
          static_cast<std::size_t>(mesh.faceCount()), Eigen::Index{-1})
{
    mesh.validate();
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(
        static_cast<std::size_t>(mesh.ownedCellCount()) +
        2U * static_cast<std::size_t>(mesh.faceCount()));
    for (Index cell : mesh.owned_cells) {
        const Index row = mesh.ownedIndex(cell);
        entries.emplace_back(row, row, 1.0);
    }
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index neighbour = mesh.face_neighbour[f];
        if (neighbour == invalid_index) {
            continue;
        }
        const Index owner = mesh.face_owner[f];
        if (mesh.isOwned(owner) && mesh.isOwned(neighbour)) {
            entries.emplace_back(
                mesh.ownedIndex(owner), mesh.ownedIndex(neighbour), 1.0);
            entries.emplace_back(
                mesh.ownedIndex(neighbour), mesh.ownedIndex(owner), 1.0);
        }
    }

    matrix_.setFromTriplets(entries.begin(), entries.end());
    matrix_.makeCompressed();
    for (Index cell : mesh.owned_cells) {
        const Index row = mesh.ownedIndex(cell);
        diagonal_positions_[static_cast<std::size_t>(cell)] =
            coefficientPosition(matrix_, row, row);
    }
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index neighbour = mesh.face_neighbour[f];
        if (neighbour == invalid_index) {
            continue;
        }
        const Index owner = mesh.face_owner[f];
        if (mesh.isOwned(owner) && mesh.isOwned(neighbour)) {
            upper_positions_[f] = coefficientPosition(
                matrix_, mesh.ownedIndex(owner), mesh.ownedIndex(neighbour));
            lower_positions_[f] = coefficientPosition(
                matrix_, mesh.ownedIndex(neighbour), mesh.ownedIndex(owner));
        }
    }
}

void SparseAssembly::update(
    const Mesh* equation_mesh,
    const std::vector<double>& diagonal,
    const std::vector<double>& upper,
    const std::vector<double>& lower)
{
    if (equation_mesh != mesh_ ||
        diagonal.size() != diagonal_positions_.size() ||
        upper.size() != upper_positions_.size() ||
        lower.size() != lower_positions_.size()) {
        throw std::invalid_argument("equation coefficients do not match assembly mesh");
    }
    double* values = matrix_.valuePtr();
    for (Index cell : mesh_->owned_cells) {
        const auto c = static_cast<std::size_t>(cell);
        values[diagonal_positions_[c]] = diagonal[c];
    }
    for (std::size_t face = 0; face < upper.size(); ++face) {
        if (upper_positions_[face] < 0) {
            continue;
        }
        values[upper_positions_[face]] = upper[face];
        values[lower_positions_[face]] = lower[face];
    }
}

void SparseAssembly::update(const ScalarEquation& equation) {
    equation.validateStorage();
    update(equation.mesh, equation.diagonal, equation.upper, equation.lower);
}

void SparseAssembly::update(const VectorEquation& equation) {
    equation.validateStorage();
    update(equation.mesh, equation.diagonal, equation.upper, equation.lower);
}

Eigen::SparseMatrix<double> assembleMatrix(const ScalarEquation& equation) {
    return assembleMatrixImpl(equation);
}

Eigen::SparseMatrix<double> assembleMatrix(const VectorEquation& equation) {
    return assembleMatrixImpl(equation);
}

LinearSystem assemble(const ScalarEquation& equation) {
    LinearSystem system;
    system.A = assembleMatrix(equation);
    assembleSource(equation, system.b);
    return system;
}

void assembleSource(const ScalarEquation& equation, Eigen::VectorXd& result) {
    equation.validateStorage();
    const Mesh& mesh = equationMesh(equation.mesh);
    result.resize(mesh.ownedCellCount());
    for (Index cell : mesh.owned_cells) {
        result[mesh.ownedIndex(cell)] =
            equation.source[static_cast<std::size_t>(cell)];
    }
}

void assembleSource(
    const VectorEquation& equation,
    std::array<Eigen::VectorXd, 3>& result)
{
    equation.validateStorage();
    const Mesh& mesh = equationMesh(equation.mesh);
    for (auto& component : result) {
        component.resize(mesh.ownedCellCount());
    }
    for (Index cell : mesh.owned_cells) {
        const auto c = static_cast<std::size_t>(cell);
        const Index row = mesh.ownedIndex(cell);
        result[0][row] = equation.source[c].x;
        result[1][row] = equation.source[c].y;
        result[2][row] = equation.source[c].z;
    }
}

std::array<Eigen::VectorXd, 3> assembleSource(const VectorEquation& equation) {
    std::array<Eigen::VectorXd, 3> result;
    assembleSource(equation, result);
    return result;
}

}  // babelsim 命名空间
