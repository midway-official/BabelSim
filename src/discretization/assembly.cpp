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
    : m_mesh(&mesh),
      m_matrix(mesh.ownedCellCount(), mesh.ownedCellCount()),
      m_diagonal_positions(
          static_cast<std::size_t>(mesh.ownedCellCount()), Eigen::Index{-1}),
      m_upper_positions(
          static_cast<std::size_t>(mesh.faceCount()), Eigen::Index{-1}),
      m_lower_positions(
          static_cast<std::size_t>(mesh.faceCount()), Eigen::Index{-1})
{
    mesh.validate();
    m_coupled_faces.reserve(mesh.owned_faces.size());
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(
        static_cast<std::size_t>(mesh.ownedCellCount()) +
        2U * mesh.owned_faces.size());
    for (Index cell : mesh.owned_cells) {
        const Index row = mesh.ownedIndex(cell);
        entries.emplace_back(row, row, 1.0);
    }
    for (Index face : mesh.owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index neighbour = mesh.face_neighbour[f];
        if (neighbour == invalid_index) {
            continue;
        }
        const Index owner = mesh.face_owner[f];
        if (mesh.isOwned(owner) && mesh.isOwned(neighbour)) {
            m_coupled_faces.push_back(face);
            entries.emplace_back(
                mesh.ownedIndex(owner), mesh.ownedIndex(neighbour), 1.0);
            entries.emplace_back(
                mesh.ownedIndex(neighbour), mesh.ownedIndex(owner), 1.0);
        }
    }

    m_matrix.setFromTriplets(entries.begin(), entries.end());
    m_matrix.makeCompressed();
    for (Index cell : mesh.owned_cells) {
        const Index row = mesh.ownedIndex(cell);
        m_diagonal_positions[static_cast<std::size_t>(row)] =
            coefficientPosition(m_matrix, row, row);
    }
    for (Index face : m_coupled_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index neighbour = mesh.face_neighbour[f];
        if (neighbour == invalid_index) {
            continue;
        }
        const Index owner = mesh.face_owner[f];
        if (mesh.isOwned(owner) && mesh.isOwned(neighbour)) {
            m_upper_positions[f] = coefficientPosition(
                m_matrix, mesh.ownedIndex(owner), mesh.ownedIndex(neighbour));
            m_lower_positions[f] = coefficientPosition(
                m_matrix, mesh.ownedIndex(neighbour), mesh.ownedIndex(owner));
        }
    }
}

void SparseAssembly::update(
    const Mesh* equation_mesh,
    const std::vector<double>& diagonal,
    const std::vector<double>& upper,
    const std::vector<double>& lower)
{
    if (equation_mesh != m_mesh ||
        diagonal.size() != static_cast<std::size_t>(m_mesh->cellCount()) ||
        upper.size() != m_upper_positions.size() ||
        lower.size() != m_lower_positions.size()) {
        throw std::invalid_argument("equation coefficients do not match assembly mesh");
    }
    double* values = m_matrix.valuePtr();
    for (Index cell : m_mesh->owned_cells) {
        const auto c = static_cast<std::size_t>(cell);
        const auto row = static_cast<std::size_t>(m_mesh->ownedIndex(cell));
        values[m_diagonal_positions[row]] = diagonal[c];
    }
    for (Index face_index : m_coupled_faces) {
        const std::size_t face = static_cast<std::size_t>(face_index);
        values[m_upper_positions[face]] = upper[face];
        values[m_lower_positions[face]] = lower[face];
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
