#include "babelsim/assembly.h"
#include "internal/mesh_access.h"

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
Eigen::SparseMatrix<double> assembleMatrixImpl(const DiscreteEquation<T>& equation) {
    SparseAssembly assembly(equationMesh(equation.mesh));
    assembly.update(equation);
    return assembly.matrix();
}

}  // 匿名命名空间

SparseAssembly::SparseAssembly(const Mesh& mesh)
    : m_mesh(&mesh),
      m_matrix(detail::ownedCellCount(mesh), detail::ownedCellCount(mesh)),
      m_diagonal_positions(
          static_cast<std::size_t>(detail::ownedCellCount(mesh)), Eigen::Index{-1}),
      m_upper_positions(
          static_cast<std::size_t>(mesh.faceCount()), Eigen::Index{-1}),
      m_lower_positions(
          static_cast<std::size_t>(mesh.faceCount()), Eigen::Index{-1})
{
    mesh.validate();
    m_coupled_faces.reserve(detail::meshData(mesh).owned_faces.size());
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(
        static_cast<std::size_t>(detail::ownedCellCount(mesh)) +
        2U * detail::meshData(mesh).owned_faces.size());
    for (Index cell : detail::meshData(mesh).owned_cells) {
        const Index row = detail::ownedIndex(mesh, cell);
        entries.emplace_back(row, row, 1.0);
    }
    for (Index face : detail::meshData(mesh).owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        if (neighbour == invalid_index) {
            continue;
        }
        const Index owner = detail::meshData(mesh).face_owner[f];
        if (detail::isOwned(mesh, owner) && detail::isOwned(mesh, neighbour)) {
            m_coupled_faces.push_back(face);
            entries.emplace_back(
                detail::ownedIndex(mesh, owner), detail::ownedIndex(mesh, neighbour), 1.0);
            entries.emplace_back(
                detail::ownedIndex(mesh, neighbour), detail::ownedIndex(mesh, owner), 1.0);
        }
    }

    m_matrix.setFromTriplets(entries.begin(), entries.end());
    m_matrix.makeCompressed();
    for (Index cell : detail::meshData(mesh).owned_cells) {
        const Index row = detail::ownedIndex(mesh, cell);
        m_diagonal_positions[static_cast<std::size_t>(row)] =
            coefficientPosition(m_matrix, row, row);
    }
    for (Index face : m_coupled_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        if (neighbour == invalid_index) {
            continue;
        }
        const Index owner = detail::meshData(mesh).face_owner[f];
        if (detail::isOwned(mesh, owner) && detail::isOwned(mesh, neighbour)) {
            m_upper_positions[f] = coefficientPosition(
                m_matrix, detail::ownedIndex(mesh, owner), detail::ownedIndex(mesh, neighbour));
            m_lower_positions[f] = coefficientPosition(
                m_matrix, detail::ownedIndex(mesh, neighbour), detail::ownedIndex(mesh, owner));
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
    for (Index cell : detail::meshData(*m_mesh).owned_cells) {
        const auto c = static_cast<std::size_t>(cell);
        const auto row = static_cast<std::size_t>(detail::ownedIndex(*m_mesh, cell));
        values[m_diagonal_positions[row]] = diagonal[c];
    }
    for (Index face_index : m_coupled_faces) {
        const std::size_t face = static_cast<std::size_t>(face_index);
        values[m_upper_positions[face]] = upper[face];
        values[m_lower_positions[face]] = lower[face];
    }
}

void SparseAssembly::update(const ScalarDiscreteEquation& equation) {
    equation.validateStorage();
    update(equation.mesh, equation.diagonal, equation.upper, equation.lower);
}

void SparseAssembly::update(const VectorDiscreteEquation& equation) {
    equation.validateStorage();
    update(equation.mesh, equation.diagonal, equation.upper, equation.lower);
}

Eigen::SparseMatrix<double> assembleMatrix(const ScalarDiscreteEquation& equation) {
    return assembleMatrixImpl(equation);
}

Eigen::SparseMatrix<double> assembleMatrix(const VectorDiscreteEquation& equation) {
    return assembleMatrixImpl(equation);
}

LinearSystem assemble(const ScalarDiscreteEquation& equation) {
    LinearSystem system;
    system.A = assembleMatrix(equation);
    assembleSource(equation, system.b);
    return system;
}

void assembleSource(const ScalarDiscreteEquation& equation, Eigen::VectorXd& result) {
    equation.validateStorage();
    const Mesh& mesh = equationMesh(equation.mesh);
    result.resize(detail::ownedCellCount(mesh));
    for (Index cell : detail::meshData(mesh).owned_cells) {
        result[detail::ownedIndex(mesh, cell)] =
            equation.source[static_cast<std::size_t>(cell)];
    }
}

void assembleSource(
    const VectorDiscreteEquation& equation,
    std::array<Eigen::VectorXd, 3>& result)
{
    equation.validateStorage();
    const Mesh& mesh = equationMesh(equation.mesh);
    for (auto& component : result) {
        component.resize(detail::ownedCellCount(mesh));
    }
    for (Index cell : detail::meshData(mesh).owned_cells) {
        const auto c = static_cast<std::size_t>(cell);
        const Index row = detail::ownedIndex(mesh, cell);
        result[0][row] = equation.source[c].x;
        result[1][row] = equation.source[c].y;
        result[2][row] = equation.source[c].z;
    }
}

std::array<Eigen::VectorXd, 3> assembleSource(const VectorDiscreteEquation& equation) {
    std::array<Eigen::VectorXd, 3> result;
    assembleSource(equation, result);
    return result;
}

}  // babelsim 命名空间
