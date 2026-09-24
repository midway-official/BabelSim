#include "internal/petsc_assembly_plan.h"

#include "internal/mesh_access.h"

#include <stdexcept>

namespace babelsim::detail {

PetscAssemblyPlan::PetscAssemblyPlan(const Mesh& mesh, const PetscIndexMap& index_map)
    : m_mesh(&mesh), m_index_map(&index_map), m_owned_cell_order(meshData(mesh).owned_cells.begin(),
                                                                  meshData(mesh).owned_cells.end()) {
    mesh.validate();
    const auto& data = meshData(mesh);
    m_rows.reserve(m_owned_cell_order.size() + 2U * data.owned_faces.size());
    m_columns.reserve(m_rows.capacity());
    m_matrix_contributions.reserve(m_rows.capacity());
    for (Index cell : m_owned_cell_order) {
        const PetscInt row = index_map.row(cell);
        m_rows.push_back(row);
        m_columns.push_back(row);
        m_matrix_contributions.push_back({ContributionKind::Diagonal, cell});
    }
    // Each rank contributes every coefficient whose row it owns. Partitioned
    // internal faces are present on both endpoint ranks with local orientation
    // adjusted, so each off-process column is inserted exactly by its row owner.
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = data.face_owner[f];
        const Index neighbour = data.face_neighbour[f];
        if (neighbour == invalid_index) continue;
        if (isOwned(mesh, owner)) {
            m_rows.push_back(index_map.row(owner));
            m_columns.push_back(index_map.globalCell(neighbour));
            m_matrix_contributions.push_back({ContributionKind::Upper, face});
        }
        if (isOwned(mesh, neighbour)) {
            m_rows.push_back(index_map.row(neighbour));
            m_columns.push_back(index_map.globalCell(owner));
            m_matrix_contributions.push_back({ContributionKind::Lower, face});
        }
    }
    if (m_rows.size() != m_columns.size() || m_rows.size() != m_matrix_contributions.size())
        throw std::logic_error("PETSc COO assembly plan buffers disagree");
}

void PetscAssemblyPlan::validateEquation(
    const Mesh* equation_mesh, const std::vector<double>& diagonal,
    const std::vector<double>& upper, const std::vector<double>& lower) const {
    if (equation_mesh != m_mesh || diagonal.size() != static_cast<std::size_t>(m_mesh->cellCount()) ||
        upper.size() != static_cast<std::size_t>(m_mesh->faceCount()) ||
        lower.size() != static_cast<std::size_t>(m_mesh->faceCount()))
        throw std::invalid_argument("equation coefficients do not match PETSc assembly plan");
}

void PetscAssemblyPlan::matrixValues(
    const std::vector<double>& diagonal, const std::vector<double>& upper,
    const std::vector<double>& lower, std::vector<PetscScalar>& values) const {
    values.resize(m_matrix_contributions.size());
    for (std::size_t i = 0; i < m_matrix_contributions.size(); ++i) {
        const Contribution contribution = m_matrix_contributions[i];
        switch (contribution.kind) {
        case ContributionKind::Diagonal:
            values[i] = static_cast<PetscScalar>(diagonal[static_cast<std::size_t>(contribution.index)]);
            break;
        case ContributionKind::Upper:
            values[i] = static_cast<PetscScalar>(upper[static_cast<std::size_t>(contribution.index)]);
            break;
        case ContributionKind::Lower:
            values[i] = static_cast<PetscScalar>(lower[static_cast<std::size_t>(contribution.index)]);
            break;
        }
    }
}

void PetscAssemblyPlan::matrixValues(
    const ScalarDiscreteEquation& equation, std::vector<PetscScalar>& values) const {
    equation.validateStorage();
    validateEquation(equation.mesh, equation.diagonal, equation.upper, equation.lower);
    matrixValues(equation.diagonal, equation.upper, equation.lower, values);
}

void PetscAssemblyPlan::matrixValues(
    const VectorDiscreteEquation& equation, std::vector<PetscScalar>& values) const {
    equation.validateStorage();
    validateEquation(equation.mesh, equation.diagonal, equation.upper, equation.lower);
    matrixValues(equation.diagonal, equation.upper, equation.lower, values);
}

void PetscAssemblyPlan::rhsValues(const ScalarDiscreteEquation& equation,
                                  std::vector<PetscScalar>& values) const {
    equation.validateStorage();
    if (equation.mesh != m_mesh) throw std::invalid_argument("scalar RHS mesh mismatch");
    values.resize(m_owned_cell_order.size());
    for (std::size_t i = 0; i < m_owned_cell_order.size(); ++i)
        values[i] = static_cast<PetscScalar>(equation.source[static_cast<std::size_t>(m_owned_cell_order[i])]);
}

void PetscAssemblyPlan::rhsValues(const VectorDiscreteEquation& equation,
                                  std::size_t component,
                                  std::vector<PetscScalar>& values) const {
    equation.validateStorage();
    if (equation.mesh != m_mesh || component >= 3)
        throw std::invalid_argument("vector RHS mesh or component is invalid");
    values.resize(m_owned_cell_order.size());
    for (std::size_t i = 0; i < m_owned_cell_order.size(); ++i) {
        const Vec3& source = equation.source[static_cast<std::size_t>(m_owned_cell_order[i])];
        values[i] = static_cast<PetscScalar>(component == 0 ? source.x : component == 1 ? source.y : source.z);
    }
}

}  // namespace babelsim::detail
