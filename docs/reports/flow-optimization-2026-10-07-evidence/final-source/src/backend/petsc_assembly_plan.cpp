#include "internal/petsc_assembly_plan.h"

#include "internal/mesh_access.h"

#include <stdexcept>

namespace babelsim::detail {

PetscAssemblyPlan::PetscAssemblyPlan(const Mesh& mesh, const PetscIndexMap& index_map)
    : m_mesh(&mesh), m_index_map(&index_map) {
    mesh.validate();
    const auto& data = meshData(mesh);
    std::size_t upper_count = 0, lower_count = 0;
    for (Index face : data.owned_faces) {
        if (mesh.neighbour(face) == invalid_index) continue;
        upper_count += isOwned(mesh, mesh.owner(face));
        lower_count += isOwned(mesh, mesh.neighbour(face));
    }
    m_upper_faces.reserve(upper_count);
    m_lower_faces.reserve(lower_count);
    for (Index face : data.owned_faces) {
        if (mesh.neighbour(face) == invalid_index) continue;
        if (isOwned(mesh, mesh.owner(face))) m_upper_faces.push_back(face);
        if (isOwned(mesh, mesh.neighbour(face))) m_lower_faces.push_back(face);
    }
}

PetscCount PetscAssemblyPlan::nonzeroContributions() const {
    return static_cast<PetscCount>(meshData(*m_mesh).owned_cells.size() +
                                   m_upper_faces.size() + m_lower_faces.size());
}

void PetscAssemblyPlan::cooPattern(std::vector<PetscInt>& rows,
                                   std::vector<PetscInt>& columns) const {
    const auto count = static_cast<std::size_t>(nonzeroContributions());
    rows.resize(count);
    columns.resize(count);
    std::size_t i = 0;
    for (Index cell : meshData(*m_mesh).owned_cells) {
        rows[i] = columns[i] = m_index_map->row(cell);
        ++i;
    }
    for (Index face : m_upper_faces) {
        rows[i] = m_index_map->row(m_mesh->owner(face));
        columns[i++] = m_index_map->globalCell(m_mesh->neighbour(face));
    }
    for (Index face : m_lower_faces) {
        rows[i] = m_index_map->row(m_mesh->neighbour(face));
        columns[i++] = m_index_map->globalCell(m_mesh->owner(face));
    }
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
    values.resize(static_cast<std::size_t>(nonzeroContributions()));
    std::size_t i = 0;
    for (Index cell : meshData(*m_mesh).owned_cells)
        values[i++] = static_cast<PetscScalar>(diagonal[static_cast<std::size_t>(cell)]);
    for (Index face : m_upper_faces)
        values[i++] = static_cast<PetscScalar>(upper[static_cast<std::size_t>(face)]);
    for (Index face : m_lower_faces)
        values[i++] = static_cast<PetscScalar>(lower[static_cast<std::size_t>(face)]);
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
    values.resize(meshData(*m_mesh).owned_cells.size());
    for (std::size_t i = 0; i < meshData(*m_mesh).owned_cells.size(); ++i)
        values[i] = static_cast<PetscScalar>(equation.source[static_cast<std::size_t>(meshData(*m_mesh).owned_cells[i])]);
}

void PetscAssemblyPlan::rhsValues(const VectorDiscreteEquation& equation,
                                  std::size_t component,
                                  std::vector<PetscScalar>& values) const {
    equation.validateStorage();
    if (equation.mesh != m_mesh || component >= 3)
        throw std::invalid_argument("vector RHS mesh or component is invalid");
    values.resize(meshData(*m_mesh).owned_cells.size());
    for (std::size_t i = 0; i < meshData(*m_mesh).owned_cells.size(); ++i) {
        const Vec3& source = equation.source[static_cast<std::size_t>(meshData(*m_mesh).owned_cells[i])];
        values[i] = static_cast<PetscScalar>(component == 0 ? source.x : component == 1 ? source.y : source.z);
    }
}

}  // namespace babelsim::detail
