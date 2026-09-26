#pragma once

#include "babelsim/discrete_equation.h"
#include "internal/petsc_index_map.h"

#include <petscmat.h>

#include <vector>

namespace babelsim::detail {

class PetscAssemblyPlan {
public:
    PetscAssemblyPlan(const Mesh& mesh, const PetscIndexMap& index_map);

    const std::vector<PetscInt>& rows() const { return m_rows; }
    const std::vector<PetscInt>& columns() const { return m_columns; }
    PetscCount nonzeroContributions() const { return static_cast<PetscCount>(m_rows.size()); }

    void matrixValues(const ScalarDiscreteEquation& equation, std::vector<PetscScalar>& values) const;
    void matrixValues(const VectorDiscreteEquation& equation, std::vector<PetscScalar>& values) const;
    void rhsValues(const ScalarDiscreteEquation& equation, std::vector<PetscScalar>& values) const;
    void rhsValues(const VectorDiscreteEquation& equation, std::size_t component,
                  std::vector<PetscScalar>& values) const;

private:
    enum class ContributionKind : unsigned char { Diagonal, Upper, Lower };
    struct Contribution {
        ContributionKind kind;
        Index index;
    };

    void validateEquation(const Mesh* equation_mesh, const std::vector<double>& diagonal,
                          const std::vector<double>& upper, const std::vector<double>& lower) const;
    void matrixValues(const std::vector<double>& diagonal, const std::vector<double>& upper,
                      const std::vector<double>& lower, std::vector<PetscScalar>& values) const;

    const Mesh* m_mesh;
    const PetscIndexMap* m_index_map;
    std::vector<PetscInt> m_rows;
    std::vector<PetscInt> m_columns;
    std::vector<Contribution> m_matrix_contributions;
    std::vector<Index> m_owned_cell_order;
};

}  // namespace babelsim::detail
