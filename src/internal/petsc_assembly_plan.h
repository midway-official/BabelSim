#pragma once

#include "babelsim/discrete_equation.h"
#include "internal/petsc_index_map.h"

#include <petscmat.h>

#include <vector>

namespace babelsim::detail {

class PetscAssemblyPlan {
public:
    PetscAssemblyPlan(const Mesh& mesh, const PetscIndexMap& index_map);

    void cooPattern(std::vector<PetscInt>& rows, std::vector<PetscInt>& columns) const;
    PetscCount nonzeroContributions() const;

    void matrixValues(const ScalarDiscreteEquation& equation,
                      std::vector<PetscScalar>& values) const;
    void matrixValues(const VectorDiscreteEquation& equation,
                      std::vector<PetscScalar>& values) const;
    void rhsValues(const ScalarDiscreteEquation& equation, std::vector<PetscScalar>& values) const;
    void rhsValues(const VectorDiscreteEquation& equation, std::size_t component,
                   std::vector<PetscScalar>& values) const;

private:
    void validateEquation(const Mesh* equation_mesh, const std::vector<double>& diagonal,
                          const std::vector<double>& upper, const std::vector<double>& lower) const;
    void matrixValues(const std::vector<double>& diagonal, const std::vector<double>& upper,
                      const std::vector<double>& lower, std::vector<PetscScalar>& values) const;

    const Mesh* m_mesh;
    const PetscIndexMap* m_index_map;
    // COO coordinates are needed only when a PETSc Mat is created. Keep the
    // compact coefficient sources; derive coordinates into temporary buffers.
    std::vector<Index> m_upper_faces;
    std::vector<Index> m_lower_faces;
};

} // namespace babelsim::detail
