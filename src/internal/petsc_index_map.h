#pragma once

#include "babelsim/mesh.h"
#include "babelsim/parallel.h"

#include <petscsys.h>

#include <vector>

namespace babelsim::detail {

// Maps BabelSim's local cell ordering and possibly sparse source IDs to PETSc's
// contiguous row ownership. Only ghost IDs referenced by this partition are
// exchanged; no rank builds a copy of the global cell table.
class PetscIndexMap {
public:
    PetscIndexMap(const Mesh& mesh, const ParallelContext& parallel);

    PetscInt globalCell(Index local_cell) const;
    PetscInt row(Index owned_cell) const;
    PetscInt localSize() const { return m_local_size; }
    PetscInt globalSize() const { return m_global_size; }
    PetscInt rowBegin() const { return m_row_begin; }
    PetscInt rowEnd() const { return m_row_end; }

private:
    const Mesh* m_mesh;
    PetscInt m_local_size = 0;
    PetscInt m_global_size = 0;
    PetscInt m_row_begin = 0;
    PetscInt m_row_end = 0;
    std::vector<PetscInt> m_local_to_global;
};

}  // namespace babelsim::detail
