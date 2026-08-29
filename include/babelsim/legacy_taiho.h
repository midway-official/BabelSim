#pragma once

#include "babelsim/field.h"

#include <filesystem>
#include <vector>

namespace babelsim {

struct ImportedPatchConditions {
    BoundaryCondition<Vec3> velocity;
    BoundaryCondition<double> pressure;
};

struct ImportedTaihoMesh {
    Mesh mesh;
    std::vector<ImportedPatchConditions> conditions;
};

// Imports TaihoCFD's outer boundary-cell format as a face-boundary mesh. The
// solved cells are the original [1,n-1) cells and the mesh is extruded once in z.
ImportedTaihoMesh readTaihoMesh(const std::filesystem::path& directory);

void applyImportedBoundaryConditions(
    const ImportedTaihoMesh& imported,
    VectorField& velocity,
    ScalarField& pressure);

}  // namespace babelsim

