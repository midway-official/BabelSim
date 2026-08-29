#include "babelsim/legacy_taiho.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    ImportedTaihoMesh imported = readTaihoMesh("tests/data/taiho_5x5");
    const Mesh& mesh = imported.mesh;
    require(mesh.dimensions == std::array<Index, 3>{{3, 3, 1}},
            "legacy boundary cells were not removed");
    require(mesh.cellCount() == 9, "legacy interior cell count is incorrect");
    require(mesh.patches.size() == 5, "legacy patch grouping is incorrect");
    for (double volume : mesh.cell_volumes) {
        require(near(volume, 1.0), "legacy cell volume changed during import");
    }
    require(
        near(mesh.cell_centres[static_cast<std::size_t>(mesh.cellId(0, 0, 0))],
             {1.5, 1.5, 0.5}),
        "legacy interior cell centre is incorrect");

    int walls = 0;
    int inlets = 0;
    int outlets = 0;
    int symmetry = 0;
    for (const BoundaryPatch& patch : mesh.patches) {
        walls += patch.kind == PatchKind::Wall;
        inlets += patch.kind == PatchKind::Inlet;
        outlets += patch.kind == PatchKind::Outlet;
        symmetry += patch.kind == PatchKind::Symmetry;
    }
    require(
        walls == 1 && inlets == 1 && outlets == 1 && symmetry == 2,
        "legacy physical patch roles are incorrect");

    VectorField velocity(mesh, FieldLocation::Cell, "U");
    ScalarField pressure(mesh, FieldLocation::Cell, "p");
    applyImportedBoundaryConditions(imported, velocity, pressure);
    for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        const PatchKind kind = mesh.patches[static_cast<std::size_t>(patch)].kind;
        if (kind == PatchKind::Inlet) {
            require(
                velocity.boundary(patch).type == BoundaryType::FixedValue &&
                near(velocity.boundary(patch).value, {1.0, 0.0, 0.0}),
                "legacy inlet value is incorrect");
        } else if (kind == PatchKind::Outlet) {
            require(
                pressure.boundary(patch).type == BoundaryType::FixedValue,
                "legacy pressure outlet is incorrect");
        } else if (kind == PatchKind::Symmetry) {
            require(
                velocity.boundary(patch).type == BoundaryType::Symmetry,
                "extrusion patch is not a symmetry boundary");
        }
    }
    std::cout << "legacy_import_test: cells=" << mesh.cellCount()
              << " patches=" << mesh.patches.size() << '\n';
}

