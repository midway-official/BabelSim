#include "babelsim/mesh_io.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    const Mesh mesh = readMeshFile("tests/data/babelsim_channel.mesh");
    require(mesh.dimensions == std::array<Index, 3>{{3, 3, 1}},
            "native mesh dimensions are incorrect");
    require(mesh.cellCount() == 9, "native mesh cell count is incorrect");
    require(mesh.patches.size() == 6, "native mesh patch count is incorrect");
    for (double volume : mesh.cell_volumes) {
        require(near(volume, 1.0), "native mesh cell volume is incorrect");
    }
    require(
        near(mesh.cell_centres[static_cast<std::size_t>(mesh.cellId(0, 0, 0))],
             {0.5, 0.5, 0.5}),
        "native mesh cell centre is incorrect");

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
        walls == 2 && inlets == 1 && outlets == 1 && symmetry == 2,
        "native physical patch roles are incorrect");

    const Mesh vertices = readMeshFile("tests/data/babelsim_vertices.mesh");
    require(vertices.cellCount() == 1 && near(vertices.cell_volumes[0], 1.0),
            "native explicit-vertex geometry is incorrect");

    std::cout << "mesh_file_test: cells=" << mesh.cellCount()
              << " patches=" << mesh.patches.size() << '\n';
}
