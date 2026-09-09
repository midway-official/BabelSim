#include "internal/mesh_access.h"
#include "babelsim/mesh.h"

#include "test_util.h"

#include <iostream>
#include <stdexcept>
#include <utility>

using namespace babelsim;

namespace {

Mesh arbitraryHexMesh() {
    std::vector<Vec3> vertices{{
        {1, 1, 1}, {0, 0, 0}, {2, 0, 1}, {1, 0, 0}, {0, 1, 0}, {2, 1, 1},
        {1, 0, 1}, {0, 1, 1}, {2, 0, 0}, {1, 1, 0}, {0, 0, 1}, {2, 1, 0},
    }};
    // 单元和顶点均故意不按空间顺序存储；唯一要求是每个 Hex 的局部 VTK 顺序。
    std::vector<std::array<Index, 8>> cells{{
        {{3, 8, 11, 9, 6, 2, 5, 0}},
        {{1, 3, 9, 4, 10, 6, 0, 7}},
    }};
    std::vector<BoundaryFaceSpec> boundaries{{
        {{{1, 10, 7, 4}}, 0}, {{{1, 3, 6, 10}}, 0}, {{{4, 7, 0, 9}}, 0},
        {{{1, 4, 9, 3}}, 0}, {{{10, 6, 0, 7}}, 0},
        {{{8, 11, 5, 2}}, 1}, {{{3, 8, 2, 6}}, 1}, {{{9, 0, 5, 11}}, 1},
        {{{3, 9, 11, 8}}, 1}, {{{6, 2, 5, 0}}, 1},
    }};
    return Mesh::unstructured(std::move(vertices), std::move(cells),
                              {{"left", PatchKind::Generic}, {"right", PatchKind::Wall}},
                              std::move(boundaries));
}

template <typename Build>
void requireRejected(Build&& build, const char* message) {
    bool rejected = false;
    try {
        build();
    } catch (const std::exception&) {
        rejected = true;
    }
    require(rejected, message);
}

}  // namespace

int main() {
    const Mesh mesh = arbitraryHexMesh();
    require(mesh.cellCount() == 2 && mesh.globalCellCount() == 2, "explicit Hex cell IDs are invalid");
    require(mesh.vertexCount() == 12 && mesh.faceCount() == 11, "explicit Hex face matching failed");
    require(mesh.patchCount() == 2 && mesh.patchName(1) == "right" &&
                mesh.patchKind(1) == PatchKind::Wall,
            "explicit Hex patches were not preserved");
    Index internal = 0;
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        if (!mesh.boundaryFace(face)) {
            ++internal;
            require(near(mesh.faceArea(face), 1.0), "internal non-structured face area is incorrect");
        }
    }
    require(internal == 1 && near(mesh.cellVolume(0), 1.0) && near(mesh.cellVolume(1), 1.0),
            "explicit Hex geometry cache is invalid");

    requireRejected([] {
        std::vector<BoundaryFaceSpec> incomplete{{
            {{{0, 4, 7, 3}}, 0}, {{{1, 2, 6, 5}}, 0}, {{{0, 1, 5, 4}}, 0},
            {{{3, 7, 6, 2}}, 0}, {{{0, 3, 2, 1}}, 0},
        }};
        Mesh::unstructured(
            {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
             {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
            {{{0, 1, 2, 3, 4, 5, 6, 7}}}, {{"boundary", PatchKind::Generic}},
            std::move(incomplete));
    }, "missing exterior face was accepted");
    requireRejected([] {
        const std::vector<Vec3> vertices{{
            {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
        const std::array<Index, 8> cell{{0, 1, 2, 3, 4, 5, 6, 7}};
        Mesh::unstructured(vertices, {cell, cell, cell}, {{"boundary", PatchKind::Generic}}, {});
    }, "non-manifold face was accepted");
    requireRejected([] {
        Mesh::unstructured(
            {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
             {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
            {{{1, 0, 3, 2, 5, 4, 7, 6}}}, {{"boundary", PatchKind::Generic}}, {});
    }, "negative Hex orientation was accepted");
    requireRejected([] {
        Mesh::unstructured(
            {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
             {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
            {{{0, 1, 2, 3, 4, 5, 6, 6}}}, {{"boundary", PatchKind::Generic}}, {});
    }, "degenerate Hex was accepted");

    std::cout << "unstructured_mesh_test: arbitrary order, patches, geometry and rejection paths passed\n";
}
