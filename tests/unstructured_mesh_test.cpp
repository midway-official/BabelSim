#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/mesh.h"
#include "babelsim/operators.h"

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
    return meshFromHexInput(std::move(vertices), std::move(cells),
                              {{"left", PatchKind::Generic}, {"right", PatchKind::Wall}},
                              std::move(boundaries));
}

Mesh pentagonalPrismMesh() {
    std::vector<Vec3> vertices;
    vertices.reserve(10);
    for (Index layer = 0; layer < 2; ++layer) {
        for (Index vertex = 0; vertex < 5; ++vertex) {
            const double angle = 2.0 * 3.14159265358979323846 * vertex / 5.0;
            vertices.push_back({std::cos(angle), std::sin(angle), static_cast<double>(layer)});
        }
    }
    std::vector<PolyhedralFaceSpec> faces{
        {{{0, 4, 3, 2, 1}}, 0, invalid_index, 0},
        {{{5, 6, 7, 8, 9}}, 0, invalid_index, 0},
    };
    for (Index vertex = 0; vertex < 5; ++vertex) {
        const Index next = (vertex + 1) % 5;
        faces.push_back({{{vertex, next, next + 5, vertex + 5}}, 0, invalid_index, 0});
    }
    return Mesh::polyhedral(std::move(vertices), std::move(faces),
                            {{"wall", PatchKind::Wall}});
}

Mesh tetrahedronMesh() {
    return Mesh::polyhedral(
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
        {
            {{{0, 2, 1}}, 0, invalid_index, 0},
            {{{0, 1, 3}}, 0, invalid_index, 0},
            {{{1, 2, 3}}, 0, invalid_index, 0},
            {{{2, 0, 3}}, 0, invalid_index, 0}},
        {{"boundary", PatchKind::Generic}});
}

Mesh triangularPrismMesh() {
    return Mesh::polyhedral(
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0},
         {0, 0, 1}, {1, 0, 1}, {0, 1, 1}},
        {
            {{{0, 2, 1}}, 0, invalid_index, 0},
            {{{3, 4, 5}}, 0, invalid_index, 0},
            {{{0, 1, 4, 3}}, 0, invalid_index, 0},
            {{{1, 2, 5, 4}}, 0, invalid_index, 0},
            {{{2, 0, 3, 5}}, 0, invalid_index, 0}},
        {{"boundary", PatchKind::Generic}});
}

Mesh squarePyramidMesh() {
    return Mesh::polyhedral(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}},
        {
            {{{0, 3, 2, 1}}, 0, invalid_index, 0},
            {{{0, 1, 4}}, 0, invalid_index, 0},
            {{{1, 2, 4}}, 0, invalid_index, 0},
            {{{2, 3, 4}}, 0, invalid_index, 0},
            {{{3, 0, 4}}, 0, invalid_index, 0}},
        {{"boundary", PatchKind::Generic}});
}

Mesh concavePrismMesh() {
    const std::vector<Vec3> base{
        {0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}, {1, 2, 0}, {0, 2, 0}};
    std::vector<Vec3> vertices = base;
    for (const Vec3& point : base) vertices.push_back({point.x, point.y, 1.0});
    std::vector<PolyhedralFaceSpec> faces{
        {{{0, 5, 4, 3, 2, 1}}, 0, invalid_index, 0},
        {{{6, 7, 8, 9, 10, 11}}, 0, invalid_index, 0}};
    for (Index vertex = 0; vertex < 6; ++vertex) {
        const Index next = (vertex + 1) % 6;
        faces.push_back({{{vertex, next, next + 6, vertex + 6}}, 0, invalid_index, 0});
    }
    return Mesh::polyhedral(std::move(vertices), std::move(faces),
                            {{"boundary", PatchKind::Generic}});
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
    const Mesh tetrahedron = tetrahedronMesh();
    const Mesh triangular_prism = triangularPrismMesh();
    const Mesh square_pyramid = squarePyramidMesh();
    const Mesh concave_prism = concavePrismMesh();
    require(tetrahedron.faceCount() == 4 && near(tetrahedron.cellVolume(0), 1.0 / 6.0),
            "tetrahedral polyhedral geometry is invalid");
    require(triangular_prism.faceCount() == 5 && near(triangular_prism.cellVolume(0), 0.5),
            "triangular-prism polyhedral geometry is invalid");
    require(square_pyramid.faceCount() == 5 && near(square_pyramid.cellVolume(0), 1.0 / 3.0),
            "pyramid polyhedral geometry is invalid");
    require(concave_prism.faceCount() == 8 && near(concave_prism.cellVolume(0), 3.0),
            "concave polyhedral face triangulation is invalid");
    for (const Mesh* mesh : {&tetrahedron, &triangular_prism, &square_pyramid, &concave_prism}) {
        ScalarField constant(*mesh, FieldLocation::Cell, "constant", 4.0);
        VectorField constant_gradient(*mesh, FieldLocation::Cell, "constantGradient");
        gradient(constant, constant_gradient, GradientMethod::LeastSquares);
        require(norm(detail::fieldData(constant_gradient)[0]) < 1e-12,
                "polyhedral least-squares gradient was not constant-preserving");
    }

    const Mesh polyhedron = pentagonalPrismMesh();
    require(polyhedron.cellCount() == 1 && polyhedron.faceCount() == 7,
            "variable-face polyhedron topology has the wrong size");
    require(polyhedron.cellFaces(0).size() == 7 && polyhedron.facePoints(0).size() == 5,
            "variable-face polyhedron connectivity was truncated");
    require(polyhedron.cellVertices(0).size() == 10 && near(polyhedron.cellVolume(0), 2.377641290737884),
            "variable-face polyhedron geometry is invalid");

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
        meshFromHexInput(
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
        meshFromHexInput(vertices, {cell, cell, cell}, {{"boundary", PatchKind::Generic}}, {});
    }, "non-manifold face was accepted");
    requireRejected([] {
        meshFromHexInput(
            {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
             {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
            {{{1, 0, 3, 2, 5, 4, 7, 6}}}, {{"boundary", PatchKind::Generic}}, {});
    }, "negative Hex orientation was accepted");
    requireRejected([] {
        meshFromHexInput(
            {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
             {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
            {{{0, 1, 2, 3, 4, 5, 6, 6}}}, {{"boundary", PatchKind::Generic}}, {});
    }, "degenerate Hex was accepted");
    requireRejected([] {
        Mesh::polyhedral(
            {{0, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 0}},
            {{{{0, 1, 2, 3}}, 0, invalid_index, 0}},
            {{"boundary", PatchKind::Generic}});
    }, "self-intersecting polygon was accepted");

    std::cout << "unstructured_mesh_test: arbitrary order, patches, geometry and rejection paths passed\n";
}
