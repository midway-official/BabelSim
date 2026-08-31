#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/field.h"
#include "internal/boundary_evaluation.h"

#include "test_util.h"

#include <iostream>
#include <utility>
#include <vector>

using namespace babelsim;

int main() {
    const Mesh mesh = Mesh::cartesian({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
    ScalarField scalar(mesh, FieldLocation::Cell, "phi", 2.0);
    require(scalar.size() == 1 && detail::fieldData(scalar) != nullptr, "cell field is not contiguous");
    scalar.setBoundary(
        static_cast<Index>(Side::XMin),
        BoundaryCondition<double>::fixedValue(5.0));
    scalar.setBoundary(
        static_cast<Index>(Side::XMax),
        BoundaryCondition<double>::fixedGradient(3.0));
    scalar.setBoundary(
        static_cast<Index>(Side::YMin),
        BoundaryCondition<double>::zeroGradient());
    scalar.setBoundary(
        static_cast<Index>(Side::YMax),
        BoundaryCondition<double>::inletOutlet(7.0));
    scalar.setBoundary(
        static_cast<Index>(Side::ZMin),
        BoundaryCondition<double>::symmetry());

    const Index xmin = detail::meshData(mesh).patches[static_cast<std::size_t>(Side::XMin)].faces.front();
    const Index xmax = detail::meshData(mesh).patches[static_cast<std::size_t>(Side::XMax)].faces.front();
    const Index ymin = detail::meshData(mesh).patches[static_cast<std::size_t>(Side::YMin)].faces.front();
    const Index ymax = detail::meshData(mesh).patches[static_cast<std::size_t>(Side::YMax)].faces.front();
    const Index zmin = detail::meshData(mesh).patches[static_cast<std::size_t>(Side::ZMin)].faces.front();
    require(near(boundaryFaceValue(scalar, xmin), 5.0), "fixedValue failed");
    require(near(boundaryFaceValue(scalar, xmax), 3.5), "fixedGradient failed");
    require(near(boundaryFaceValue(scalar, ymin), 2.0), "zeroGradient failed");
    require(near(boundaryFaceValue(scalar, ymax, -1.0), 7.0), "inlet failed");
    require(near(boundaryFaceValue(scalar, ymax, 1.0), 2.0), "outlet failed");
    require(near(boundaryFaceValue(scalar, zmin), 2.0), "scalar symmetry failed");

    ScalarField correction(mesh, FieldLocation::Cell, "phiPrime");
    const bool has_fixed_value = setHomogeneousCorrectionBoundaries(correction, scalar);
    require(has_fixed_value, "fixed scalar boundary was not detected");
    require(
        correction.boundary(static_cast<Index>(Side::XMin)).type == BoundaryType::FixedValue &&
        near(correction.boundary(static_cast<Index>(Side::XMin)).value, 0.0),
        "fixed correction boundary is not homogeneous");
    require(
        correction.boundary(static_cast<Index>(Side::ZMin)).type == BoundaryType::Symmetry,
        "symmetry correction boundary changed type");
    require(
        correction.boundary(static_cast<Index>(Side::XMax)).type == BoundaryType::ZeroGradient &&
        correction.boundary(static_cast<Index>(Side::YMax)).type == BoundaryType::ZeroGradient,
        "non-value correction boundary is not zeroGradient");

    VectorField vector(mesh, FieldLocation::Cell, "U", {1.0, 2.0, 3.0});
    vector.setBoundary(
        static_cast<Index>(Side::ZMax),
        BoundaryCondition<Vec3>::symmetry());
    const Index zmax = detail::meshData(mesh).patches[static_cast<std::size_t>(Side::ZMax)].faces.front();
    require(
        near(boundaryFaceValue(vector, zmax), {1.0, 2.0, 0.0}),
        "vector symmetry did not remove the normal component");

    ScalarField face(mesh, FieldLocation::Face, "flux");
    VectorField vertex(mesh, FieldLocation::Vertex, "point displacement");
    require(
        face.size() == static_cast<std::size_t>(mesh.faceCount()) &&
        vertex.size() == static_cast<std::size_t>(mesh.vertexCount()),
        "face or vertex field size is incorrect");

    std::vector<Vec3> skewed_points;
    for (Index k = 0; k <= 1; ++k) {
        for (Index j = 0; j <= 1; ++j) {
            for (Index i = 0; i <= 1; ++i) {
                skewed_points.push_back({
                    static_cast<double>(i) + 0.4 * j,
                    static_cast<double>(j),
                    static_cast<double>(k),
                });
            }
        }
    }
    const Mesh skewed = Mesh::structured({1, 1, 1}, std::move(skewed_points));
    ScalarField skewed_scalar(skewed, FieldLocation::Cell, "skewedPhi", 2.0);
    skewed_scalar.setBoundary(
        static_cast<Index>(Side::XMax),
        BoundaryCondition<double>::fixedGradient(3.0));
    const Index skewed_xmax =
        detail::meshData(skewed).patches[static_cast<std::size_t>(Side::XMax)].faces.front();
    const double normal_distance = boundaryNormalDistance(skewed, skewed_xmax);
    require(
        near(
            boundaryFaceValue(skewed_scalar, skewed_xmax),
            2.0 + 3.0 * normal_distance),
        "non-orthogonal fixedGradient used centre distance instead of normal distance");
    std::cout << "field_boundary_test: scalar/vector and mirror conditions passed\n";
}
