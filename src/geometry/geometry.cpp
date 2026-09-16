#include "babelsim/geometry.h"

#include "internal/field_access.h"

namespace babelsim::geometry {
namespace {

template <typename T, typename Value>
Field<T> values(const Mesh& mesh, FieldLocation location, const char* name,
                Value value) {
    Field<T> result(mesh, location, name);
    T* output = detail::fieldData(result);
    const Index count = location == FieldLocation::Cell ? mesh.cellCount() : mesh.faceCount();
    for (Index index = 0; index < count; ++index)
        output[index] = value(index);
    detail::markHaloValid(result);
    if (location == FieldLocation::Cell)
        result.useCalculatedBoundary();
    return result;
}

}  // namespace

ScalarField cellVolumes(const Mesh& mesh) {
    return values<double>(mesh, FieldLocation::Cell, "cellVolumes",
                          [&mesh](Index cell) { return mesh.cellVolume(cell); });
}

VectorField cellCentres(const Mesh& mesh) {
    return values<Vec3>(mesh, FieldLocation::Cell, "cellCentres",
                        [&mesh](Index cell) { return mesh.cellCentre(cell); });
}

ScalarField faceAreas(const Mesh& mesh) {
    return values<double>(mesh, FieldLocation::Face, "faceAreas",
                          [&mesh](Index face) { return mesh.faceArea(face); });
}

VectorField faceAreaVectors(const Mesh& mesh) {
    return values<Vec3>(mesh, FieldLocation::Face, "faceAreaVectors",
                        [&mesh](Index face) { return mesh.faceAreaVector(face); });
}

VectorField faceCentres(const Mesh& mesh) {
    return values<Vec3>(mesh, FieldLocation::Face, "faceCentres",
                        [&mesh](Index face) { return mesh.faceCentre(face); });
}

VectorField faceUnitNormals(const Mesh& mesh) {
    return values<Vec3>(mesh, FieldLocation::Face, "faceUnitNormals",
                        [&mesh](Index face) { return mesh.faceNormal(face); });
}

}  // namespace babelsim::geometry
