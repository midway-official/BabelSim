#pragma once

#include "internal/mesh_access.h"
#include "internal/field_access.h"

namespace babelsim {
inline double boundaryNormalDistance(const Mesh& mesh, Index face) {
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    return dot(
        detail::meshData(mesh).face_centres[f] - detail::meshData(mesh).cell_centres[static_cast<std::size_t>(owner)],
        mesh.faceNormal(face));
}

inline double boundaryFaceValue(
    const ScalarField& field,
    Index face,
    double outward_flux = 0.0)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const Index patch = detail::meshData(mesh).face_patch[f];
    const auto& condition = field.boundary(patch);
    const double owner_value = detail::fieldData(field)[owner];
    const double distance = boundaryNormalDistance(mesh, face);
    switch (condition.type) {
        case BoundaryType::FixedValue:
            return condition.value;
        case BoundaryType::FixedGradient:
            return owner_value + condition.value * distance;
        case BoundaryType::ZeroGradient:
        case BoundaryType::Symmetry:
            return owner_value;
        case BoundaryType::InletOutlet:
            return outward_flux >= 0.0 ? owner_value : condition.value;
    }
    throw std::invalid_argument("unknown scalar boundary condition");
}

inline Vec3 boundaryFaceValue(
    const VectorField& field,
    Index face,
    double outward_flux = 0.0)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const Index patch = detail::meshData(mesh).face_patch[f];
    const auto& condition = field.boundary(patch);
    const Vec3 owner_value = detail::fieldData(field)[owner];
    const double distance = boundaryNormalDistance(mesh, face);
    switch (condition.type) {
        case BoundaryType::FixedValue:
            return condition.value;
        case BoundaryType::FixedGradient:
            return owner_value + condition.value * distance;
        case BoundaryType::ZeroGradient:
            return owner_value;
        case BoundaryType::InletOutlet:
            return outward_flux >= 0.0 ? owner_value : condition.value;
        case BoundaryType::Symmetry: {
            const Vec3 normal = mesh.faceNormal(face);
            return owner_value - dot(owner_value, normal) * normal;
        }
    }
    throw std::invalid_argument("unknown vector boundary condition");
}

}  // babelsim 命名空间
