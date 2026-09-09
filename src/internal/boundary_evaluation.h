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

template <typename T>
inline T boundaryFaceValue(const Field<T>& field, Index face,
    double outward_flux = std::numeric_limits<double>::quiet_NaN())
{
    return detail::FieldAccess::trace(field, face, outward_flux);
}

}  // babelsim namespace
