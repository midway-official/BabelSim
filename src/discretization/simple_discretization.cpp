#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "internal/simple_discretization.h"
#include "babelsim/operators.h"

namespace babelsim {


namespace detail {

// Rhie–Chow 修正是领域专用离散，不放进通用 operators.cpp，也不持有 SIMPLE 状态。
void applyMomentumInterpolation(
    const ScalarField& pressure,
    const VectorField& pressure_gradient,
    const ScalarField& face_mobility,
    const VectorField& face_pressure_response,
    ScalarField& predicted_flux,
    DiffusionMethod method)
{
    const Mesh& mesh = pressure.mesh();
    const auto require_field = [&mesh](const auto& field, FieldLocation location) {
        field.validateStorage();
        if (&field.mesh() != &mesh || field.location() != location)
            throw std::invalid_argument("momentum interpolation field location or mesh is invalid");
    };
    require_field(pressure, FieldLocation::Cell);
    require_field(pressure_gradient, FieldLocation::Cell);
    require_field(face_mobility, FieldLocation::Face);
    require_field(face_pressure_response, FieldLocation::Face);
    require_field(predicted_flux, FieldLocation::Face);
    for (Index face : detail::meshData(mesh).owned_faces) {
        const std::size_t index = static_cast<std::size_t>(face);
        if (detail::meshData(mesh).face_neighbour[index] == invalid_index) continue;
        detail::fieldData(predicted_flux)[face] +=
            dot(detail::fieldData(face_pressure_response)[face], detail::meshData(mesh).face_area_vectors[index]) -
            detail::fieldData(face_mobility)[face] * integratedNormalGradient(
                pressure, pressure_gradient, face, method);
    }
}

}  // detail 命名空间

}  // babelsim 命名空间
