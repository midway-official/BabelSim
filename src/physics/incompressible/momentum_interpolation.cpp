#include "babelsim/incompressible_operators.h"

#include <stdexcept>

namespace babelsim {

void MomentumInterpolation::apply(
    const Mesh& mesh,
    const VectorField& velocity,
    const ScalarField& pressure,
    const ScalarField& mobility,
    const VectorField& pressure_gradient,
    ScalarField& face_flux)
{
    if (&velocity.mesh() != &mesh || &pressure.mesh() != &mesh ||
        &mobility.mesh() != &mesh || &pressure_gradient.mesh() != &mesh ||
        &face_flux.mesh() != &mesh || velocity.location() != FieldLocation::Cell ||
        pressure.location() != FieldLocation::Cell ||
        mobility.location() != FieldLocation::Cell ||
        pressure_gradient.location() != FieldLocation::Cell ||
        face_flux.location() != FieldLocation::Face) {
        throw std::invalid_argument("momentum interpolation fields do not match mesh");
    }
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        if (neighbour == invalid_index) {
            const Vec3 boundary_velocity = boundaryFaceValue(
                velocity, face, face_flux[face]);
            face_flux[face] = dot(boundary_velocity, mesh.face_area_vectors[f]);
            continue;
        }
        const double weight = mesh.face_owner_weights[f];
        const double face_mobility =
            weight * mobility[owner] + (1.0 - weight) * mobility[neighbour];
        const Vec3 face_velocity =
            weight * velocity[owner] + (1.0 - weight) * velocity[neighbour];
        const Vec3 interpolated_pressure_response =
            weight * mobility[owner] * pressure_gradient[owner] +
            (1.0 - weight) * mobility[neighbour] * pressure_gradient[neighbour];
        face_flux[face] =
            dot(face_velocity, mesh.face_area_vectors[f]) +
            dot(interpolated_pressure_response, mesh.face_area_vectors[f]) -
            face_mobility * mesh.face_orthogonal_coefficients[f] *
                (pressure[neighbour] - pressure[owner]);
    }
}

}  // babelsim 命名空间
