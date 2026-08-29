#include "babelsim/incompressible_operators.h"

#include "babelsim/operators.h"

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
    mesh.validate();
    apply(
        mesh, velocity, pressure, mobility, pressure_gradient, face_flux,
        InterpolationMethod::Linear, GradientMethod::GreenGauss,
        DiffusionMethod::Orthogonal);
}

void MomentumInterpolation::apply(
    const Mesh& mesh,
    const VectorField& velocity,
    const ScalarField& pressure,
    const ScalarField& mobility,
    const VectorField& pressure_gradient,
    ScalarField& face_flux,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
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
    velocity.validateStorage();
    pressure.validateStorage();
    mobility.validateStorage();
    pressure_gradient.validateStorage();
    face_flux.validateStorage();

    // H/a 在真实面中心重构；压力响应仍保留为独立项，使直接面压力差能够抑制
    // 同位网格的棋盘格压力。
    flux(velocity, face_flux, interpolation_method, gradient_method);
    VectorField pressure_response(mesh, FieldLocation::Cell, "rAUGradP");
    VectorField face_pressure_response(mesh, FieldLocation::Face, "rAUGradPFace");
    ScalarField face_mobility(mesh, FieldLocation::Face, "rAUFace");
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        pressure_response[cell] = mobility[cell] * pressure_gradient[cell];
    }
    interpolate(
        pressure_response, face_pressure_response,
        interpolation_method, gradient_method);
    interpolate(mobility, face_mobility, interpolation_method, gradient_method);

    for (Index face : mesh.owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index neighbour = mesh.face_neighbour[f];
        if (neighbour == invalid_index) {
            continue;
        }
        face_flux[face] +=
            dot(face_pressure_response[face], mesh.face_area_vectors[f]) -
            face_mobility[face] * integratedNormalGradient(
                pressure, pressure_gradient, face, diffusion_method);
    }
}

}  // babelsim 命名空间
