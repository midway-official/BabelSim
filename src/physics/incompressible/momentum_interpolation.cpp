#include "babelsim/incompressible_operators.h"

#include <stdexcept>

namespace babelsim {

void MomentumInterpolation::apply(
    RunTime& run_time,
    const Mesh& mesh,
    const VectorField& velocity,
    const ScalarField& pressure,
    const ScalarField& mobility,
    const VectorField& pressure_gradient,
    ScalarField& face_flux,
    MomentumInterpolationWorkspace& workspace)
{
    mesh.validate();
    apply(
        run_time, mesh, velocity, pressure, mobility, pressure_gradient, face_flux, workspace,
        InterpolationMethod::Linear, GradientMethod::GreenGauss,
        DiffusionMethod::Orthogonal);
}

void MomentumInterpolation::apply(
    RunTime& run_time,
    const Mesh& mesh,
    const VectorField& velocity,
    const ScalarField& pressure,
    const ScalarField& mobility,
    const VectorField& pressure_gradient,
    ScalarField& face_flux,
    MomentumInterpolationWorkspace& workspace,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    if (&velocity.mesh() != &mesh || &pressure.mesh() != &mesh ||
        &mobility.mesh() != &mesh || &pressure_gradient.mesh() != &mesh ||
        &face_flux.mesh() != &mesh || &workspace.pressure_response.mesh() != &mesh ||
        &workspace.face_pressure_response.mesh() != &mesh ||
        &workspace.face_mobility.mesh() != &mesh ||
        velocity.location() != FieldLocation::Cell ||
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
    workspace.pressure_response.validateStorage();
    workspace.face_pressure_response.validateStorage();
    workspace.face_mobility.validateStorage();
    if (&run_time.mesh() != &mesh) {
        throw std::invalid_argument("momentum interpolation mesh does not belong to run time");
    }

    // H/a 在真实面中心重构；压力响应仍保留为独立项，使直接面压力差能够抑制
    // 同位网格的棋盘格压力。
    if (interpolation_method != run_time.methods().interpolation ||
        gradient_method != run_time.methods().gradient) {
        throw std::invalid_argument("momentum interpolation methods must match the run time");
    }
    run_time.evaluate(fvc::flux(velocity), face_flux);
    run_time.synchronize(const_cast<ScalarField&>(mobility));
    run_time.synchronize(const_cast<VectorField&>(pressure_gradient));
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        workspace.pressure_response[cell] = mobility[cell] * pressure_gradient[cell];
    }
    run_time.evaluate(
        fvc::interpolate(workspace.pressure_response), workspace.face_pressure_response);
    run_time.evaluate(fvc::interpolate(mobility), workspace.face_mobility);

    for (Index face : mesh.owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index neighbour = mesh.face_neighbour[f];
        if (neighbour == invalid_index) {
            continue;
        }
        face_flux[face] +=
            dot(workspace.face_pressure_response[face], mesh.face_area_vectors[f]) -
            workspace.face_mobility[face] * fvc::integratedNormalGradient(
                pressure, pressure_gradient, face, diffusion_method);
    }
    run_time.synchronize(face_flux);
}

}  // babelsim 命名空间
