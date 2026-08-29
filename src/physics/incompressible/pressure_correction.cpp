#include "babelsim/incompressible_operators.h"
#include "babelsim/operators.h"

#include <stdexcept>

namespace babelsim {

void PressureCorrection::assemble(
    ScalarEquation& equation,
    ScalarField& correction,
    const Mesh& mesh,
    const ScalarField& face_flux,
    const ScalarField& mobility,
    const ScalarField& pressure,
    bool has_fixed_pressure,
    const ParallelContext& parallel)
{
    mesh.validate();
    assemble(
        equation, correction, mesh, face_flux, mobility, pressure, nullptr,
        DiffusionMethod::Orthogonal, has_fixed_pressure, parallel);
}

void PressureCorrection::assemble(
    ScalarEquation& equation,
    ScalarField& correction,
    const Mesh& mesh,
    const ScalarField& face_flux,
    const ScalarField& mobility,
    const ScalarField& pressure,
    const VectorField* correction_gradient,
    DiffusionMethod diffusion_method,
    bool has_fixed_pressure,
    const ParallelContext& parallel)
{
    if (equation.mesh != &mesh || &correction.mesh() != &mesh ||
        &face_flux.mesh() != &mesh || &mobility.mesh() != &mesh ||
        &pressure.mesh() != &mesh || correction.location() != FieldLocation::Cell ||
        face_flux.location() != FieldLocation::Face ||
        mobility.location() != FieldLocation::Cell ||
        pressure.location() != FieldLocation::Cell ||
        (diffusion_method != DiffusionMethod::Orthogonal &&
         (correction_gradient == nullptr ||
          &correction_gradient->mesh() != &mesh ||
          correction_gradient->location() != FieldLocation::Cell))) {
        throw std::invalid_argument("pressure-correction fields do not match mesh");
    }
    correction.validateStorage();
    face_flux.validateStorage();
    mobility.validateStorage();
    pressure.validateStorage();
    if (correction_gradient != nullptr) correction_gradient->validateStorage();
    equation.reset();
    for (Index face : mesh.owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        equation.source[static_cast<std::size_t>(owner)] -= face_flux[face];
        if (neighbour != invalid_index) {
            equation.source[static_cast<std::size_t>(neighbour)] += face_flux[face];
            const double weight = mesh.face_owner_weights[f];
            const double face_mobility =
                weight * mobility[owner] + (1.0 - weight) * mobility[neighbour];
            const double coefficient =
                face_mobility * mesh.face_orthogonal_coefficients[f];
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.diagonal[static_cast<std::size_t>(neighbour)] += coefficient;
            equation.upper[f] -= coefficient;
            equation.lower[f] -= coefficient;
            if (correction_gradient != nullptr) {
                const double orthogonal_part =
                    mesh.face_orthogonal_coefficients[f] *
                    (correction[neighbour] - correction[owner]);
                const double explicit_correction = face_mobility *
                    (integratedNormalGradient(
                         correction, *correction_gradient, face,
                         diffusion_method) - orthogonal_part);
                equation.source[static_cast<std::size_t>(owner)] +=
                    explicit_correction;
                equation.source[static_cast<std::size_t>(neighbour)] -=
                    explicit_correction;
            }
        } else if (pressure.boundary(mesh.face_patch[f]).type ==
                   BoundaryType::FixedValue) {
            const double coefficient =
                mobility[owner] * mesh.face_orthogonal_coefficients[f];
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            if (correction_gradient != nullptr) {
                const double orthogonal_part =
                    mesh.face_orthogonal_coefficients[f] * (-correction[owner]);
                equation.source[static_cast<std::size_t>(owner)] += mobility[owner] *
                    (integratedNormalGradient(
                         correction, *correction_gradient, face,
                         diffusion_method) - orthogonal_part);
            }
        }
    }
    if (!has_fixed_pressure && parallel.rank == 0) {
        if (mesh.owned_cells.empty()) {
            throw std::runtime_error("closed domain has no pressure reference cell");
        }
        const auto reference = static_cast<std::size_t>(mesh.owned_cells.front());
        if (!(equation.diagonal[reference] > 0.0)) {
            throw std::runtime_error("pressure reference diagonal is invalid");
        }
        equation.diagonal[reference] += equation.diagonal[reference];
    }
}

void PressureCorrection::apply(
    const Mesh& mesh,
    double pressure_relaxation,
    ScalarField& pressure,
    VectorField& velocity,
    ScalarField& face_flux,
    const ScalarField& correction,
    const VectorField& correction_gradient,
    const ScalarField& mobility)
{
    apply(
        mesh, pressure_relaxation, pressure, velocity, face_flux, correction,
        correction_gradient, mobility, DiffusionMethod::Orthogonal);
}

void PressureCorrection::apply(
    const Mesh& mesh,
    double pressure_relaxation,
    ScalarField& pressure,
    VectorField& velocity,
    ScalarField& face_flux,
    const ScalarField& correction,
    const VectorField& correction_gradient,
    const ScalarField& mobility,
    DiffusionMethod diffusion_method)
{
    if (!(pressure_relaxation > 0.0 && pressure_relaxation <= 1.0) ||
        &pressure.mesh() != &mesh || &velocity.mesh() != &mesh ||
        &face_flux.mesh() != &mesh || &correction.mesh() != &mesh ||
        &correction_gradient.mesh() != &mesh || &mobility.mesh() != &mesh) {
        throw std::invalid_argument("pressure-correction inputs are invalid");
    }
    for (Index cell : mesh.owned_cells) {
        pressure[cell] += pressure_relaxation * correction[cell];
        velocity[cell] -= mobility[cell] * correction_gradient[cell];
    }
    for (Index face : mesh.owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        if (neighbour != invalid_index) {
            const double weight = mesh.face_owner_weights[f];
            const double face_mobility =
                weight * mobility[owner] + (1.0 - weight) * mobility[neighbour];
            face_flux[face] -= face_mobility * integratedNormalGradient(
                correction, correction_gradient, face, diffusion_method);
        } else if (pressure.boundary(mesh.face_patch[f]).type ==
                   BoundaryType::FixedValue) {
            face_flux[face] -= mobility[owner] * integratedNormalGradient(
                correction, correction_gradient, face, diffusion_method);
        }
    }
}

}  // babelsim 命名空间
