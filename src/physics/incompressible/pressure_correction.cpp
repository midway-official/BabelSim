#include "babelsim/incompressible_operators.h"

#include "babelsim/fvc.h"
#include "babelsim/fvm.h"

#include <stdexcept>

namespace babelsim {
namespace {

void requireCell(const Field<double>& field, const Mesh& mesh, const char* label) {
    field.validateStorage();
    if (&field.mesh() != &mesh || field.location() != FieldLocation::Cell) {
        throw std::invalid_argument(std::string(label) + " must be a cell scalar field");
    }
}

void requireCell(const Field<Vec3>& field, const Mesh& mesh, const char* label) {
    field.validateStorage();
    if (&field.mesh() != &mesh || field.location() != FieldLocation::Cell) {
        throw std::invalid_argument(std::string(label) + " must be a cell vector field");
    }
}

void requireFace(const ScalarField& field, const Mesh& mesh, const char* label) {
    field.validateStorage();
    if (&field.mesh() != &mesh || field.location() != FieldLocation::Face) {
        throw std::invalid_argument(std::string(label) + " must be a face scalar field");
    }
}

}  // 匿名命名空间

SolveResult PressureCorrection::solve(
    RunTime& run_time,
    ScalarField& correction,
    const Mesh& mesh,
    const ScalarField& face_flux,
    const ScalarField& mobility,
    const ScalarField& pressure,
    ScalarField& face_mobility,
    ScalarField& divergence,
    bool has_fixed_pressure,
    bool reset_correction)
{
    if (&run_time.mesh() != &mesh) {
        throw std::invalid_argument("pressure correction mesh does not belong to run time");
    }
    requireCell(correction, mesh, "pressure correction");
    requireFace(face_flux, mesh, "pressure flux");
    requireCell(mobility, mesh, "pressure mobility");
    requireCell(pressure, mesh, "pressure");
    requireFace(face_mobility, mesh, "face mobility");
    requireCell(divergence, mesh, "pressure divergence");

    // 只有 SIMPLE 特有的“以 rAU 构造压力方程”保留在本组件；面插值、
    // 散度、非正交修正和分布式线性求解均由通用 fvc/fvm/RunTime 完成。
    run_time.evaluate(fvc::interpolate(mobility), face_mobility);
    run_time.evaluate(fvc::div(face_flux), divergence);
    if (reset_correction) correction.fill(0.0);
    ScalarEquationControl control;
    control.fix_reference = !has_fixed_pressure;
    return run_time.solve(
        -fvm::laplacian(face_mobility, correction) == -fvm::source(divergence),
        control);
}

void PressureCorrection::apply(
    RunTime& run_time,
    const Mesh& mesh,
    double pressure_relaxation,
    ScalarField& pressure,
    VectorField& velocity,
    ScalarField& face_flux,
    const ScalarField& correction,
    const ScalarField& mobility,
    VectorField& correction_gradient,
    DiffusionMethod diffusion_method)
{
    if (&run_time.mesh() != &mesh ||
        !(pressure_relaxation > 0.0 && pressure_relaxation <= 1.0)) {
        throw std::invalid_argument("pressure correction inputs are invalid");
    }
    requireCell(pressure, mesh, "pressure");
    requireCell(velocity, mesh, "velocity");
    requireFace(face_flux, mesh, "face flux");
    requireCell(correction, mesh, "pressure correction");
    requireCell(mobility, mesh, "pressure mobility");
    requireCell(correction_gradient, mesh, "pressure-correction gradient");

    run_time.evaluate(fvc::grad(correction), correction_gradient);
    for (Index cell : mesh.owned_cells) {
        pressure[cell] += pressure_relaxation * correction[cell];
        velocity[cell] -= mobility[cell] * correction_gradient[cell];
    }
    for (Index face : mesh.owned_faces) {
        const std::size_t index = static_cast<std::size_t>(face);
        const Index owner = mesh.face_owner[index];
        const Index neighbour = mesh.face_neighbour[index];
        if (neighbour != invalid_index) {
            const double owner_weight = mesh.face_owner_weights[index];
            const double face_mobility = owner_weight * mobility[owner] +
                (1.0 - owner_weight) * mobility[neighbour];
            face_flux[face] -= face_mobility * fvc::integratedNormalGradient(
                correction, correction_gradient, face, diffusion_method);
        } else if (pressure.boundary(mesh.face_patch[index]).type ==
                   BoundaryType::FixedValue) {
            face_flux[face] -= mobility[owner] * fvc::integratedNormalGradient(
                correction, correction_gradient, face, diffusion_method);
        }
    }
    run_time.synchronize(pressure);
    run_time.synchronize(velocity);
    run_time.synchronize(face_flux);
}

}  // babelsim 命名空间
