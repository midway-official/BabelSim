#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/runtime.h"
#include "babelsim/operators.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    const Mesh mesh = Mesh::cartesian({2, 1, 1}, {0, 0, 0}, {2, 1, 1});
    RuntimeControl control;
    control.methods.interpolation = InterpolationMethod::Linear;
    control.methods.gradient = GradientMethod::GreenGauss;
    control.methods.convection = ConvectionMethod::Upwind;
    control.methods.diffusion = DiffusionMethod::Orthogonal;
    RunTime run_time = RunTime::forMesh(mesh, control);

    ScalarField scalar(mesh, FieldLocation::Cell, "T");
    detail::fieldData(scalar)[mesh.cellId(0, 0, 0)] = 1.0;
    detail::fieldData(scalar)[mesh.cellId(1, 0, 0)] = 3.0;
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        scalar.boundary(patch) = zeroGradient();
    }
    VectorField scalar_gradient(mesh, FieldLocation::Cell, "gradT");
    ScalarField face_scalar(mesh, FieldLocation::Face, "Tf");
    ScalarField reconstructed_scalar(mesh, FieldLocation::Face, "Trec");
    ScalarField unit_laplacian(mesh, FieldLocation::Cell, "lapT");
    ScalarField doubled_laplacian(mesh, FieldLocation::Cell, "lap2T");
    ScalarField cell_diffusivity(mesh, FieldLocation::Cell, "k", 2.0);
    ScalarField field_laplacian(mesh, FieldLocation::Cell, "lapkT");
    math::evaluate(math::grad(scalar), scalar_gradient);
    math::evaluate(math::interpolate(scalar), face_scalar);
    math::evaluate(math::reconstruct(scalar, scalar_gradient), reconstructed_scalar);
    math::evaluate(math::laplacian(scalar), unit_laplacian);
    math::evaluate(math::laplacian(2.0, scalar), doubled_laplacian);
    math::evaluate(math::laplacian(cell_diffusivity, scalar), field_laplacian);
    for (Index cell : detail::meshData(mesh).owned_cells) {
        require(
            near(detail::fieldData(doubled_laplacian)[cell], 2.0 * detail::fieldData(unit_laplacian)[cell]),
            "constant explicit diffusivity is inconsistent");
        require(
            near(detail::fieldData(field_laplacian)[cell], detail::fieldData(doubled_laplacian)[cell]),
            "cell explicit diffusivity was not interpolated consistently");
    }

    VectorField velocity_correction(mesh, FieldLocation::Cell, "Ucorrected", {1.0, 1.0, 1.0});
    math::subtract(cell_diffusivity, math::grad(scalar), velocity_correction);
    for (Index cell : detail::meshData(mesh).owned_cells) {
        require(
            near(
                detail::fieldData(velocity_correction)[cell],
                Vec3{1.0, 1.0, 1.0} - 2.0 * detail::fieldData(scalar_gradient)[cell]),
            "high-level gradient correction is inconsistent");
    }

    ScalarField flux_correction(mesh, FieldLocation::Face, "pFlux", 0.0);
    math::subtract(math::flux(cell_diffusivity, scalar), flux_correction);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const double expected = -2.0 * integratedNormalGradient(
            scalar, scalar_gradient, face, control.methods.diffusion);
        require(
            near(detail::fieldData(flux_correction)[face], expected),
            "high-level diffusion-flux correction is inconsistent");
    }

    // 面矢量直接投影；给定梯度的扩散通量复用重构，不把 SIMPLE 公式放回框架。
    VectorField face_response(mesh, FieldLocation::Face, "response", {2.0, -1.0, 0.5});
    ScalarField face_coefficient(mesh, FieldLocation::Face, "kf", 1.7);
    ScalarField correction(mesh, FieldLocation::Face, "correction", 3.0);
    ScalarField projected(mesh, FieldLocation::Face, "projected");
    math::evaluate(math::flux(face_response), projected);
    math::add(math::flux(face_response), correction, math::FaceRegion::Interior);
    math::subtract(math::flux(face_coefficient, math::reconstruct(scalar, scalar_gradient)),
                  correction, math::FaceRegion::Interior);
    for (Index face : detail::meshData(mesh).owned_faces) {
        const double projection = dot(detail::fieldData(face_response)[face], mesh.faceAreaVector(face));
        require(near(detail::fieldData(projected)[face], projection), "face flux interpolated an already face-centred field");
        double expected = 3.0;
        if (detail::meshData(mesh).face_neighbour[face] != invalid_index)
            expected += projection - 1.7 * integratedNormalGradient(
                scalar, scalar_gradient, face, control.methods.diffusion);
        require(near(detail::fieldData(correction)[face], expected), "interior flux update changed the boundary or its sign");
    }
    math::add(math::flux(face_response), projected);
    for (Index face : detail::meshData(mesh).owned_faces)
        require(near(detail::fieldData(projected)[face], 2 * dot(
            detail::fieldData(face_response)[face], mesh.faceAreaVector(face))), "all-face addition omitted a boundary");

    const auto rejects = [](const auto& operation) {
        try { operation(); } catch (const std::invalid_argument&) { return; }
        throw std::runtime_error("invalid public flux operation was accepted");
    };
    rejects([&] { math::add(math::flux(face_response), scalar); });
    rejects([&] { math::subtract(math::flux(face_coefficient, scalar), face_coefficient); });
    rejects([&] { math::evaluate(math::flux(face_coefficient, math::reconstruct(scalar, face_response)), projected); });
    rejects([&] { math::add(math::flux(face_response), projected, static_cast<math::FaceRegion>(-1)); });
    const Mesh other_mesh = Mesh::cartesian({2, 1, 1}, {0, 0, 0}, {2, 1, 1});
    VectorField other_face(other_mesh, FieldLocation::Face);
    rejects([&] { math::evaluate(math::flux(other_face), projected); });

    VectorField velocity(mesh, FieldLocation::Cell, "U");
    TensorField velocity_gradient(mesh, FieldLocation::Cell, "gradU");
    VectorField reconstructed_velocity(mesh, FieldLocation::Face, "Urec");
    ScalarField flux(mesh, FieldLocation::Face, "phi");
    ScalarField flux_divergence(mesh, FieldLocation::Cell, "divPhi");
    ScalarField vector_divergence(mesh, FieldLocation::Cell, "divU");
    ScalarField scalar_convection(mesh, FieldLocation::Cell, "divPhiT");
    VectorField vector_convection(mesh, FieldLocation::Cell, "divPhiU");
    math::evaluate(math::flux(velocity), flux);
    math::evaluate(math::grad(velocity), velocity_gradient);
    math::evaluate(
        math::reconstruct(velocity, velocity_gradient), reconstructed_velocity);
    math::evaluate(math::div(flux), flux_divergence);
    math::evaluate(math::div(velocity), vector_divergence);
    math::evaluate(math::div(flux, scalar), scalar_convection);
    math::evaluate(math::div(flux, velocity), vector_convection);
    for (Index cell : detail::meshData(mesh).owned_cells) {
        require(near(detail::fieldData(flux_divergence)[cell], 0.0), "zero face flux has divergence");
        require(near(detail::fieldData(vector_divergence)[cell], 0.0), "zero velocity has divergence");
        require(near(detail::fieldData(scalar_convection)[cell], 0.0), "zero flux has scalar convection");
        require(near(detail::fieldData(vector_convection)[cell], {}), "zero flux has vector convection");
    }
    std::cout << "math_runtime_test: explicit Field operations passed\n";
}
