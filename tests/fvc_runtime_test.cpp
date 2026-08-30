#include "babelsim/runtime.h"

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
    scalar[mesh.cellId(0, 0, 0)] = 1.0;
    scalar[mesh.cellId(1, 0, 0)] = 3.0;
    for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        scalar.boundary(patch) = zeroGradient();
    }
    VectorField scalar_gradient(mesh, FieldLocation::Cell, "gradT");
    ScalarField face_scalar(mesh, FieldLocation::Face, "Tf");
    ScalarField reconstructed_scalar(mesh, FieldLocation::Face, "Trec");
    ScalarField unit_laplacian(mesh, FieldLocation::Cell, "lapT");
    ScalarField doubled_laplacian(mesh, FieldLocation::Cell, "lap2T");
    ScalarField cell_diffusivity(mesh, FieldLocation::Cell, "k", 2.0);
    ScalarField field_laplacian(mesh, FieldLocation::Cell, "lapkT");
    run_time.evaluate(fvc::grad(scalar), scalar_gradient);
    run_time.evaluate(fvc::interpolate(scalar), face_scalar);
    run_time.evaluate(fvc::reconstruct(scalar, scalar_gradient), reconstructed_scalar);
    run_time.evaluate(fvc::laplacian(scalar), unit_laplacian);
    run_time.evaluate(fvc::laplacian(2.0, scalar), doubled_laplacian);
    run_time.evaluate(fvc::laplacian(cell_diffusivity, scalar), field_laplacian);
    for (Index cell : mesh.owned_cells) {
        require(
            near(doubled_laplacian[cell], 2.0 * unit_laplacian[cell]),
            "constant explicit diffusivity is inconsistent");
        require(
            near(field_laplacian[cell], doubled_laplacian[cell]),
            "cell explicit diffusivity was not interpolated consistently");
    }

    VectorField velocity(mesh, FieldLocation::Cell, "U");
    TensorField velocity_gradient(mesh, FieldLocation::Cell, "gradU");
    VectorField reconstructed_velocity(mesh, FieldLocation::Face, "Urec");
    ScalarField flux(mesh, FieldLocation::Face, "phi");
    ScalarField flux_divergence(mesh, FieldLocation::Cell, "divPhi");
    ScalarField vector_divergence(mesh, FieldLocation::Cell, "divU");
    ScalarField scalar_convection(mesh, FieldLocation::Cell, "divPhiT");
    VectorField vector_convection(mesh, FieldLocation::Cell, "divPhiU");
    run_time.evaluate(fvc::flux(velocity), flux);
    run_time.evaluate(fvc::grad(velocity), velocity_gradient);
    run_time.evaluate(
        fvc::reconstruct(velocity, velocity_gradient), reconstructed_velocity);
    run_time.evaluate(fvc::div(flux), flux_divergence);
    run_time.evaluate(fvc::div(velocity), vector_divergence);
    run_time.evaluate(fvc::div(flux, scalar), scalar_convection);
    run_time.evaluate(fvc::div(flux, velocity), vector_convection);
    for (Index cell : mesh.owned_cells) {
        require(near(flux_divergence[cell], 0.0), "zero face flux has divergence");
        require(near(vector_divergence[cell], 0.0), "zero velocity has divergence");
        require(near(scalar_convection[cell], 0.0), "zero flux has scalar convection");
        require(near(vector_convection[cell], {}), "zero flux has vector convection");
    }
    std::cout << "fvc_runtime_test: explicit Field operations passed\n";
}
