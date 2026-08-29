#include "babelsim/operators.h"

#include "test_util.h"

#include <Eigen/Core>

#include <algorithm>
#include <iostream>

using namespace babelsim;

namespace {

Mesh affineSkewedMesh() {
    const std::array<Index, 3> dimensions{{5, 5, 5}};
    std::vector<Vec3> points;
    for (Index k = 0; k <= dimensions[2]; ++k) {
        for (Index j = 0; j <= dimensions[1]; ++j) {
            for (Index i = 0; i <= dimensions[0]; ++i) {
                points.push_back({
                    static_cast<double>(i) + 0.25 * j + 0.10 * k,
                    static_cast<double>(j) + 0.15 * k,
                    static_cast<double>(k),
                });
            }
        }
    }
    return Mesh::structured(dimensions, std::move(points));
}

double linearValue(const Vec3& point) {
    return 2.0 * point.x - 3.0 * point.y + 0.5 * point.z + 1.0;
}

Vec3 linearVectorValue(const Vec3& point) {
    return {
        point.x + 2.0 * point.y - 0.5 * point.z,
        -point.x + 3.0 * point.z,
        0.25 * point.x - point.y + 2.0 * point.z,
    };
}

}  // 匿名命名空间

int main() {
    const Mesh mesh = Mesh::cartesian({3, 3, 3}, {0, 0, 0}, {3, 3, 3});
    ScalarField linear(mesh, FieldLocation::Cell, "linear");
    VectorField velocity(mesh, FieldLocation::Cell, "U");
    ScalarField quadratic(mesh, FieldLocation::Cell, "quadratic");
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        const Vec3& centre = mesh.cell_centres[static_cast<std::size_t>(cell)];
        linear[cell] = linearValue(centre);
        velocity[cell] = centre;
        quadratic[cell] = squaredNorm(centre);
    }

    const Index centre_cell = mesh.cellId(1, 1, 1);
    VectorField gg(mesh, FieldLocation::Cell, "gradGG");
    VectorField ls(mesh, FieldLocation::Cell, "gradLS");
    gradient(linear, gg, GradientMethod::GreenGauss);
    gradient(linear, ls, GradientMethod::LeastSquares);
    require(near(gg[centre_cell], {2.0, -3.0, 0.5}), "Green-Gauss gradient failed");
    require(near(ls[centre_cell], {2.0, -3.0, 0.5}), "least-squares gradient failed");

    ScalarField div(mesh, FieldLocation::Cell, "divU");
    divergence(velocity, div);
    require(near(div[centre_cell], 3.0), "vector divergence failed");
    ScalarField face_flux(mesh, FieldLocation::Face, "phi");
    ScalarField div_flux(mesh, FieldLocation::Cell, "divPhi");
    flux(velocity, face_flux);
    divergence(face_flux, div_flux);
    require(near(div_flux[centre_cell], 3.0), "flux divergence failed");

    ScalarField lap(mesh, FieldLocation::Cell, "laplacian");
    laplacian(
        quadratic, lap, GradientMethod::GreenGauss,
        DiffusionMethod::Orthogonal);
    require(near(lap[centre_cell], 6.0), "orthogonal Laplacian failed");

    const Mesh skewed = affineSkewedMesh();
    ScalarField skewed_linear(skewed, FieldLocation::Cell, "linear");
    for (Index cell = 0; cell < skewed.cellCount(); ++cell) {
        skewed_linear[cell] = linearValue(
            skewed.cell_centres[static_cast<std::size_t>(cell)]);
    }
    const Index skewed_centre = skewed.cellId(2, 2, 2);
    VectorField skewed_gradient(skewed, FieldLocation::Cell, "grad");
    gradient(skewed_linear, skewed_gradient, GradientMethod::LeastSquares);
    require(
        near(skewed_gradient[skewed_centre], {2.0, -3.0, 0.5}, 1e-10),
        "least-squares gradient is not exact on an affine skew mesh");
    ScalarField corrected(skewed, FieldLocation::Cell, "correctedLap");
    laplacian(
        skewed_linear, corrected, GradientMethod::LeastSquares,
        DiffusionMethod::Corrected);
    require(
        std::abs(corrected[skewed_centre]) < 1e-10,
        "non-orthogonal correction does not preserve a linear field");
    ScalarField skewed_quadratic(skewed, FieldLocation::Cell, "quadratic");
    for (Index cell = 0; cell < skewed.cellCount(); ++cell) {
        const Vec3& point =
            skewed.cell_centres[static_cast<std::size_t>(cell)];
        skewed_quadratic[cell] = point.x * point.y;
    }
    ScalarField orthogonal_quadratic(
        skewed, FieldLocation::Cell, "orthogonalQuadraticLap");
    ScalarField corrected_quadratic(
        skewed, FieldLocation::Cell, "correctedQuadraticLap");
    laplacian(
        skewed_quadratic, orthogonal_quadratic,
        GradientMethod::LeastSquares, DiffusionMethod::Orthogonal);
    laplacian(
        skewed_quadratic, corrected_quadratic,
        GradientMethod::LeastSquares, DiffusionMethod::Corrected);
    require(
        near(corrected_quadratic[skewed_centre], 0.0, 1e-10),
        "corrected skew-mesh Laplacian failed the quadratic manufactured solution");
    require(
        std::abs(orthogonal_quadratic[skewed_centre]) > 1e-3,
        "manufactured solution does not exercise non-orthogonal correction");

    VectorField skewed_vector(skewed, FieldLocation::Cell, "linearVector");
    for (Index cell = 0; cell < skewed.cellCount(); ++cell) {
        skewed_vector[cell] = linearVectorValue(
            skewed.cell_centres[static_cast<std::size_t>(cell)]);
    }
    TensorField skewed_tensor(skewed, FieldLocation::Cell, "gradVector");
    gradient(skewed_vector, skewed_tensor, GradientMethod::LeastSquares);
    require(
        near(skewed_tensor[skewed_centre][0], {1.0, 2.0, -0.5}, 1e-10) &&
            near(skewed_tensor[skewed_centre][1], {-1.0, 0.0, 3.0}, 1e-10) &&
            near(skewed_tensor[skewed_centre][2], {0.25, -1.0, 2.0}, 1e-10),
        "least-squares vector gradient is not exact on an affine skew mesh");
    VectorEquation vector_diffusion(skewed);
    addDiffusion(
        vector_diffusion, 1.0, skewed_vector,
        GradientMethod::LeastSquares, DiffusionMethod::Corrected);
    Vec3 vector_residual =
        vector_diffusion.diagonal[static_cast<std::size_t>(skewed_centre)] *
            skewed_vector[skewed_centre] -
        vector_diffusion.source[static_cast<std::size_t>(skewed_centre)];
    for (Index face : skewed.cell_faces[static_cast<std::size_t>(skewed_centre)]) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = skewed.face_owner[f];
        const Index neighbour = skewed.face_neighbour[f];
        if (neighbour == invalid_index) {
            continue;
        }
        vector_residual += owner == skewed_centre
            ? vector_diffusion.upper[f] * skewed_vector[neighbour]
            : vector_diffusion.lower[f] * skewed_vector[owner];
    }
    require(
        norm(vector_residual) < 1e-10,
        "corrected vector diffusion does not preserve an affine linear field");

    ScalarField reconstructed_scalar(skewed, FieldLocation::Face, "linearFace");
    VectorField reconstructed_vector(skewed, FieldLocation::Face, "vectorFace");
    reconstruct(skewed_linear, skewed_gradient, reconstructed_scalar);
    reconstruct(skewed_vector, skewed_tensor, reconstructed_vector);
    for (Index face : skewed.cell_faces[static_cast<std::size_t>(skewed_centre)]) {
        const auto f = static_cast<std::size_t>(face);
        if (skewed.face_neighbour[f] == invalid_index) {
            continue;
        }
        require(
            near(reconstructed_scalar[face], linearValue(skewed.face_centres[f]), 1e-10),
            "skew-corrected scalar reconstruction is not affine exact");
        require(
            near(reconstructed_vector[face], linearVectorValue(skewed.face_centres[f]), 1e-10),
            "skew-corrected vector reconstruction is not affine exact");
    }

    ScalarEquation time_equation(mesh);
    addTimeDerivative(
        time_equation, linear, 0.25, 2.0, TimeMethod::BDF2, &quadratic);
    const double coefficient = 2.0 * mesh.cell_volumes[
        static_cast<std::size_t>(centre_cell)] / 0.25;
    require(
        near(time_equation.diagonal[static_cast<std::size_t>(centre_cell)],
             1.5 * coefficient),
        "BDF2 diagonal is incorrect");
    VectorField older_velocity(
        mesh, FieldLocation::Cell, "olderU", {0.5, -0.25, 0.75});
    VectorEquation vector_time_equation(mesh);
    addTimeDerivative(
        vector_time_equation, velocity, 0.25, 2.0,
        TimeMethod::BDF2, &older_velocity);
    require(
        near(
            vector_time_equation.source[static_cast<std::size_t>(centre_cell)],
            coefficient *
                (2.0 * velocity[centre_cell] - 0.5 * older_velocity[centre_cell])),
        "vector BDF2 source is incorrect");

    VectorField uniform_velocity(mesh, FieldLocation::Cell, "uniformU", {1.0, 0.0, 0.0});
    ScalarField uniform_flux(mesh, FieldLocation::Face, "uniformPhi");
    flux(uniform_velocity, uniform_flux);
    ScalarField constant(mesh, FieldLocation::Cell, "constant", 2.0);
    ScalarEquation convection(mesh);
    addConvection(convection, uniform_flux, constant, ConvectionMethod::Upwind);
    double maximum_constant_residual = 0.0;
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        double row_value = convection.diagonal[static_cast<std::size_t>(cell)] * 2.0;
        for (Index face : mesh.cell_faces[static_cast<std::size_t>(cell)]) {
            const auto f = static_cast<std::size_t>(face);
            if (mesh.face_neighbour[f] == invalid_index) {
                continue;
            }
            row_value += (mesh.face_owner[f] == cell
                ? convection.upper[f]
                : convection.lower[f]) * 2.0;
        }
        row_value -= convection.source[static_cast<std::size_t>(cell)];
        maximum_constant_residual = std::max(
            maximum_constant_residual, std::abs(row_value));
    }
    require(
        maximum_constant_residual < 1e-11,
        "upwind convection does not preserve a constant field");
    ScalarEquation central_convection(mesh);
    addConvection(
        central_convection, uniform_flux, constant,
        ConvectionMethod::Central);
    maximum_constant_residual = 0.0;
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        double row_value =
            central_convection.diagonal[static_cast<std::size_t>(cell)] * 2.0;
        for (Index face : mesh.cell_faces[static_cast<std::size_t>(cell)]) {
            const auto f = static_cast<std::size_t>(face);
            if (mesh.face_neighbour[f] == invalid_index) {
                continue;
            }
            row_value += (mesh.face_owner[f] == cell
                ? central_convection.upper[f]
                : central_convection.lower[f]) * 2.0;
        }
        row_value -= central_convection.source[static_cast<std::size_t>(cell)];
        maximum_constant_residual = std::max(
            maximum_constant_residual, std::abs(row_value));
    }
    require(
        maximum_constant_residual < 1e-11,
        "central convection does not preserve a constant field");

    std::cout << "operators_test: scalar/vector grad, corrected diffusion, "
                 "skew reconstruction, div/laplacian/time/convection passed\n";
}
