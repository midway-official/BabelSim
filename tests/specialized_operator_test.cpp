#include "babelsim/incompressible_operators.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    const Mesh mesh = Mesh::cartesian({2, 2, 1}, {0, 0, 0}, {1, 1, 1});
    VectorField velocity(mesh, FieldLocation::Cell, "U");
    ScalarField pressure(mesh, FieldLocation::Cell, "p");
    ScalarField mobility(mesh, FieldLocation::Cell, "rAU", 0.25);
    VectorField gradient(mesh, FieldLocation::Cell, "gradP");
    ScalarField flux(mesh, FieldLocation::Face, "phi");
    for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        velocity.setBoundary(patch, BoundaryCondition<Vec3>::zeroGradient());
        pressure.setBoundary(patch, BoundaryCondition<double>::zeroGradient());
    }

    MomentumInterpolation::apply(mesh, velocity, pressure, mobility, gradient, flux);
    for (double value : flux.values()) require(near(value, 0.0), "zero state gained flux");

    ScalarEquation equation(mesh);
    ScalarField correction(mesh, FieldLocation::Cell, "pPrime");
    PressureCorrection::assemble(
        equation, correction, mesh, flux, mobility, pressure, false, ParallelContext{});
    for (Index cell : mesh.owned_cells) {
        require(equation.diagonal[static_cast<std::size_t>(cell)] > 0.0,
                "pressure correction has no positive diagonal");
    }
    VectorField correction_gradient(mesh, FieldLocation::Cell, "gradPPrime");
    PressureCorrection::apply(
        mesh, 0.3, pressure, velocity, flux, correction, correction_gradient, mobility);
    for (Index cell : mesh.owned_cells) {
        require(near(pressure[cell], 0.0) && near(velocity[cell], {}),
                "zero pressure correction changed the field");
    }
    std::cout << "specialized_operator_test: momentum interpolation and pressure correction passed\n";
}
