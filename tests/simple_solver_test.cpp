#include "babelsim/incompressible.h"
#include "babelsim/legacy_taiho.h"

#include "test_util.h"

#include <cmath>
#include <iostream>

using namespace babelsim;

int main() {
    ImportedTaihoMesh imported = readTaihoMesh("tests/data/taiho_5x5");
    IncompressibleFields fields(imported.mesh);
    applyImportedBoundaryConditions(imported, fields.velocity, fields.pressure);

    Methods methods;
    methods.gradient = GradientMethod::GreenGauss;
    methods.convection = ConvectionMethod::Upwind;
    methods.diffusion = DiffusionMethod::Orthogonal;
    SimpleControl control;
    control.max_iterations = 500;
    control.velocity_relaxation = 0.5;
    control.pressure_relaxation = 0.3;
    control.continuity_tolerance = 1e-9;
    control.velocity_tolerance = 1e-7;
    control.velocity_solver.absolute_tolerance = 1e-13;
    control.velocity_solver.relative_tolerance = 1e-10;
    control.pressure_solver.absolute_tolerance = 1e-13;
    control.pressure_solver.relative_tolerance = 1e-10;

    SimpleSolver solver(fields, {1.0, 0.1}, methods, control);
    SimpleIterationResult result;
    int iterations = 0;
    for (int iteration = 1; iteration <= control.max_iterations; ++iteration) {
        result = solver.iterate();
        iterations = iteration;
        require(result.healthy, "SIMPLE produced a numerical failure");
        if (result.converged) {
            break;
        }
    }
    require(result.converged, "SIMPLE did not converge on the imported channel");
    require(
        result.continuity.relative <= control.continuity_tolerance,
        "SIMPLE converged without satisfying continuity");

    double maximum_z_velocity = 0.0;
    for (const Vec3& velocity : fields.velocity.values()) {
        maximum_z_velocity = std::max(maximum_z_velocity, std::abs(velocity.z));
        require(isFinite(velocity), "SIMPLE velocity contains a non-finite value");
    }
    require(maximum_z_velocity < 1e-13, "nz=1 SIMPLE generated a z velocity");

    IncompressibleFields transient_fields(imported.mesh);
    applyImportedBoundaryConditions(
        imported, transient_fields.velocity, transient_fields.pressure);
    VectorField previous(imported.mesh, FieldLocation::Cell, "Uold");
    Methods transient_methods = methods;
    transient_methods.time = TimeMethod::Euler;
    SimpleSolver transient_solver(
        transient_fields, {1.0, 0.1}, transient_methods, control);
    const SimpleIterationResult transient = transient_solver.iterate(
        {0.1, &previous, nullptr});
    require(
        transient.healthy,
        "SIMPLE ignored or failed the configured Euler time method");

    std::cout << "simple_solver_test: iterations=" << iterations
              << " mass=" << result.continuity.relative
              << " dU=" << result.relative_velocity_change << '\n';
}
