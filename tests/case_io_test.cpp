#include "babelsim/case.h"
#include "babelsim/numerics_io.h"
#include "babelsim/field_io.h"
#include "babelsim/mesh_io.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    const CaseDefinition cavity = readCase("cases/cavity");
    require(cavity.solver == "simple", "cavity solver selection is incorrect");
    const Parameters physics(cavity.physics_file);
    const Parameters solution(cavity.solution_file);
    require(near(physics.number("density"), 1.0), "cavity density is incorrect");
    require(near(physics.number("dynamicViscosity"), 0.01), "cavity viscosity is incorrect");
    require(solution.integer("maxIterations", 0) == 5000, "cavity iteration limit is incorrect");

    const Mesh mesh = readMeshFile(cavity.mesh_file);
    ScalarField pressure(mesh, FieldLocation::Cell, "p");
    VectorField velocity(mesh, FieldLocation::Cell, "U");
    readFieldFile(cavity.fields_directory / "U.field", velocity);
    readFieldFile(cavity.fields_directory / "p.field", pressure);
    require(velocity.boundary(static_cast<Index>(3)).type == BoundaryType::FixedValue,
            "lid velocity condition was not read");
    require(near(velocity.boundary(static_cast<Index>(3)).value, {1.0, 0.0, 0.0}),
            "lid velocity value was not read");
    require(pressure.boundary(static_cast<Index>(4)).type == BoundaryType::Symmetry,
            "scalar symmetry condition was not read");

    const CaseDefinition channel = readCase("cases/poiseuille");
    const Mesh channel_mesh = readMeshFile(channel.mesh_file);
    VectorField channel_velocity(channel_mesh, FieldLocation::Cell, "U");
    readFieldFile(channel.fields_directory / "U.field", channel_velocity);
    require(channel_velocity.boundary(static_cast<Index>(0)).type == BoundaryType::FixedValue,
            "inlet condition was not read");
    require(near(channel_velocity.boundary(static_cast<Index>(0)).value, {1.0, 0.0, 0.0}),
            "inlet velocity value was not read");

    const CaseDefinition transport = readCase("cases/transport");
    const Parameters transport_physics(transport.physics_file);
    const Methods methods = readMethodsFile(transport.methods_file);
    require(transport.solver == "transport", "transport solver selection is incorrect");
    require(near(transport_physics.number("diffusivity"), 0.01), "transport diffusivity is incorrect");
    require(methods.convectionFor("C") == ConvectionMethod::Upwind,
            "Field-specific convection method was not read");
    const Methods overrides = readMethodsFile("tests/data/methods.bs");
    require(overrides.interpolationFor("T") == InterpolationMethod::Corrected &&
            overrides.gradientFor("T") == GradientMethod::LeastSquares &&
            overrides.convectionFor("T") == ConvectionMethod::LinearUpwind &&
            overrides.diffusionFor("T") == DiffusionMethod::LimitedCorrected,
            "method overrides were mixed between operator types");
    require(overrides.interpolationFor("C") == InterpolationMethod::Linear &&
            overrides.gradientFor("C") == GradientMethod::GreenGauss &&
            overrides.convectionFor("C") == ConvectionMethod::Upwind &&
            overrides.diffusionFor("C") == DiffusionMethod::Orthogonal &&
            overrides.time == TimeMethod::Euler,
        "method overrides changed the default or time method");

    ConfigLine amg_line;
    amg_line.number = 1;
    amg_line.tokens = {
        "scalarSolver", "bicgstab", "amg", "1e-14", "1e-9", "400",
        "amgMaxLevels=9", "amgCoarseSize=24",
        "amgSmoothingSteps=3", "amgRefreshInterval=5"};
    LinearSolverConfig amg_config;
    readLinearSolverLine("tests/data/solution.bs", amg_line, amg_config);
    amg_config.validate();
    require(
        amg_config.solver == LinearSolverType::BiCGSTAB &&
            amg_config.preconditioner == PreconditionerType::AlgebraicMultigrid &&
            amg_config.amg_max_levels == 9 && amg_config.amg_coarse_size == 24 &&
            amg_config.amg_smoothing_steps == 3 &&
            amg_config.amg_refresh_interval == 5,
        "BiCGSTAB/AMG configuration was not read");

    bool rejected_gmres = false;
    ConfigLine retired_line{2, {"scalarSolver", "gmres", "ilut", "1e-14", "1e-9", "400"}};
    try {
        readLinearSolverLine("tests/data/solution.bs", retired_line, amg_config);
    } catch (const std::exception&) {
        rejected_gmres = true;
    }
    require(rejected_gmres, "retired GMRES configuration was accepted");

    ConfigLine no_preconditioner_line{
        3, {"scalarSolver", "cg", "none", "1e-14", "1e-9", "400"}};
    LinearSolverConfig no_preconditioner_config;
    readLinearSolverLine("tests/data/solution.bs", no_preconditioner_line,
                         no_preconditioner_config);
    no_preconditioner_config.validate();
    require(no_preconditioner_config.preconditioner == PreconditionerType::None,
            "unpreconditioned solver configuration was not read");

    bool rejected_standalone_amg = false;
    retired_line.tokens = {"scalarSolver", "amg", "none", "1e-14", "1e-9", "400"};
    try {
        readLinearSolverLine("tests/data/solution.bs", retired_line, amg_config);
    } catch (const std::exception&) {
        rejected_standalone_amg = true;
    }
    require(rejected_standalone_amg, "standalone AMG configuration was accepted");
    std::cout << "case_io_test: SIMPLE, heat-compatible and transport dictionaries passed\n";
}
