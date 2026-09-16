#include "babelsim/case.h"
#include "babelsim/numerics_io.h"
#include "babelsim/field_io.h"
#include "babelsim/mesh_io.h"

#include "test_util.h"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <limits>

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

    // Typed dictionary access keeps Physics code away from raw tokens and
    // exposes whether a value came from the file without consuming it.
    const auto parameter_path = std::filesystem::temp_directory_path() /
        "babelsim_parameter_api_test.bs";
    {
        std::ofstream file(parameter_path);
        file << "name kOmega\nactive yes\niterations 12\nstrength 2.5\n";
    }
    const Parameters typed(parameter_path);
    const auto initial_state = typed.inspect("name");
    require(initial_state.configured && !initial_state.consumed && initial_state.line == 1,
            "parameter inspection should be read-only");
    require(typed.word("name") == "kOmega" && typed.boolean("active") &&
            typed.integer("iterations") == 12 && near(typed.number("strength"), 2.5),
            "typed parameter access returned an unexpected value");
    require(typed.inspect("name").consumed, "parameter consumption was not tracked");
    require(typed.word("missing", "fallback") == "fallback" &&
            typed.boolean("missingBool", false) == false &&
            typed.integer("missingInt", 3) == 3 &&
            near(typed.nonnegative("missingNonnegative", 0.25), 0.25),
            "typed fallback access failed");
    bool invalid_default = false;
    try { (void)typed.number("missingNumber", std::numeric_limits<double>::quiet_NaN()); }
    catch (const std::invalid_argument&) { invalid_default = true; }
    require(invalid_default, "non-finite numeric defaults were accepted");
    typed.requireAllUsed();
    std::filesystem::remove(parameter_path);

    const auto output_path = std::filesystem::temp_directory_path() /
        "babelsim_output_selection_test.bs";
    {
        std::ofstream file(output_path);
        file << "directory results\ntimeName final\nwriteInterval 4\n"
             << "writeFields T p\nexcludeFields p\n";
    }
    CaseDefinition output_definition = cavity;
    output_definition.output_file = output_path;
    const OutputControl output = readOutputControl(output_definition);
    require(output.write_interval == 4 && output.write_fields == std::vector<std::string>{"T", "p"} &&
            output.exclude_fields == std::vector<std::string>{"p"},
            "output field selection was not parsed");
    std::filesystem::remove(output_path);

    // A field-specific linear entry is selected by the field overload, while
    // the base scalarSolver entry remains the runtime default.
    const auto override_case = std::filesystem::temp_directory_path() /
        "babelsim_linear_override_case";
    std::filesystem::remove_all(override_case);
    std::filesystem::copy("cases/heat", override_case,
                          std::filesystem::copy_options::recursive);
    {
        std::ofstream file(override_case / "numerics/solution.bs", std::ios::app);
        file << "scalarSolver.T cg incompleteCholesky 1e-13 2e-9 321\n";
    }
    {
        Case custom(override_case);
        const auto& temperature = custom.scalarField("T");
        const auto selected = readLinearControl(custom, temperature);
        require(selected.solver == LinearSolverType::ConjugateGradient &&
                selected.preconditioner == PreconditionerType::IncompleteCholesky &&
                near(selected.absolute_tolerance, 1e-13) &&
                near(selected.relative_tolerance, 2e-9) && selected.max_iterations == 321,
                "field-specific linear control was not applied");
    }
    std::filesystem::remove_all(override_case);

    const auto output_case = std::filesystem::temp_directory_path() /
        "babelsim_output_selection_case";
    std::filesystem::remove_all(output_case);
    std::filesystem::copy("cases/heat", output_case,
                          std::filesystem::copy_options::recursive);
    {
        std::ofstream file(output_case / "output.bs");
        file << "directory results\ntimeName final\nwriteInterval 1\n"
             << "writeFields T derived\n";
    }
    {
        Case custom(output_case);
        custom.physics().positive("density");
        custom.physics().positive("heatCapacity");
        custom.physics().nonnegative("conductivity");
        custom.physics().number("source");
        custom.scalarField("T");
        custom.createScalarField("derived", 2.0);
        custom.write();
        const auto result_directory = output_case / "results/0/rank-0000";
        require(std::filesystem::exists(result_directory / "T.csv") &&
                std::filesystem::exists(result_directory / "derived.csv") &&
                !std::filesystem::exists(result_directory / "U.csv"),
                "output.bs field selection did not control the written cell fields");
    }
    std::filesystem::remove_all(output_case);

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
