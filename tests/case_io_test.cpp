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
    const auto transport_options = methods.equationOptions("transport");
    require(transport_options.convection == ConvectionMethod::Upwind,
            "transport equation convection method was not read");
    const Methods overrides = readMethodsFile("tests/data/methods.bs");
    const auto temperature_options = overrides.equationOptions("temperature");
    const auto concentration_options = overrides.equationOptions("concentration");
    require(temperature_options.interpolation == InterpolationMethod::Corrected &&
            temperature_options.gradient == GradientMethod::LeastSquares &&
            temperature_options.convection == ConvectionMethod::LinearUpwind &&
            temperature_options.diffusion == DiffusionMethod::LimitedCorrected,
            "equation numerical methods were mixed between operator types");
    require(concentration_options.interpolation == InterpolationMethod::Linear &&
            concentration_options.gradient == GradientMethod::GreenGauss &&
            concentration_options.convection == ConvectionMethod::Upwind &&
            concentration_options.diffusion == DiffusionMethod::Orthogonal &&
            overrides.time == TimeMethod::Euler,
        "equation methods changed the configured or time method");
    overrides.requireAllUsed();

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

    // A named equation owns its complete linear configuration.
    const auto override_case = std::filesystem::temp_directory_path() /
        "babelsim_linear_override_case";
    std::filesystem::remove_all(override_case);
    std::filesystem::copy("cases/heat", override_case,
                          std::filesystem::copy_options::recursive);
    {
        std::ofstream file(override_case / "numerics/solution.bs");
        file << "equation.temperature.solver cg\n"
             << "equation.temperature.preconditioner incompleteCholesky\n"
             << "equation.temperature.absoluteTolerance 1e-13\n"
             << "equation.temperature.relativeTolerance 2e-9\n"
             << "equation.temperature.maxIterations 321\n";
    }
    {
        Case custom(override_case);
        const auto& temperature = custom.scalarField("T");
        const auto selected = readEquationControl(custom, "temperature", temperature);
        require(selected.linear.solver == LinearSolverType::ConjugateGradient &&
                selected.linear.preconditioner == PreconditionerType::IncompleteCholesky &&
                near(selected.linear.absolute_tolerance, 1e-13) &&
                near(selected.linear.relative_tolerance, 2e-9) && selected.linear.max_iterations == 321,
                "named equation linear control was not applied");
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
        readEquationControl(custom, "temperature", custom.scalarField("T"));
        custom.createScalarField("derived", 2.0);
        custom.write();
        const auto result_directory = output_case / "results/0/rank-0000";
        require(std::filesystem::exists(result_directory / "T.csv") &&
                std::filesystem::exists(result_directory / "derived.csv") &&
                !std::filesystem::exists(result_directory / "U.csv"),
                "output.bs field selection did not control the written cell fields");
    }
    std::filesystem::remove_all(output_case);

    const auto old_methods = std::filesystem::temp_directory_path() /
        "babelsim_old_methods_syntax.bs";
    {
        std::ofstream file(old_methods);
        file << "interpolation linear\ngradient leastSquares\ntime euler\n";
    }
    bool rejected_old_methods = false;
    try { (void)readMethodsFile(old_methods); }
    catch (const std::exception&) { rejected_old_methods = true; }
    require(rejected_old_methods, "legacy global methods syntax was accepted");
    std::filesystem::remove(old_methods);
    std::cout << "case_io_test: SIMPLE, heat-compatible and transport dictionaries passed\n";
}
