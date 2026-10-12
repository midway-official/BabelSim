#include "babelsim/builtin_descriptions.h"
#include "babelsim/case.h"
#include "babelsim/model_description.h"
#include "internal/field_access.h"
#include "test_util.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace babelsim;
namespace {
struct TemporaryDirectory {
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("babelsim-module-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    TemporaryDirectory() { std::filesystem::create_directory(path); }
    ~TemporaryDirectory() { std::filesystem::remove_all(path); }
};
void write(const std::filesystem::path& path, const std::string& content) {
    std::ofstream output(path);
    output << content;
}
std::string read(const std::filesystem::path& path) {
    std::ifstream input(path);
    std::ostringstream content; content << input.rdbuf(); return content.str();
}
template <typename Function>
void rejected(Function&& function, const std::string& text) {
    try { function(); }
    catch (const std::exception& error) {
        require(std::string(error.what()).find(text) != std::string::npos,
            "wrong module diagnostic: " + std::string(error.what()) + "; expected " + text);
        return;
    }
    throw std::runtime_error("module accepted invalid input; expected " + text);
}
void declareHeat(Case& problem) {
    auto& temperature = problem.load.scalarField("T");
    (void)readEquationControl(problem, "temperature", temperature, {"diffusion"});
}
} // namespace

int main() {
    TemporaryDirectory temporary;
    for (const auto* description : std::vector<const ModuleDescription*>{&builtin::heat(), &builtin::transport(), &builtin::simple(),
        &builtin::transientSimple(), &builtin::piso(), &builtin::spalartAllmaras(),
        &builtin::kOmega(), &builtin::kEpsilon()}) {
        validateModuleDescription(*description);
        require(!description->boundary_types.empty(), "static registration lost boundary metadata");
        const auto json = moduleJson(*description);
        const auto markdown = moduleMarkdown(*description);
        for (const auto& parameter : description->parameters) {
            require(json.find("\"" + parameter.name + "\"") != std::string::npos,
                "JSON omitted a declared parameter");
            require(markdown.find(parameter.name) != std::string::npos,
                "documentation omitted a declared parameter");
        }
    }
    const auto heat = temporary.path / "heat";
    writeCaseTemplate(builtin::heat(), heat);
    const auto original = read(heat / "physics/model.bs");
    require(std::filesystem::is_regular_file(heat / "module.schema.json") &&
        std::filesystem::is_regular_file(heat / "README.md"), "template omitted generated description artifacts");
    rejected([&] { writeCaseTemplate(builtin::heat(), heat); }, "new directory");
    require(read(heat / "physics/model.bs") == original, "existing template was overwritten");
    {
        Case problem(heat, "module-check");
        validateModuleConfiguration(problem, builtin::heat());
        const auto parameters = problem.read.parameters(builtin::heat());
        require(near(parameters.physics.number("density"), 1.0) &&
            near(parameters.physics.number("conductivity"), 0.1), "typed heat values do not match template");
        rejected([&] { (void)parameters.physics.integer("density"); }, "not an integer");
        rejected([&] { (void)problem.fields.scalarField("T"); }, "not been declared");
        require(!problem.read.methods().named.at("equation.temperature").used,
            "description check consumed a numerical profile");
        require(!problem.read.solution().inspect("equation.temperature.kspType").consumed,
            "description check consumed a linear control");
        declareHeat(problem);
        problem.finalizeInitialization();
    }
    {
        Case problem(heat, "typed-config");
        const auto config = problem.read.config(builtin::heat());
        require(near(config.density, 1.0) && near(config.heatCapacity, 1.0) &&
            near(config.conductivity, 0.1) && near(config.source, 1.0), "configuration member bindings differ from template");
        rejected([&] { problem.finalizeInitialization(); }, "required field is not declared: T");
        auto& T = problem.load.scalarField("T");
        rejected([&] { problem.finalizeInitialization(); }, "equation is not bound: temperature");
        auto& wrong = problem.create.scalarField("wrong");
        rejected([&] { (void)equ::createEquation(problem, builtin::heat().equation("temperature"), wrong); }, "requires Case cell unknown T");
        const ScalarField copy = T;
        (void)readEquationControl(problem, "temperature", copy, {"diffusion"});
        rejected([&] { problem.finalizeInitialization(); }, "equation unknown differs: temperature");
    }
    {
        Case problem(heat, "finalized-config");
        (void)problem.read.config(builtin::heat());
        auto& T = problem.load.scalarField("T");
        auto clock = time::makeStepper(problem);
        rejected([&] { clock.advance(); }, "call finalizeInitialization()");
        rejected([&] { problem.write(); }, "call finalizeInitialization()");
        auto equation = equ::createEquation(problem, builtin::heat().equation("temperature"), T);
        const double initialTime = problem.time();
        problem.finalizeInitialization();
        require(problem.time() == initialTime && problem.step() == 0 && clock.step() == 0,
            "finalization advanced physical time");
        require(&problem.fields.scalarField("T") == &T, "finalization invalidated fields");
        rejected([&] { problem.create.scalarField("tooLate"); }, "declarations are closed");
        rejected([&] { (void)problem.read.config(builtin::heat()); }, "declarations are closed");
        rejected([&] { (void)equ::createEquation(problem, "lateEquation", T); }, "declarations are closed");
        clock.advance();
        require(problem.step() == 1, "finalized clock cannot advance");
    }
    {
        struct Config { struct Material { double diffusivity; } material; bool active; int count; std::string method; std::optional<int> optionalCount; };
        ConfigurationDescription<Config> description;
        description.name = "typedScalar";
        ParameterDescription diffusion;
        diffusion.name = "diffusivity"; diffusion.required = false; diffusion.default_value = 0.25;
        description.bind(diffusion, &Config::material, &Config::Material::diffusivity);
        ParameterDescription enabled; enabled.name = "active"; enabled.required = false; enabled.default_value = true;
        description.bind(enabled, &Config::active);
        auto count = enabled; count.name = "count"; count.default_value = 3;
        description.bind(count, &Config::count);
        auto method = enabled; method.name = "method"; method.default_value = std::string("implicit");
        description.bind(method, &Config::method);
        auto optional = count; optional.name = "optionalCount"; optional.default_value.reset();
        description.bind(optional, &Config::optionalCount);
        Case problem(heat, "member-types");
        description.parameters[0].name = "diffusionCoefficient"; // copied metadata remains the single key source
        const auto config = problem.read.config(description);
        require(near(config.material.diffusivity, 0.25) && config.active && config.count == 3 &&
            config.method == "implicit" && !config.optionalCount, "typed optional/default binding failed");
        require(description.parameters[1].type == ParameterType::Boolean &&
            description.parameters[2].type == ParameterType::Integer && description.parameters[3].type == ParameterType::Word,
            "member types did not supply schema types");
        rejected([&] { description.bind(optional, &Config::count); }, "std::optional member");
        description.parameters.push_back(diffusion);
        rejected([&] { (void)description.configuration({}); }, "declared with bind()");
    }
    write(heat / "physics/model.bs", "density 0\nheatCapacity 1\nconductivity 0.1\nsource 1\n");
    {
        Case problem(heat, "range-check");
        rejected([&] { (void)problem.read.parameters(builtin::heat()); }, "model.bs:1: density");
    }
    write(heat / "physics/model.bs", original + "unrecognized 123\n");
    {
        Case problem(heat, "unused-check");
        validateModuleConfiguration(problem, builtin::heat());
        declareHeat(problem);
        rejected([&] { problem.finalizeInitialization(); }, "unrecognized");
    }
    write(heat / "physics/model.bs", original);
    const auto field = read(heat / "fields/initial/T.field");
    write(heat / "fields/initial/T.field", "field T\n{ type vector location cell internal uniform (0 0 0) }\n");
    {
        Case problem(heat, "field-check");
        rejected([&] { validateModuleConfiguration(problem, builtin::heat()); }, "input header differs");
    }
    write(heat / "fields/initial/T.field", field);
    {
        auto obsolete = field;
        const auto position = obsolete.find("zeroGradient");
        require(position != std::string::npos, "template omitted canonical boundary name");
        obsolete.replace(position, std::string("zeroGradient").size(), "mirror");
        write(heat / "fields/initial/T.field", obsolete);
        Case problem(heat, "obsolete-boundary-check");
        rejected([&] { validateModuleConfiguration(problem, builtin::heat()); }, "does not support boundary type mirror");
        rejected([&] { (void)problem.load.scalarField("T"); }, "unknown boundary condition mirror");
        write(heat / "fields/initial/T.field", field);
    }
    for (const std::string filename : {"type", "location"}) {
        auto nonuniform = field;
        const auto position = nonuniform.find("internal uniform (0)");
        require(position != std::string::npos, "template initial value differs from test fixture");
        nonuniform.replace(position, std::string("internal uniform (0)").size(), "internal file " + filename);
        write(heat / "fields/initial/T.field", nonuniform);
        write(heat / "fields/initial" / filename, "0 2.5\n");
        Case problem(heat, "keyword-filename-check");
        validateModuleConfiguration(problem, builtin::heat());
        declareHeat(problem);
        require(near(detail::fieldData(problem.fields.scalarField("T"))[0], 2.5),
            "header precheck changed nonuniform field reading");
        problem.finalizeInitialization();
    }
    write(heat / "fields/initial/T.field", field);
    {
        Case problem(heat, "boundary-check");
        auto limited = builtin::heat(); limited.boundary_types = {"fixedValue"}; limited.supports_case_template = false;
        rejected([&] { validateModuleConfiguration(problem, limited); }, "does not support boundary type zeroGradient");
    }
    auto methods = read(heat / "numerics/methods.bs");
    write(heat / "numerics/methods.bs", "time steady\n" + methods.substr(methods.find('\n') + 1));
    {
        Case problem(heat, "time-check");
        rejected([&] { validateModuleConfiguration(problem, builtin::heat()); }, "does not support time steady");
    }
    write(heat / "numerics/methods.bs", methods);
    {
        Case problem(heat, "default-check");
        const auto sa = readModuleParameters(problem, builtin::spalartAllmaras());
        const double expected = 0.1355 / (0.41 * 0.41) + (1.0 + 0.622) / (2.0 / 3.0);
        require(near(sa.physics.number("saCw1"), expected), "SA dependent default changed");
        require(near(sa.solution.number("turbulenceRelaxation"), 0.7), "model fraction default changed");
        rejected([&] { validateModuleConfiguration(problem, builtin::spalartAllmaras()); }, "incompatible with solver");
    }
    write(heat / "physics/model.bs", original + "saCb1 0.2\nsaKappa 0.5\n");
    {
        Case problem(heat, "dependent-default-check");
        const auto sa = readModuleParameters(problem, builtin::spalartAllmaras());
        require(near(sa.physics.number("saCw1"), 0.2/(0.5*0.5) + (1.0+0.622)/(2.0/3.0)),
            "dependent default ignored configured model constants");
    }
    write(heat / "physics/model.bs", original);
    {
        write(heat / "physics/model.bs", "density 1\ndynamicViscosity 0.01\nturbulenceModel laminar\n");
        Case problem(heat, "obsolete-model-check");
        rejected([&] { (void)problem.read.parameters(builtin::piso()); }, "outside declared choices");
        write(heat / "physics/model.bs", original);
    }
    {
        auto custom = builtin::heat(); custom.name = "customHeat";
        auto count = ParameterDescription{};
        count.name = "minimumCorrections"; count.dictionary = ParameterDictionary::Solution;
        count.type = ParameterType::Integer; count.required = false; count.default_value = 3;
        auto limit = count; limit.name = "maximumCorrections"; limit.default_value = 2;
        custom.parameters.push_back(count); custom.parameters.push_back(limit);
        custom.parameter_constraints.push_back({ParameterDictionary::Solution, count.name,
            ParameterRelation::LessEqual, limit.name, "The limit must cover requested corrections."});
        Case problem(heat, "relation-check");
        rejected([&] { (void)problem.read.parameters(custom); }, "minimumCorrections <= maximumCorrections required");
    }
    {
        auto malformed = builtin::heat(); malformed.parameters.push_back(malformed.parameters.front());
        rejected([&] { validateModuleDescription(malformed); }, "duplicate module parameter");
        malformed = builtin::heat(); malformed.equations[0].unknown = "undeclared";
        rejected([&] { validateModuleDescription(malformed); }, "undeclared unknown");
        malformed = builtin::heat(); malformed.boundary_types = {"fixedValue"};
        rejected([&] { writeCaseTemplate(malformed, temporary.path / "unsupported"); }, "requires zeroGradient");
        require(!std::filesystem::exists(temporary.path / "unsupported"), "invalid tutorial created partial output");
        malformed = builtin::heat(); malformed.fields[0].template_value = "not-a-number";
        rejected([&] { writeCaseTemplate(malformed, temporary.path / "invalid-field"); }, "finite components");
        require(!std::filesystem::exists(temporary.path / "invalid-field"), "invalid field tutorial created partial output");
    }
    ModuleDescription custom = builtin::heat(); custom.name = "describedExternal";
    {
        ModuleRegistration registration(custom);
        const auto registered = registeredModuleDescriptions();
        require(std::find(registered.begin(), registered.end(), &custom) != registered.end(), "external model was not discoverable");
        {
            ModuleRegistration duplicate(custom);
            rejected([&] { (void)registeredModuleDescriptions(); }, "duplicate module registration");
        }
    }
    const auto registered = registeredModuleDescriptions();
    require(std::find(registered.begin(), registered.end(), &custom) == registered.end(),
        "destroyed registration remained discoverable");
    {
        auto multiline = builtin::heat();
        multiline.parameters[0].explanation = "Density | material\nUse SI units.";
        multiline.parameters[0].unit = "kg/m^3\r\nSI";
        require(moduleMarkdown(multiline).find("Density \\| material<br>Use SI units.") != std::string::npos,
            "Markdown table did not escape multiline descriptions");
        const auto directory = temporary.path / "multiline";
        writeCaseTemplate(multiline, directory);
        Case problem(directory, "multiline-check");
        validateModuleConfiguration(problem, multiline);
        declareHeat(problem);
        problem.finalizeInitialization();
    }
    std::cout << "module_description_test: generated templates, typed defaults, frontend checks, and registry passed\n";
}
