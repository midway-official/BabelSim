#include "babelsim/model_description.h"
#include "module_description_support.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace babelsim {
using namespace detail::module_description;

namespace {

std::string singleLine(std::string text) {
    std::replace(text.begin(), text.end(), '\n', ' ');
    std::replace(text.begin(), text.end(), '\r', ' ');
    return text;
}
void writeFile(const std::filesystem::path& path, const std::string& content) {
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("cannot write " + path.string());
    output << content;
    output.close();
    if (!output)
        throw std::runtime_error("cannot finish writing " + path.string());
}

} // namespace

void writeCaseTemplate(const ModuleDescription& description,
                       const std::filesystem::path& directory) {
    validateModuleDescription(description);
    if (!description.supports_case_template)
        throw std::invalid_argument("module " + description.name + " has no generic case template");
    if (directory.empty() || std::filesystem::exists(directory) ||
        std::filesystem::is_symlink(directory))
        throw std::invalid_argument("case template output must be a new directory: " +
                                    directory.string());
    if (!directory.parent_path().empty())
        std::filesystem::create_directories(directory.parent_path());
    if (!std::filesystem::create_directory(directory))
        throw std::invalid_argument("case template output already exists: " + directory.string());
    try {
        for (const auto* child : {"mesh", "fields/initial", "physics", "numerics"})
            std::filesystem::create_directories(directory / child);
        writeFile(directory / "case.bs",
                  "solver " + description.name +
                      "\nmesh mesh/domain.mesh\nfields fields/initial\nphysics physics/model.bs\n"
                      "methods numerics/methods.bs\nsolution numerics/solution.bs\ncontrol "
                      "control.bs\noutput output.bs\n");
        writeFile(
            directory / "mesh/domain.mesh",
            "BABELSIM_MESH 3\nvertices 8\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n0 0 1\n1 0 1\n1 1 1\n0 1 1\n"
            "faces 6\nface 4 0 -1 0 0 4 7 3\nface 4 0 -1 1 1 2 6 5\nface 4 0 -1 2 0 1 5 4\n"
            "face 4 0 -1 3 3 7 6 2\nface 4 0 -1 4 0 3 2 1\nface 4 0 -1 5 4 5 6 7\n"
            "patches 6\npatch xmin generic\npatch xmax generic\npatch ymin generic\npatch ymax "
            "generic\n"
            "patch zmin generic\npatch zmax generic\nend\n");
        std::ostringstream physics, solution, methods;
        ModuleParameters tutorial;
        for (const auto& parameter : description.parameters) {
            auto& values = valuesFor(tutorial, parameter.dictionary);
            std::optional<ParameterValue> value =
                parameter.template_value ? parameter.template_value : parameter.default_value;
            if (!value && parameter.default_factory)
                value = parameter.default_factory(values);
            if (!value)
                continue;
            checkValue(parameter, *value, "tutorial " + parameter.name);
            values.insert(parameter.name, *value);
            auto& output =
                parameter.dictionary == ParameterDictionary::Physics ? physics : solution;
            output << "# " << singleLine(parameter.explanation + " [" + parameter.unit + "]")
                   << '\n'
                   << parameter.name << ' ' << valueText(*value) << '\n';
        }
        checkRelations(description, tutorial, "tutorial " + description.name);
        methods << "time " << timeName(description.time_methods.front()) << '\n';
        for (const auto& equation : description.equations) {
            methods << "equation." << equation.name << ".interpolation linear\n"
                    << "equation." << equation.name << ".gradient greenGauss\n"
                    << "equation." << equation.name << ".convection upwind\n"
                    << "equation." << equation.name << ".diffusion orthogonal\n";
            solution << "equation." << equation.name << ".kspType bcgs\n"
                     << "equation." << equation.name << ".pcType jacobi\n"
                     << "equation." << equation.name << ".absoluteTolerance 1e-14\n"
                     << "equation." << equation.name << ".relativeTolerance 1e-10\n"
                     << "equation." << equation.name << ".maxIterations 1000\n";
        }
        writeFile(directory / "physics/model.bs", physics.str());
        writeFile(directory / "numerics/solution.bs", solution.str());
        writeFile(directory / "numerics/methods.bs", methods.str());
        writeFile(directory / "control.bs", "startTime 0\nendTime 0.05\ndeltaT 0.01\n");
        writeFile(directory / "output.bs", "directory results\ntimeName final\nwriteInterval 1\n");
        for (const auto& field : description.fields) {
            if (field.initialization != FieldInitialization::File)
                continue;
            std::ostringstream output;
            output << "field " << field.name << "\n{\n    type " << fieldTypeName(field.type)
                   << "\n    location cell\n    internal uniform (" << field.template_value
                   << ")\n    boundary\n    {\n";
            for (const auto* patch : {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"})
                output << "        " << patch << " { type zeroGradient }\n";
            output << "    }\n}\n";
            writeFile(directory / "fields/initial" / (field.name + ".field"), output.str());
        }
        writeFile(directory / "README.md",
                  moduleMarkdown(description) +
                      "\nThis unit-cube tutorial tests configuration and development workflow. "
                      "It is not a physical validation case. Run with `babelsim-solve -case "
                      "<this-directory>`.\n");
        writeFile(directory / "module.schema.json", moduleJson(description));
    } catch (...) {
        std::filesystem::remove_all(directory);
        throw;
    }
}
} // namespace babelsim
