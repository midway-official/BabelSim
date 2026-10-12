#pragma once

#include "babelsim/config.h"
#include "babelsim/methods.h"

#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace babelsim {
class Case;

enum class ParameterDictionary { Physics, Solution };
enum class ParameterType { Number, Integer, Boolean, Word };
using ParameterValue = std::variant<double, int, bool, std::string>;

// Parsed values have a stable value type. Reading them does not parse files,
// consume additional configuration, or change fields and time state.
class ParameterValues {
public:
    bool contains(const std::string& name) const;
    double number(const std::string& name) const;
    int integer(const std::string& name) const;
    bool boolean(const std::string& name) const;
    const std::string& word(const std::string& name) const;
    const ParameterValue& value(const std::string& name) const;
    void insert(const std::string& name, ParameterValue value);
private:
    std::map<std::string, ParameterValue> values_;
};

struct ParameterDescription {
    std::string name;
    ParameterDictionary dictionary = ParameterDictionary::Physics;
    ParameterType type = ParameterType::Number;
    std::string unit = "1";
    std::string explanation;
    bool required = true;
    std::optional<ParameterValue> default_value;
    // Required parameters may have a tutorial value without making it a default.
    std::optional<ParameterValue> template_value;
    std::optional<double> minimum;
    std::optional<double> maximum;
    bool exclusive_minimum = false;
    bool exclusive_maximum = false;
    std::vector<std::string> choices;
    // Evaluate an optional default from previously declared values in the SAME
    // dictionary. This preserves model defaults that depend on other constants.
    ParameterValue (*default_factory)(const ParameterValues&) = nullptr;
    std::string default_explanation;
};

enum class DescribedFieldType { Scalar, Vector, Tensor };
enum class DescribedFieldLocation { Cell, Face };
enum class FieldRole { Input, Derived };
struct FieldDescription {
    std::string name;
    DescribedFieldType type = DescribedFieldType::Scalar;
    DescribedFieldLocation location = DescribedFieldLocation::Cell;
    FieldRole role = FieldRole::Input;
    std::string unit = "1";
    std::string explanation;
    // Values use the existing .field literal syntax without parentheses.
    std::string template_value = "0";
};

struct EquationDescription {
    std::string name;
    std::string unknown;
    std::vector<std::string> terms;
    std::string explanation;
    std::vector<std::string> numerical_options{
        "interpolation", "gradient", "convection", "diffusion"};
};

enum class ParameterRelation { LessEqual, GreaterEqual, Equal };
struct ParameterConstraint {
    ParameterDictionary dictionary = ParameterDictionary::Solution;
    std::string left;
    ParameterRelation relation = ParameterRelation::GreaterEqual;
    std::string right;
    std::string explanation;
};

struct ModuleDescription {
    std::string name;
    std::string explanation;
    std::vector<FieldDescription> fields;
    std::vector<ParameterDescription> parameters;
    std::vector<EquationDescription> equations;
    std::vector<TimeMethod> time_methods;
    std::vector<std::string> boundary_types;
    std::vector<std::string> compatible_solvers;
    std::vector<ParameterConstraint> parameter_constraints;
    std::vector<std::string> limitations;
    // The generic tutorial writer is opt-in: it creates a unit cube with the
    // input fields, complete numerical profiles, and the declared defaults.
    bool supports_case_template = false;
};

struct ModuleParameters {
    ParameterValues physics;
    ParameterValues solution;
};

// Optional metadata-only registration for composable models without a top-level
// solver function. The description must outlive this registration.
class ModuleRegistration {
public:
    explicit ModuleRegistration(const ModuleDescription&) noexcept;
    ~ModuleRegistration() noexcept;
    ModuleRegistration(const ModuleRegistration&) = delete;
    ModuleRegistration& operator=(const ModuleRegistration&) = delete;
private:
    const ModuleDescription* description_;
    mutable const ModuleRegistration* next_;
    static const ModuleRegistration*& first() noexcept;
    friend std::vector<const ModuleDescription*> registeredModuleDescriptions();
};
std::vector<const ModuleDescription*> registeredModuleDescriptions();

void validateModuleDescription(const ModuleDescription&);
ModuleParameters readModuleParameters(const Case&, const ModuleDescription&);
// Frontend checks only: parameter values, time/solver compatibility, input
// headers, and declared profiles. It never loads/creates fields or marks
// numerical profiles used. Case::validate() still enforces actual consumption.
void validateModuleConfiguration(const Case&, const ModuleDescription&);
std::string moduleMarkdown(const ModuleDescription&);
std::string moduleJson(const ModuleDescription&);
enum class ModuleOutputFormat { Markdown, Json };
struct ModuleListing { std::string name; bool solver; bool described; };
// Explicit CLI output remains owned by IO. Pure string generators above are
// suitable for embedded applications, files, and GUI consumers.
void writeModuleDescription(const ModuleDescription&, ModuleOutputFormat);
void writeModuleList(const std::vector<ModuleListing>&);
void writeCaseTemplate(const ModuleDescription&, const std::filesystem::path& directory);

} // namespace babelsim
