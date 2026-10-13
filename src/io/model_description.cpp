#include "babelsim/model_description.h"
#include "babelsim/case.h"
#include "module_description_support.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace babelsim {
namespace detail::module_description {

const char* dictionaryName(ParameterDictionary dictionary) {
    return dictionary == ParameterDictionary::Physics ? "physics" : "solution";
}
const char* typeName(ParameterType type) {
    switch (type) {
    case ParameterType::Number:
        return "number";
    case ParameterType::Integer:
        return "integer";
    case ParameterType::Boolean:
        return "boolean";
    case ParameterType::Word:
        return "string";
    }
    throw std::invalid_argument("invalid parameter type");
}
const char* fieldRoleName(FieldRole role) {
    switch (role) {
    case FieldRole::Input:
        return "input";
    case FieldRole::Output:
        return "output";
    case FieldRole::Work:
        return "work";
    }
    throw std::invalid_argument("invalid field role");
}
const char* initializationName(FieldInitialization source) {
    switch (source) {
    case FieldInitialization::File:
        return "file";
    case FieldInitialization::Provided:
        return "provided";
    case FieldInitialization::Constant:
        return "constant";
    case FieldInitialization::Computed:
        return "computed";
    case FieldInitialization::HomogeneousCorrection:
        return "homogeneousCorrection";
    }
    throw std::invalid_argument("invalid field initialization source");
}
const char* fieldTypeName(DescribedFieldType type) {
    switch (type) {
    case DescribedFieldType::Scalar:
        return "scalar";
    case DescribedFieldType::Vector:
        return "vector";
    case DescribedFieldType::Tensor:
        return "tensor";
    }
    throw std::invalid_argument("invalid field type");
}
const char* timeName(TimeMethod method) {
    switch (method) {
    case TimeMethod::Steady:
        return "steady";
    case TimeMethod::Euler:
        return "euler";
    case TimeMethod::BDF2:
        return "bdf2";
    }
    throw std::invalid_argument("invalid time method");
}
const char* relationName(ParameterRelation relation) {
    switch (relation) {
    case ParameterRelation::LessEqual:
        return "<=";
    case ParameterRelation::GreaterEqual:
        return ">=";
    case ParameterRelation::Equal:
        return "==";
    }
    throw std::invalid_argument("invalid parameter relation");
}
bool safeName(const std::string& name) {
    if (name.empty())
        return false;
    return std::all_of(name.begin(), name.end(),
                       [](unsigned char character) {
                           return std::isalnum(character) || character == '_' || character == '.' ||
                                  character == '-';
                       }) &&
           name != "." && name != "..";
}
void requireName(const std::string& name, const std::string& context) {
    if (!safeName(name))
        throw std::invalid_argument("invalid " + context + " name: " + name);
}

bool hasType(const ParameterValue& value, ParameterType type) {
    return (type == ParameterType::Number && std::holds_alternative<double>(value)) ||
           (type == ParameterType::Integer && std::holds_alternative<int>(value)) ||
           (type == ParameterType::Boolean && std::holds_alternative<bool>(value)) ||
           (type == ParameterType::Word && std::holds_alternative<std::string>(value));
}
double numericValue(const ParameterValue& value) {
    if (const auto* number = std::get_if<double>(&value))
        return *number;
    if (const auto* integer = std::get_if<int>(&value))
        return *integer;
    throw std::invalid_argument("parameter relation needs number/integer values");
}
void checkValue(const ParameterDescription& parameter, const ParameterValue& value,
                const std::string& context) {
    if (!hasType(value, parameter.type))
        throw std::invalid_argument(context + ": expected " + typeName(parameter.type));
    if (parameter.type == ParameterType::Number || parameter.type == ParameterType::Integer) {
        const double number = numericValue(value);
        if (!std::isfinite(number))
            throw std::invalid_argument(context + ": expected a finite number");
        if (parameter.minimum && (parameter.exclusive_minimum ? number <= *parameter.minimum
                                                              : number < *parameter.minimum))
            throw std::invalid_argument(context + ": value below declared minimum " +
                                        valueText(*parameter.minimum) +
                                        (parameter.exclusive_minimum ? " (exclusive)" : ""));
        if (parameter.maximum && (parameter.exclusive_maximum ? number >= *parameter.maximum
                                                              : number > *parameter.maximum))
            throw std::invalid_argument(context + ": value above declared maximum " +
                                        valueText(*parameter.maximum) +
                                        (parameter.exclusive_maximum ? " (exclusive)" : ""));
    } else if (parameter.type == ParameterType::Word) {
        const auto& word = std::get<std::string>(value);
        if (word.empty() || std::any_of(word.begin(), word.end(), [](unsigned char character) {
                return std::isspace(character) || character == '#';
            }))
            throw std::invalid_argument(context + ": expected one nonempty word");
        if (!parameter.choices.empty() &&
            std::find(parameter.choices.begin(), parameter.choices.end(), word) ==
                parameter.choices.end())
            throw std::invalid_argument(context + ": value outside declared choices");
    }
}
ParameterValue readValue(const Parameters& input, const ParameterDescription& parameter) {
    switch (parameter.type) {
    case ParameterType::Number:
        return input.number(parameter.name);
    case ParameterType::Integer:
        return input.integer(parameter.name);
    case ParameterType::Boolean:
        return input.boolean(parameter.name);
    case ParameterType::Word:
        return input.word(parameter.name);
    }
    throw std::invalid_argument("invalid parameter type");
}
ParameterValues& valuesFor(ModuleParameters& values, ParameterDictionary dictionary) {
    return dictionary == ParameterDictionary::Physics ? values.physics : values.solution;
}
const ParameterValues& valuesFor(const ModuleParameters& values, ParameterDictionary dictionary) {
    return dictionary == ParameterDictionary::Physics ? values.physics : values.solution;
}
void checkRelations(const ModuleDescription& description, const ModuleParameters& values,
                    const std::string& context) {
    for (const auto& constraint : description.parameter_constraints) {
        const auto& dictionary = valuesFor(values, constraint.dictionary);
        // An optional unset setting remains solver-controlled. A relation is
        // checked as soon as both effective values are available.
        if (!dictionary.contains(constraint.left) || !dictionary.contains(constraint.right))
            continue;
        const double left = numericValue(dictionary.value(constraint.left));
        const double right = numericValue(dictionary.value(constraint.right));
        const bool valid =
            constraint.relation == ParameterRelation::LessEqual
                ? left <= right
                : (constraint.relation == ParameterRelation::GreaterEqual ? left >= right
                                                                          : left == right);
        if (!valid)
            throw std::invalid_argument(context + ": " + constraint.left + " " +
                                        relationName(constraint.relation) + " " + constraint.right +
                                        " required; " + constraint.explanation);
    }
}
void checkTutorialLiteral(const FieldDescription& field) {
    const int components = field.type == DescribedFieldType::Scalar
                               ? 1
                               : (field.type == DescribedFieldType::Vector ? 3 : 9);
    std::istringstream input(field.template_value);
    for (int component = 0; component < components; ++component) {
        double value = 0.0;
        if (!(input >> value) || !std::isfinite(value))
            throw std::invalid_argument("field tutorial requires " + std::to_string(components) +
                                        " finite components: " + field.name);
    }
    std::string extra;
    if (input >> extra)
        throw std::invalid_argument("unexpected field tutorial component: " + field.name);
}
void checkFieldHeader(const std::filesystem::path& path, const FieldDescription& field,
                      const std::vector<std::string>& boundaries) {
    const auto tokens = readConfigTokens(path);
    const auto expectedLocation = field.location == DescribedFieldLocation::Cell ? "cell" : "face";
    // Only inspect the top-level header. The existing field reader remains the
    // authority for values, patches, aliases, and the complete field grammar.
    bool type = false, location = false;
    int depth = 0;
    for (std::size_t i = 2; i < tokens.size(); ++i) {
        if (tokens[i].text == "{")
            ++depth;
        else if (tokens[i].text == "}")
            --depth;
        else if (depth == 1 && i + 1 < tokens.size()) {
            if (tokens[i].text == "type")
                type = tokens[++i].text == fieldTypeName(field.type);
            else if (tokens[i].text == "location")
                location = tokens[++i].text == expectedLocation;
            else if (tokens[i].text == "internal" && tokens[i + 1].text == "file" &&
                     i + 2 < tokens.size())
                i += 2; // A data filename may itself be a header keyword.
        } else if (depth == 3 && tokens[i].text == "type" && i + 1 < tokens.size() &&
                   !boundaries.empty()) {
            const auto& token = tokens[++i];
            if (std::find(boundaries.begin(), boundaries.end(), token.text) == boundaries.end())
                throw std::invalid_argument(path.string() + ":" + std::to_string(token.line) +
                                            ": module does not support boundary type " +
                                            token.text);
        }
    }
    if (tokens.size() < 8 || tokens[0].text != "field" || tokens[1].text != field.name ||
        tokens[2].text != "{" || !type || !location)
        throw std::invalid_argument(path.string() + ": input header differs from module field " +
                                    field.name);
}
} // namespace detail::module_description

using namespace detail::module_description;

const ModuleRegistration*& ModuleRegistration::first() noexcept {
    static const ModuleRegistration* head = nullptr;
    return head;
}
ModuleRegistration::ModuleRegistration(const ModuleDescription& description) noexcept
    : description_(&description), next_(first()) {
    first() = this;
}
ModuleRegistration::~ModuleRegistration() noexcept {
    const ModuleRegistration** link = &first();
    while (*link && *link != this)
        link = &(*link)->next_;
    if (*link)
        *link = next_;
}
std::vector<const ModuleDescription*> registeredModuleDescriptions() {
    std::vector<const ModuleDescription*> result;
    std::set<std::string> names;
    for (auto* entry = ModuleRegistration::first(); entry; entry = entry->next_) {
        validateModuleDescription(*entry->description_);
        if (!names.insert(entry->description_->name).second)
            throw std::invalid_argument("duplicate module registration: " +
                                        entry->description_->name);
        result.push_back(entry->description_);
    }
    return result;
}

const EquationDescription& ModuleDescription::equation(const std::string& equationName) const {
    const auto found =
        std::find_if(equations.begin(), equations.end(),
                     [&](const EquationDescription& value) { return value.name == equationName; });
    if (found == equations.end())
        throw std::invalid_argument("undeclared module equation: " + name + "." + equationName);
    return *found;
}

bool ParameterValues::contains(const std::string& name) const {
    return values_.count(name) != 0;
}
const ParameterValue& ParameterValues::value(const std::string& name) const {
    const auto found = values_.find(name);
    if (found == values_.end())
        throw std::invalid_argument("undeclared or unset module parameter: " + name);
    return found->second;
}
double ParameterValues::number(const std::string& name) const {
    const auto* result = std::get_if<double>(&value(name));
    if (!result)
        throw std::invalid_argument("module parameter is not a number: " + name);
    return *result;
}
int ParameterValues::integer(const std::string& name) const {
    const auto* result = std::get_if<int>(&value(name));
    if (!result)
        throw std::invalid_argument("module parameter is not an integer: " + name);
    return *result;
}
bool ParameterValues::boolean(const std::string& name) const {
    const auto* result = std::get_if<bool>(&value(name));
    if (!result)
        throw std::invalid_argument("module parameter is not a boolean: " + name);
    return *result;
}
const std::string& ParameterValues::word(const std::string& name) const {
    const auto* result = std::get_if<std::string>(&value(name));
    if (!result)
        throw std::invalid_argument("module parameter is not a word: " + name);
    return *result;
}
void ParameterValues::insert(const std::string& name, ParameterValue value) {
    if (!values_.emplace(name, std::move(value)).second)
        throw std::invalid_argument("duplicate parsed module parameter: " + name);
}

void validateModuleDescription(const ModuleDescription& description) {
    requireName(description.name, "module");
    std::set<std::string> names;
    for (const auto& field : description.fields) {
        requireName(field.name, "field");
        if (!names.insert(field.name).second)
            throw std::invalid_argument("duplicate module field: " + field.name);
        fieldTypeName(field.type);
        fieldRoleName(field.role);
        initializationName(field.initialization);
        if (field.initialization == FieldInitialization::HomogeneousCorrection &&
            std::none_of(description.fields.begin(), description.fields.end(),
                         [&](const FieldDescription& source) {
                             return source.name == field.source_field &&
                                    source.name != field.name &&
                                    source.type == DescribedFieldType::Scalar &&
                                    source.location == DescribedFieldLocation::Cell;
                         }))
            throw std::invalid_argument("correction needs a declared scalar source field: " +
                                        field.name);
        if (description.supports_case_template &&
            field.initialization == FieldInitialization::File &&
            field.location != DescribedFieldLocation::Cell)
            throw std::invalid_argument("generic case templates require cell input fields: " +
                                        field.name);
        if (field.template_value.find_first_of("\r\n{}#") != std::string::npos)
            throw std::invalid_argument("invalid field tutorial literal: " + field.name);
        if (description.supports_case_template && field.initialization == FieldInitialization::File)
            checkTutorialLiteral(field);
    }
    std::set<std::string> parameters;
    for (const auto& parameter : description.parameters) {
        requireName(parameter.name, "parameter");
        const std::string key =
            std::string(dictionaryName(parameter.dictionary)) + "." + parameter.name;
        if (!parameters.insert(key).second)
            throw std::invalid_argument("duplicate module parameter: " + key);
        typeName(parameter.type);
        if (parameter.required && (parameter.default_value || parameter.default_factory))
            throw std::invalid_argument("required module parameter cannot have a default: " + key);
        if (parameter.default_value && parameter.default_factory)
            throw std::invalid_argument("module parameter has two defaults: " + key);
        if ((parameter.minimum && !std::isfinite(*parameter.minimum)) ||
            (parameter.maximum && !std::isfinite(*parameter.maximum)) ||
            (parameter.minimum && parameter.maximum &&
             (*parameter.minimum > *parameter.maximum ||
              (*parameter.minimum == *parameter.maximum &&
               (parameter.exclusive_minimum || parameter.exclusive_maximum)))))
            throw std::invalid_argument("invalid module parameter range: " + key);
        if ((parameter.minimum || parameter.maximum) && parameter.type != ParameterType::Number &&
            parameter.type != ParameterType::Integer)
            throw std::invalid_argument("range requires a numeric module parameter: " + key);
        if (!parameter.choices.empty() && parameter.type != ParameterType::Word)
            throw std::invalid_argument("choices require a word module parameter: " + key);
        if (parameter.default_value)
            checkValue(parameter, *parameter.default_value, key + " default");
        if (parameter.template_value)
            checkValue(parameter, *parameter.template_value, key + " tutorial");
        if (description.supports_case_template && parameter.required && !parameter.template_value)
            throw std::invalid_argument("required module parameter needs a tutorial value: " + key);
    }
    std::set<std::string> equations;
    for (const auto& equation : description.equations) {
        requireName(equation.name, "equation");
        if (!equations.insert(equation.name).second)
            throw std::invalid_argument("duplicate module equation: " + equation.name);
        if (!names.count(equation.unknown))
            throw std::invalid_argument("module equation has an undeclared unknown: " +
                                        equation.unknown);
        std::set<std::string> terms;
        for (const auto& term : equation.terms) {
            requireName(term, "term");
            if (!terms.insert(term).second)
                throw std::invalid_argument("duplicate module term: " + term);
        }
    }
    for (const auto method : description.time_methods)
        timeName(method);
    if (description.supports_case_template && description.time_methods.empty())
        throw std::invalid_argument("module tutorial requires a time method");
    if (description.supports_case_template &&
        std::find(description.boundary_types.begin(), description.boundary_types.end(),
                  "zeroGradient") == description.boundary_types.end())
        throw std::invalid_argument(
            "generic module tutorial requires zeroGradient boundary capability");
    for (const auto& constraint : description.parameter_constraints) {
        relationName(constraint.relation);
        for (const auto& name : {constraint.left, constraint.right}) {
            const auto found = std::find_if(
                description.parameters.begin(), description.parameters.end(),
                [&](const ParameterDescription& parameter) {
                    return parameter.dictionary == constraint.dictionary && parameter.name == name;
                });
            if (found == description.parameters.end() ||
                (found->type != ParameterType::Number && found->type != ParameterType::Integer))
                throw std::invalid_argument(
                    "module relation references an undeclared numeric parameter: " + name);
        }
    }
}

ModuleParameters readModuleParameters(const Case& problem, const ModuleDescription& description) {
    validateModuleDescription(description);
    ModuleParameters result;
    for (const auto& parameter : description.parameters) {
        const auto& input = parameter.dictionary == ParameterDictionary::Physics
                                ? problem.read.physics()
                                : problem.read.solution();
        auto& values = valuesFor(result, parameter.dictionary);
        std::optional<ParameterValue> value;
        if (input.contains(parameter.name))
            value = readValue(input, parameter);
        else if (parameter.default_value)
            value = parameter.default_value;
        else if (parameter.default_factory)
            value = parameter.default_factory(values);
        else if (parameter.required)
            value = readValue(input, parameter); // Existing file/line diagnostics.
        if (value) {
            const auto line = input.inspect(parameter.name).line;
            checkValue(parameter, *value,
                       input.sourcePath().string() + ":" + std::to_string(line) + ": " +
                           parameter.name);
            values.insert(parameter.name, *value);
        }
    }
    checkRelations(description, result, "module " + description.name);
    return result;
}

void validateModuleConfiguration(const Case& problem, const ModuleDescription& description) {
    (void)readModuleParameters(problem, description);
    if (!description.time_methods.empty() &&
        std::find(description.time_methods.begin(), description.time_methods.end(),
                  problem.read.methods().time) == description.time_methods.end())
        throw std::invalid_argument("module " + description.name + " does not support time " +
                                    timeName(problem.read.methods().time));
    if (!description.compatible_solvers.empty() &&
        std::find(description.compatible_solvers.begin(), description.compatible_solvers.end(),
                  problem.solver()) == description.compatible_solvers.end())
        throw std::invalid_argument("module " + description.name + " is incompatible with solver " +
                                    problem.solver());
    for (const auto& field : description.fields) {
        if (field.initialization == FieldInitialization::File)
            checkFieldHeader(problem.definition().fields_directory / (field.name + ".field"), field,
                             description.boundary_types);
    }
    for (const auto& equation : description.equations) {
        const std::string key = "equation." + equation.name;
        const auto found = problem.read.methods().named.find(key);
        if (found == problem.read.methods().named.end())
            throw std::invalid_argument("module " + description.name +
                                        " needs numerical configuration " + key);
        found->second.options.requireComplete(key);
        for (const auto& option :
             {"kspType", "pcType", "absoluteTolerance", "relativeTolerance", "maxIterations"})
            if (!problem.read.solution().contains(key + "." + option))
                throw std::invalid_argument(problem.read.solution().sourcePath().string() +
                                            ": missing entry " + key + "." + option);
    }
}

} // namespace babelsim
