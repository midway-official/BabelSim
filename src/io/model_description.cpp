#include "babelsim/model_description.h"
#include "babelsim/case.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace babelsim {
namespace {

const char* dictionaryName(ParameterDictionary dictionary) {
    return dictionary == ParameterDictionary::Physics ? "physics" : "solution";
}
const char* typeName(ParameterType type) {
    switch (type) {
    case ParameterType::Number: return "number";
    case ParameterType::Integer: return "integer";
    case ParameterType::Boolean: return "boolean";
    case ParameterType::Word: return "string";
    }
    throw std::invalid_argument("invalid parameter type");
}
const char* fieldTypeName(DescribedFieldType type) {
    switch (type) {
    case DescribedFieldType::Scalar: return "scalar";
    case DescribedFieldType::Vector: return "vector";
    case DescribedFieldType::Tensor: return "tensor";
    }
    throw std::invalid_argument("invalid field type");
}
const char* timeName(TimeMethod method) {
    switch (method) {
    case TimeMethod::Steady: return "steady";
    case TimeMethod::Euler: return "euler";
    case TimeMethod::BDF2: return "bdf2";
    }
    throw std::invalid_argument("invalid time method");
}
const char* relationName(ParameterRelation relation) {
    switch (relation) {
    case ParameterRelation::LessEqual: return "<=";
    case ParameterRelation::GreaterEqual: return ">=";
    case ParameterRelation::Equal: return "==";
    }
    throw std::invalid_argument("invalid parameter relation");
}
bool safeName(const std::string& name) {
    if (name.empty()) return false;
    return std::all_of(name.begin(), name.end(), [](unsigned char character) {
        return std::isalnum(character) || character == '_' || character == '.' || character == '-';
    }) && name != "." && name != "..";
}
void requireName(const std::string& name, const std::string& context) {
    if (!safeName(name)) throw std::invalid_argument("invalid " + context + " name: " + name);
}

std::string escaped(const std::string& text) {
    std::ostringstream result;
    for (unsigned char character : text) {
        switch (character) {
        case '"': result << "\\\""; break;
        case '\\': result << "\\\\"; break;
        case '\n': result << "\\n"; break;
        case '\r': result << "\\r"; break;
        case '\t': result << "\\t"; break;
        default:
            if (character < 0x20) result << "\\u" << std::hex << std::setw(4)
                << std::setfill('0') << static_cast<unsigned int>(character) << std::dec;
            else result << character;
        }
    }
    return result.str();
}
std::string quote(const std::string& text) { return "\"" + escaped(text) + "\""; }
std::string singleLine(std::string text) {
    std::replace(text.begin(), text.end(), '\n', ' ');
    std::replace(text.begin(), text.end(), '\r', ' ');
    return text;
}
std::string tableCell(const std::string& text) {
    std::string result;
    for (char character : text) {
        if (character == '|') result += "\\|";
        else if (character == '\n') result += "<br>";
        else if (character != '\r') result += character;
    }
    return result;
}
std::string valueText(const ParameterValue& value, bool json = false) {
    std::ostringstream text;
    text << std::setprecision(17);
    if (const auto* number = std::get_if<double>(&value)) text << *number;
    else if (const auto* integer = std::get_if<int>(&value)) text << *integer;
    else if (const auto* boolean = std::get_if<bool>(&value)) text << (*boolean ? "true" : "false");
    else text << (json ? quote(std::get<std::string>(value)) : std::get<std::string>(value));
    return text.str();
}
bool hasType(const ParameterValue& value, ParameterType type) {
    return (type == ParameterType::Number && std::holds_alternative<double>(value)) ||
           (type == ParameterType::Integer && std::holds_alternative<int>(value)) ||
           (type == ParameterType::Boolean && std::holds_alternative<bool>(value)) ||
           (type == ParameterType::Word && std::holds_alternative<std::string>(value));
}
double numericValue(const ParameterValue& value) {
    if (const auto* number = std::get_if<double>(&value)) return *number;
    if (const auto* integer = std::get_if<int>(&value)) return *integer;
    throw std::invalid_argument("parameter relation needs number/integer values");
}
void checkValue(const ParameterDescription& parameter, const ParameterValue& value,
                const std::string& context) {
    if (!hasType(value, parameter.type))
        throw std::invalid_argument(context + ": expected " + typeName(parameter.type));
    if (parameter.type == ParameterType::Number || parameter.type == ParameterType::Integer) {
        const double number = numericValue(value);
        if (!std::isfinite(number)) throw std::invalid_argument(context + ": expected a finite number");
        if (parameter.minimum && (parameter.exclusive_minimum
            ? number <= *parameter.minimum : number < *parameter.minimum))
            throw std::invalid_argument(context + ": value below declared minimum " +
                valueText(*parameter.minimum) + (parameter.exclusive_minimum ? " (exclusive)" : ""));
        if (parameter.maximum && (parameter.exclusive_maximum
            ? number >= *parameter.maximum : number > *parameter.maximum))
            throw std::invalid_argument(context + ": value above declared maximum " +
                valueText(*parameter.maximum) + (parameter.exclusive_maximum ? " (exclusive)" : ""));
    } else if (parameter.type == ParameterType::Word) {
        const auto& word = std::get<std::string>(value);
        if (word.empty() || std::any_of(word.begin(), word.end(), [](unsigned char character) {
            return std::isspace(character) || character == '#';
        })) throw std::invalid_argument(context + ": expected one nonempty word");
        if (!parameter.choices.empty() && std::find(parameter.choices.begin(),
            parameter.choices.end(), word) == parameter.choices.end())
            throw std::invalid_argument(context + ": value outside declared choices");
    }
}
ParameterValue readValue(const Parameters& input, const ParameterDescription& parameter) {
    switch (parameter.type) {
    case ParameterType::Number: return input.number(parameter.name);
    case ParameterType::Integer: return input.integer(parameter.name);
    case ParameterType::Boolean: return input.boolean(parameter.name);
    case ParameterType::Word: return input.word(parameter.name);
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
        if (!dictionary.contains(constraint.left) || !dictionary.contains(constraint.right)) continue;
        const double left = numericValue(dictionary.value(constraint.left));
        const double right = numericValue(dictionary.value(constraint.right));
        const bool valid = constraint.relation == ParameterRelation::LessEqual ? left <= right :
            (constraint.relation == ParameterRelation::GreaterEqual ? left >= right : left == right);
        if (!valid) throw std::invalid_argument(context + ": " + constraint.left + " " +
            relationName(constraint.relation) + " " + constraint.right + " required; " +
            constraint.explanation);
    }
}
void writeFile(const std::filesystem::path& path, const std::string& content) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot write " + path.string());
    output << content;
    output.close();
    if (!output) throw std::runtime_error("cannot finish writing " + path.string());
}
void jsonStrings(std::ostream& output, const std::vector<std::string>& values) {
    output << '[';
    for (std::size_t i = 0; i < values.size(); ++i) output << (i ? ", " : "") << quote(values[i]);
    output << ']';
}
void checkTutorialLiteral(const FieldDescription& field) {
    const int components = field.type == DescribedFieldType::Scalar ? 1 :
        (field.type == DescribedFieldType::Vector ? 3 : 9);
    std::istringstream input(field.template_value);
    for (int component = 0; component < components; ++component) {
        double value = 0.0;
        if (!(input >> value) || !std::isfinite(value))
            throw std::invalid_argument("field tutorial requires " + std::to_string(components) +
                " finite components: " + field.name);
    }
    std::string extra;
    if (input >> extra) throw std::invalid_argument("unexpected field tutorial component: " + field.name);
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
        if (tokens[i].text == "{") ++depth;
        else if (tokens[i].text == "}") --depth;
        else if (depth == 1 && i + 1 < tokens.size()) {
            if (tokens[i].text == "type") type = tokens[++i].text == fieldTypeName(field.type);
            else if (tokens[i].text == "location") location = tokens[++i].text == expectedLocation;
            else if (tokens[i].text == "internal" && tokens[i + 1].text == "file" && i + 2 < tokens.size())
                i += 2; // A data filename may itself be a header keyword.
        } else if (depth == 3 && tokens[i].text == "type" && i + 1 < tokens.size() && !boundaries.empty()) {
            const auto& token = tokens[++i];
            if (std::find(boundaries.begin(), boundaries.end(), token.text) == boundaries.end())
                throw std::invalid_argument(path.string() + ":" + std::to_string(token.line) +
                    ": module does not support boundary type " + token.text);
        }
    }
    if (tokens.size() < 8 || tokens[0].text != "field" || tokens[1].text != field.name ||
        tokens[2].text != "{" || !type || !location)
        throw std::invalid_argument(path.string() + ": input header differs from module field " + field.name);
}
} // namespace

const ModuleRegistration*& ModuleRegistration::first() noexcept {
    static const ModuleRegistration* head = nullptr;
    return head;
}
ModuleRegistration::ModuleRegistration(const ModuleDescription& description) noexcept
    : description_(&description), next_(first()) { first() = this; }
ModuleRegistration::~ModuleRegistration() noexcept {
    const ModuleRegistration** link = &first();
    while (*link && *link != this) link = &(*link)->next_;
    if (*link) *link = next_;
}
std::vector<const ModuleDescription*> registeredModuleDescriptions() {
    std::vector<const ModuleDescription*> result;
    std::set<std::string> names;
    for (auto* entry = ModuleRegistration::first(); entry; entry = entry->next_) {
        validateModuleDescription(*entry->description_);
        if (!names.insert(entry->description_->name).second)
            throw std::invalid_argument("duplicate module registration: " + entry->description_->name);
        result.push_back(entry->description_);
    }
    return result;
}

bool ParameterValues::contains(const std::string& name) const { return values_.count(name) != 0; }
const ParameterValue& ParameterValues::value(const std::string& name) const {
    const auto found = values_.find(name);
    if (found == values_.end()) throw std::invalid_argument("undeclared or unset module parameter: " + name);
    return found->second;
}
double ParameterValues::number(const std::string& name) const {
    const auto* result = std::get_if<double>(&value(name));
    if (!result) throw std::invalid_argument("module parameter is not a number: " + name);
    return *result;
}
int ParameterValues::integer(const std::string& name) const {
    const auto* result = std::get_if<int>(&value(name));
    if (!result) throw std::invalid_argument("module parameter is not an integer: " + name);
    return *result;
}
bool ParameterValues::boolean(const std::string& name) const {
    const auto* result = std::get_if<bool>(&value(name));
    if (!result) throw std::invalid_argument("module parameter is not a boolean: " + name);
    return *result;
}
const std::string& ParameterValues::word(const std::string& name) const {
    const auto* result = std::get_if<std::string>(&value(name));
    if (!result) throw std::invalid_argument("module parameter is not a word: " + name);
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
        if (!names.insert(field.name).second) throw std::invalid_argument("duplicate module field: " + field.name);
        fieldTypeName(field.type);
        if (description.supports_case_template && field.role == FieldRole::Input &&
            field.location != DescribedFieldLocation::Cell)
            throw std::invalid_argument("generic case templates require cell input fields: " + field.name);
        if (field.template_value.find_first_of("\r\n{}#") != std::string::npos)
            throw std::invalid_argument("invalid field tutorial literal: " + field.name);
        if (description.supports_case_template && field.role == FieldRole::Input) checkTutorialLiteral(field);
    }
    std::set<std::string> parameters;
    for (const auto& parameter : description.parameters) {
        requireName(parameter.name, "parameter");
        const std::string key = std::string(dictionaryName(parameter.dictionary)) + "." + parameter.name;
        if (!parameters.insert(key).second) throw std::invalid_argument("duplicate module parameter: " + key);
        typeName(parameter.type);
        if (parameter.required && (parameter.default_value || parameter.default_factory))
            throw std::invalid_argument("required module parameter cannot have a default: " + key);
        if (parameter.default_value && parameter.default_factory)
            throw std::invalid_argument("module parameter has two defaults: " + key);
        if ((parameter.minimum && !std::isfinite(*parameter.minimum)) ||
            (parameter.maximum && !std::isfinite(*parameter.maximum)) ||
            (parameter.minimum && parameter.maximum && (*parameter.minimum > *parameter.maximum ||
                (*parameter.minimum == *parameter.maximum &&
                    (parameter.exclusive_minimum || parameter.exclusive_maximum)))))
            throw std::invalid_argument("invalid module parameter range: " + key);
        if ((parameter.minimum || parameter.maximum) && parameter.type != ParameterType::Number &&
            parameter.type != ParameterType::Integer)
            throw std::invalid_argument("range requires a numeric module parameter: " + key);
        if (!parameter.choices.empty() && parameter.type != ParameterType::Word)
            throw std::invalid_argument("choices require a word module parameter: " + key);
        if (parameter.default_value) checkValue(parameter, *parameter.default_value, key + " default");
        if (parameter.template_value) checkValue(parameter, *parameter.template_value, key + " tutorial");
        if (description.supports_case_template && parameter.required && !parameter.template_value)
            throw std::invalid_argument("required module parameter needs a tutorial value: " + key);
    }
    std::set<std::string> equations;
    for (const auto& equation : description.equations) {
        requireName(equation.name, "equation");
        if (!equations.insert(equation.name).second)
            throw std::invalid_argument("duplicate module equation: " + equation.name);
        if (!names.count(equation.unknown))
            throw std::invalid_argument("module equation has an undeclared unknown: " + equation.unknown);
        std::set<std::string> terms;
        for (const auto& term : equation.terms) {
            requireName(term, "term");
            if (!terms.insert(term).second) throw std::invalid_argument("duplicate module term: " + term);
        }
    }
    for (const auto method : description.time_methods) timeName(method);
    if (description.supports_case_template && description.time_methods.empty())
        throw std::invalid_argument("module tutorial requires a time method");
    if (description.supports_case_template && std::find(description.boundary_types.begin(),
        description.boundary_types.end(), "zeroGradient") == description.boundary_types.end())
        throw std::invalid_argument("generic module tutorial requires zeroGradient boundary capability");
    for (const auto& constraint : description.parameter_constraints) {
        relationName(constraint.relation);
        for (const auto& name : {constraint.left, constraint.right}) {
            const auto found = std::find_if(description.parameters.begin(), description.parameters.end(),
                [&](const ParameterDescription& parameter) {
                    return parameter.dictionary == constraint.dictionary && parameter.name == name;
                });
            if (found == description.parameters.end() ||
                (found->type != ParameterType::Number && found->type != ParameterType::Integer))
                throw std::invalid_argument("module relation references an undeclared numeric parameter: " + name);
        }
    }
}

ModuleParameters readModuleParameters(const Case& problem, const ModuleDescription& description) {
    validateModuleDescription(description);
    ModuleParameters result;
    for (const auto& parameter : description.parameters) {
        const auto& input = parameter.dictionary == ParameterDictionary::Physics ? problem.physics() : problem.solution();
        auto& values = valuesFor(result, parameter.dictionary);
        std::optional<ParameterValue> value;
        if (input.contains(parameter.name)) value = readValue(input, parameter);
        else if (parameter.default_value) value = parameter.default_value;
        else if (parameter.default_factory) value = parameter.default_factory(values);
        else if (parameter.required) value = readValue(input, parameter); // Existing file/line diagnostics.
        if (value) {
            const auto line = input.inspect(parameter.name).line;
            checkValue(parameter, *value, input.sourcePath().string() + ":" + std::to_string(line) + ": " + parameter.name);
            values.insert(parameter.name, *value);
        }
    }
    checkRelations(description, result, "module " + description.name);
    return result;
}

void validateModuleConfiguration(const Case& problem, const ModuleDescription& description) {
    (void)readModuleParameters(problem, description);
    if (!description.time_methods.empty() && std::find(description.time_methods.begin(),
        description.time_methods.end(), problem.methods().time) == description.time_methods.end())
        throw std::invalid_argument("module " + description.name + " does not support time " + timeName(problem.methods().time));
    if (!description.compatible_solvers.empty() && std::find(description.compatible_solvers.begin(),
        description.compatible_solvers.end(), problem.solver()) == description.compatible_solvers.end())
        throw std::invalid_argument("module " + description.name + " is incompatible with solver " + problem.solver());
    for (const auto& field : description.fields) {
        if (field.role == FieldRole::Input)
            checkFieldHeader(problem.definition().fields_directory / (field.name + ".field"), field,
                description.boundary_types);
    }
    for (const auto& equation : description.equations) {
        const std::string key = "equation." + equation.name;
        const auto found = problem.methods().named.find(key);
        if (found == problem.methods().named.end())
            throw std::invalid_argument("module " + description.name + " needs numerical configuration " + key);
        found->second.options.requireComplete(key);
        for (const auto& option : {"kspType", "pcType", "absoluteTolerance", "relativeTolerance", "maxIterations"})
            if (!problem.solution().contains(key + "." + option))
                throw std::invalid_argument(problem.solution().sourcePath().string() +
                    ": missing entry " + key + "." + option);
    }
}

std::string moduleMarkdown(const ModuleDescription& description) {
    validateModuleDescription(description);
    std::ostringstream output;
    output << "# " << description.name << "\n\n" << description.explanation << "\n\n"
        << "## Fields\n\n| Name | Type | Location | Role | Unit | Tutorial value | Description |\n"
        << "| --- | --- | --- | --- | --- | --- | --- |\n";
    for (const auto& field : description.fields) output << "| " << field.name << " | "
        << fieldTypeName(field.type) << " | " << (field.location == DescribedFieldLocation::Cell ? "cell" : "face")
        << " | " << (field.role == FieldRole::Input ? "input" : "derived") << " | "
        << tableCell(field.unit) << " | " << (field.role == FieldRole::Input ? field.template_value : "program-owned")
        << " | " << tableCell(field.explanation) << " |\n";
    output << "\n## Parameters\n\n| Dictionary | Name | Type | Unit | Default | Range / choices | Description |\n"
        << "| --- | --- | --- | --- | --- | --- | --- |\n";
    for (const auto& parameter : description.parameters) {
        output << "| " << dictionaryName(parameter.dictionary) << " | " << parameter.name
            << " | " << typeName(parameter.type) << " | " << tableCell(parameter.unit) << " | ";
        if (parameter.default_value) output << tableCell(valueText(*parameter.default_value));
        else if (parameter.default_factory) output << tableCell(parameter.default_explanation);
        else output << (parameter.required ? "required" : "solver-controlled");
        output << " | ";
        if (parameter.minimum) output << (parameter.exclusive_minimum ? "> " : ">= ") << valueText(*parameter.minimum);
        if (parameter.minimum && parameter.maximum) output << "; ";
        if (parameter.maximum) output << (parameter.exclusive_maximum ? "< " : "<= ") << valueText(*parameter.maximum);
        for (std::size_t i = 0; i < parameter.choices.size(); ++i) output << (i ? ", " : "") << tableCell(parameter.choices[i]);
        output << " | " << tableCell(parameter.explanation) << " |\n";
    }
    output << "\n## Equations\n\n";
    for (const auto& equation : description.equations) {
        output << "- `" << equation.name << "`: unknown `" << equation.unknown << "`; terms ";
        for (std::size_t i = 0; i < equation.terms.size(); ++i) output << (i ? ", " : "") << '`' << equation.terms[i] << '`';
        output << ". " << equation.explanation << " Numerical options: ";
        for (std::size_t i = 0; i < equation.numerical_options.size(); ++i)
            output << (i ? ", " : "") << equation.numerical_options[i];
        output << ".\n";
    }
    output << "\nSupported time methods: ";
    for (std::size_t i = 0; i < description.time_methods.size(); ++i)
        output << (i ? ", " : "") << timeName(description.time_methods[i]);
    output << ".\n\nBoundary capabilities: ";
    for (std::size_t i = 0; i < description.boundary_types.size(); ++i)
        output << (i ? ", " : "") << description.boundary_types[i];
    output << ".\n";
    if (!description.compatible_solvers.empty()) {
        output << "\nCompatible solvers: ";
        for (std::size_t i = 0; i < description.compatible_solvers.size(); ++i)
            output << (i ? ", " : "") << description.compatible_solvers[i];
        output << ".\n";
    }
    for (const auto& constraint : description.parameter_constraints)
        output << "\n- Required relation: `" << constraint.left << " " << relationName(constraint.relation)
            << " " << constraint.right << "`. " << constraint.explanation << '\n';
    for (const auto& limitation : description.limitations) output << "\n- " << limitation << '\n';
    output << "\nCase template: " << (description.supports_case_template ? "available" : "not supplied") << ".\n";
    return output.str();
}

std::string moduleJson(const ModuleDescription& description) {
    validateModuleDescription(description);
    std::ostringstream output;
    output << std::setprecision(17)
        << "{\n  \"$schema\": \"https://json-schema.org/draft/2020-12/schema\",\n"
        << "  \"title\": " << quote(description.name) << ",\n  \"description\": " << quote(description.explanation)
        << ",\n  \"type\": \"object\",\n  \"properties\": {\n";
    for (const auto dictionary : {ParameterDictionary::Physics, ParameterDictionary::Solution}) {
        output << "    " << quote(dictionaryName(dictionary)) << ": {\"type\": \"object\", \"properties\": {";
        bool first = true;
        std::vector<std::string> required;
        for (const auto& parameter : description.parameters) {
            if (parameter.dictionary != dictionary) continue;
            output << (first ? "\n" : ",\n") << "      " << quote(parameter.name) << ": {\"type\": "
                << quote(typeName(parameter.type)) << ", \"description\": " << quote(parameter.explanation)
                << ", \"x-unit\": " << quote(parameter.unit);
            first = false;
            if (parameter.required) required.push_back(parameter.name);
            if (parameter.default_value) output << ", \"default\": " << valueText(*parameter.default_value, true);
            if (parameter.default_factory) output << ", \"x-computedDefault\": " << quote(parameter.default_explanation);
            if (parameter.template_value) output << ", \"examples\": [" << valueText(*parameter.template_value, true) << ']';
            if (parameter.minimum) output << ", " << quote(parameter.exclusive_minimum ? "exclusiveMinimum" : "minimum")
                << ": " << *parameter.minimum;
            if (parameter.maximum) output << ", " << quote(parameter.exclusive_maximum ? "exclusiveMaximum" : "maximum")
                << ": " << *parameter.maximum;
            if (!parameter.choices.empty()) { output << ", \"enum\": "; jsonStrings(output, parameter.choices); }
            output << '}';
        }
        output << "\n    }, \"required\": "; jsonStrings(output, required);
        // Descriptions compose; equation controls and other models use the same
        // dictionaries. Actual unknown/unused keys are checked by Case::validate.
        output << ", \"additionalProperties\": true}" << (dictionary == ParameterDictionary::Physics ? ",\n" : "\n");
    }
    std::vector<std::string> requiredDictionaries;
    for (const auto dictionary : {ParameterDictionary::Physics, ParameterDictionary::Solution})
        if (std::any_of(description.parameters.begin(), description.parameters.end(),
            [dictionary](const ParameterDescription& parameter) {
                return parameter.dictionary == dictionary && parameter.required;
            })) requiredDictionaries.emplace_back(dictionaryName(dictionary));
    output << "  },\n  \"required\": "; jsonStrings(output, requiredDictionaries);
    output << ",\n  \"x-fields\": [";
    for (std::size_t i = 0; i < description.fields.size(); ++i) {
        const auto& field = description.fields[i];
        output << (i ? ",\n" : "\n") << "    {\"name\": " << quote(field.name)
            << ", \"type\": " << quote(fieldTypeName(field.type)) << ", \"location\": "
            << quote(field.location == DescribedFieldLocation::Cell ? "cell" : "face")
            << ", \"role\": " << quote(field.role == FieldRole::Input ? "input" : "derived")
            << ", \"unit\": " << quote(field.unit) << ", \"tutorialValue\": " << quote(field.template_value)
            << ", \"description\": " << quote(field.explanation) << '}';
    }
    output << "\n  ],\n  \"x-equations\": [";
    for (std::size_t i = 0; i < description.equations.size(); ++i) {
        const auto& equation = description.equations[i];
        output << (i ? ",\n" : "\n") << "    {\"name\": " << quote(equation.name)
            << ", \"unknown\": " << quote(equation.unknown) << ", \"terms\": ";
        jsonStrings(output, equation.terms);
        output << ", \"numericalOptions\": "; jsonStrings(output, equation.numerical_options);
        output << ", \"description\": " << quote(equation.explanation) << '}';
    }
    std::vector<std::string> times;
    for (const auto method : description.time_methods) times.emplace_back(timeName(method));
    output << "\n  ],\n  \"x-timeMethods\": "; jsonStrings(output, times);
    output << ",\n  \"x-boundaryTypes\": "; jsonStrings(output, description.boundary_types);
    output << ",\n  \"x-compatibleSolvers\": "; jsonStrings(output, description.compatible_solvers);
    output << ",\n  \"x-parameterConstraints\": [";
    for (std::size_t i = 0; i < description.parameter_constraints.size(); ++i) {
        const auto& constraint = description.parameter_constraints[i];
        output << (i ? ", " : "") << "{\"dictionary\": " << quote(dictionaryName(constraint.dictionary))
            << ", \"left\": " << quote(constraint.left) << ", \"relation\": " << quote(relationName(constraint.relation))
            << ", \"right\": " << quote(constraint.right) << ", \"description\": " << quote(constraint.explanation) << '}';
    }
    output << "],\n  \"x-limitations\": "; jsonStrings(output, description.limitations);
    output << ",\n  \"x-caseTemplate\": " << (description.supports_case_template ? "true" : "false") << "\n}\n";
    return output.str();
}

void writeModuleDescription(const ModuleDescription& description, ModuleOutputFormat format) {
    std::cout << (format == ModuleOutputFormat::Json ? moduleJson(description) : moduleMarkdown(description));
}
void writeModuleList(const std::vector<ModuleListing>& modules) {
    auto sorted = modules;
    std::sort(sorted.begin(), sorted.end(), [](const ModuleListing& left, const ModuleListing& right) {
        return left.name < right.name;
    });
    for (const auto& module : sorted) {
        std::cout << module.name << (module.solver ? " (solver; description " : " (model; description ")
            << (module.described ? "available" : "not supplied") << ")\n";
    }
}

void writeCaseTemplate(const ModuleDescription& description, const std::filesystem::path& directory) {
    validateModuleDescription(description);
    if (!description.supports_case_template)
        throw std::invalid_argument("module " + description.name + " has no generic case template");
    if (directory.empty() || std::filesystem::exists(directory) || std::filesystem::is_symlink(directory))
        throw std::invalid_argument("case template output must be a new directory: " + directory.string());
    if (!directory.parent_path().empty()) std::filesystem::create_directories(directory.parent_path());
    if (!std::filesystem::create_directory(directory))
        throw std::invalid_argument("case template output already exists: " + directory.string());
    try {
        for (const auto* child : {"mesh", "fields/initial", "physics", "numerics"})
            std::filesystem::create_directories(directory / child);
        writeFile(directory / "case.bs", "solver " + description.name +
            "\nmesh mesh/domain.mesh\nfields fields/initial\nphysics physics/model.bs\n"
            "methods numerics/methods.bs\nsolution numerics/solution.bs\ncontrol control.bs\noutput output.bs\n");
        writeFile(directory / "mesh/domain.mesh",
            "BABELSIM_MESH 3\nvertices 8\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n0 0 1\n1 0 1\n1 1 1\n0 1 1\n"
            "faces 6\nface 4 0 -1 0 0 4 7 3\nface 4 0 -1 1 1 2 6 5\nface 4 0 -1 2 0 1 5 4\n"
            "face 4 0 -1 3 3 7 6 2\nface 4 0 -1 4 0 3 2 1\nface 4 0 -1 5 4 5 6 7\n"
            "patches 6\npatch xmin generic\npatch xmax generic\npatch ymin generic\npatch ymax generic\n"
            "patch zmin generic\npatch zmax generic\nend\n");
        std::ostringstream physics, solution, methods;
        ModuleParameters tutorial;
        for (const auto& parameter : description.parameters) {
            auto& values = valuesFor(tutorial, parameter.dictionary);
            std::optional<ParameterValue> value = parameter.template_value ? parameter.template_value : parameter.default_value;
            if (!value && parameter.default_factory) value = parameter.default_factory(values);
            if (!value) continue;
            checkValue(parameter, *value, "tutorial " + parameter.name);
            values.insert(parameter.name, *value);
            auto& output = parameter.dictionary == ParameterDictionary::Physics ? physics : solution;
            output << "# " << singleLine(parameter.explanation + " [" + parameter.unit + "]") << '\n'
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
            if (field.role != FieldRole::Input) continue;
            std::ostringstream output;
            output << "field " << field.name << "\n{\n    type " << fieldTypeName(field.type)
                << "\n    location cell\n    internal uniform (" << field.template_value << ")\n    boundary\n    {\n";
            for (const auto* patch : {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"})
                output << "        " << patch << " { type zeroGradient }\n";
            output << "    }\n}\n";
            writeFile(directory / "fields/initial" / (field.name + ".field"), output.str());
        }
        writeFile(directory / "README.md", moduleMarkdown(description) +
            "\nThis unit-cube tutorial tests configuration and development workflow. "
            "It is not a physical validation case. Run with `babelsim-solve -case <this-directory>`.\n");
        writeFile(directory / "module.schema.json", moduleJson(description));
    } catch (...) {
        std::filesystem::remove_all(directory);
        throw;
    }
}
} // namespace babelsim
