#include "babelsim/model_description.h"
#include "module_description_support.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace babelsim {
namespace detail::module_description {

std::string escaped(const std::string& text) {
    std::ostringstream result;
    for (unsigned char character : text) {
        switch (character) {
        case '"':
            result << "\\\"";
            break;
        case '\\':
            result << "\\\\";
            break;
        case '\n':
            result << "\\n";
            break;
        case '\r':
            result << "\\r";
            break;
        case '\t':
            result << "\\t";
            break;
        default:
            if (character < 0x20)
                result << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                       << static_cast<unsigned int>(character) << std::dec;
            else
                result << character;
        }
    }
    return result.str();
}
std::string quote(const std::string& text) {
    return "\"" + escaped(text) + "\"";
}
std::string tableCell(const std::string& text) {
    std::string result;
    for (char character : text) {
        if (character == '|')
            result += "\\|";
        else if (character == '\n')
            result += "<br>";
        else if (character != '\r')
            result += character;
    }
    return result;
}
std::string valueText(const ParameterValue& value, bool json) {
    std::ostringstream text;
    text << std::setprecision(17);
    if (const auto* number = std::get_if<double>(&value))
        text << *number;
    else if (const auto* integer = std::get_if<int>(&value))
        text << *integer;
    else if (const auto* boolean = std::get_if<bool>(&value))
        text << (*boolean ? "true" : "false");
    else
        text << (json ? quote(std::get<std::string>(value)) : std::get<std::string>(value));
    return text.str();
}
void jsonStrings(std::ostream& output, const std::vector<std::string>& values) {
    output << '[';
    for (std::size_t i = 0; i < values.size(); ++i)
        output << (i ? ", " : "") << quote(values[i]);
    output << ']';
}

} // namespace detail::module_description

using namespace detail::module_description;

std::string moduleMarkdown(const ModuleDescription& description) {
    validateModuleDescription(description);
    std::ostringstream output;
    output << "# " << description.name << "\n\n"
           << description.explanation << "\n\n"
           << "## Fields\n\n| Name | Type | Location | Role | Initialization | Unit | Tutorial "
              "value | Description |\n"
           << "| --- | --- | --- | --- | --- | --- | --- | --- |\n";
    for (const auto& field : description.fields)
        output << "| " << field.name << " | " << fieldTypeName(field.type) << " | "
               << (field.location == DescribedFieldLocation::Cell ? "cell" : "face") << " | "
               << fieldRoleName(field.role) << " | " << initializationName(field.initialization)
               << " | " << tableCell(field.unit) << " | "
               << (field.initialization == FieldInitialization::File ? field.template_value
                                                                     : "program-owned")
               << " | " << tableCell(field.explanation) << " |\n";
    output << "\n## Parameters\n\n| Dictionary | Name | Type | Unit | Default | Range / choices | "
              "Description |\n"
           << "| --- | --- | --- | --- | --- | --- | --- |\n";
    for (const auto& parameter : description.parameters) {
        output << "| " << dictionaryName(parameter.dictionary) << " | " << parameter.name << " | "
               << typeName(parameter.type) << " | " << tableCell(parameter.unit) << " | ";
        if (parameter.default_value)
            output << tableCell(valueText(*parameter.default_value));
        else if (parameter.default_factory)
            output << tableCell(parameter.default_explanation);
        else
            output << (parameter.required ? "required" : "solver-controlled");
        output << " | ";
        if (parameter.minimum)
            output << (parameter.exclusive_minimum ? "> " : ">= ") << valueText(*parameter.minimum);
        if (parameter.minimum && parameter.maximum)
            output << "; ";
        if (parameter.maximum)
            output << (parameter.exclusive_maximum ? "< " : "<= ") << valueText(*parameter.maximum);
        for (std::size_t i = 0; i < parameter.choices.size(); ++i)
            output << (i ? ", " : "") << tableCell(parameter.choices[i]);
        output << " | " << tableCell(parameter.explanation) << " |\n";
    }
    output << "\n## Equations\n\n";
    for (const auto& equation : description.equations) {
        output << "- `" << equation.name << "`: unknown `" << equation.unknown << "`; terms ";
        for (std::size_t i = 0; i < equation.terms.size(); ++i)
            output << (i ? ", " : "") << '`' << equation.terms[i] << '`';
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
        output << "\n- Required relation: `" << constraint.left << " "
               << relationName(constraint.relation) << " " << constraint.right << "`. "
               << constraint.explanation << '\n';
    for (const auto& limitation : description.limitations)
        output << "\n- " << limitation << '\n';
    output << "\nCase template: "
           << (description.supports_case_template ? "available" : "not supplied") << ".\n";
    return output.str();
}

std::string moduleJson(const ModuleDescription& description) {
    validateModuleDescription(description);
    std::ostringstream output;
    output << std::setprecision(17)
           << "{\n  \"$schema\": \"https://json-schema.org/draft/2020-12/schema\",\n"
           << "  \"title\": " << quote(description.name)
           << ",\n  \"description\": " << quote(description.explanation)
           << ",\n  \"type\": \"object\",\n  \"properties\": {\n";
    for (const auto dictionary : {ParameterDictionary::Physics, ParameterDictionary::Solution}) {
        output << "    " << quote(dictionaryName(dictionary))
               << ": {\"type\": \"object\", \"properties\": {";
        bool first = true;
        std::vector<std::string> required;
        for (const auto& parameter : description.parameters) {
            if (parameter.dictionary != dictionary)
                continue;
            output << (first ? "\n" : ",\n") << "      " << quote(parameter.name)
                   << ": {\"type\": " << quote(typeName(parameter.type))
                   << ", \"description\": " << quote(parameter.explanation)
                   << ", \"x-unit\": " << quote(parameter.unit);
            first = false;
            if (parameter.required)
                required.push_back(parameter.name);
            if (parameter.default_value)
                output << ", \"default\": " << valueText(*parameter.default_value, true);
            if (parameter.default_factory)
                output << ", \"x-computedDefault\": " << quote(parameter.default_explanation);
            if (parameter.template_value)
                output << ", \"examples\": [" << valueText(*parameter.template_value, true) << ']';
            if (parameter.minimum)
                output << ", "
                       << quote(parameter.exclusive_minimum ? "exclusiveMinimum" : "minimum")
                       << ": " << *parameter.minimum;
            if (parameter.maximum)
                output << ", "
                       << quote(parameter.exclusive_maximum ? "exclusiveMaximum" : "maximum")
                       << ": " << *parameter.maximum;
            if (!parameter.choices.empty()) {
                output << ", \"enum\": ";
                jsonStrings(output, parameter.choices);
            }
            output << '}';
        }
        output << "\n    }, \"required\": ";
        jsonStrings(output, required);
        // Descriptions compose; equation controls and other models use the same
        // dictionaries. Actual unknown/unused keys are checked by Case::finalizeInitialization().
        output << ", \"additionalProperties\": true}"
               << (dictionary == ParameterDictionary::Physics ? ",\n" : "\n");
    }
    std::vector<std::string> requiredDictionaries;
    for (const auto dictionary : {ParameterDictionary::Physics, ParameterDictionary::Solution})
        if (std::any_of(description.parameters.begin(), description.parameters.end(),
                        [dictionary](const ParameterDescription& parameter) {
                            return parameter.dictionary == dictionary && parameter.required;
                        }))
            requiredDictionaries.emplace_back(dictionaryName(dictionary));
    output << "  },\n  \"required\": ";
    jsonStrings(output, requiredDictionaries);
    output << ",\n  \"x-fields\": [";
    for (std::size_t i = 0; i < description.fields.size(); ++i) {
        const auto& field = description.fields[i];
        output << (i ? ",\n" : "\n") << "    {\"name\": " << quote(field.name)
               << ", \"type\": " << quote(fieldTypeName(field.type)) << ", \"location\": "
               << quote(field.location == DescribedFieldLocation::Cell ? "cell" : "face")
               << ", \"role\": " << quote(fieldRoleName(field.role))
               << ", \"initialization\": " << quote(initializationName(field.initialization))
               << ", \"sourceField\": " << quote(field.source_field)
               << ", \"unit\": " << quote(field.unit)
               << ", \"tutorialValue\": " << quote(field.template_value)
               << ", \"description\": " << quote(field.explanation) << '}';
    }
    output << "\n  ],\n  \"x-equations\": [";
    for (std::size_t i = 0; i < description.equations.size(); ++i) {
        const auto& equation = description.equations[i];
        output << (i ? ",\n" : "\n") << "    {\"name\": " << quote(equation.name)
               << ", \"unknown\": " << quote(equation.unknown) << ", \"terms\": ";
        jsonStrings(output, equation.terms);
        output << ", \"numericalOptions\": ";
        jsonStrings(output, equation.numerical_options);
        output << ", \"description\": " << quote(equation.explanation) << '}';
    }
    std::vector<std::string> times;
    for (const auto method : description.time_methods)
        times.emplace_back(timeName(method));
    output << "\n  ],\n  \"x-timeMethods\": ";
    jsonStrings(output, times);
    output << ",\n  \"x-boundaryTypes\": ";
    jsonStrings(output, description.boundary_types);
    output << ",\n  \"x-compatibleSolvers\": ";
    jsonStrings(output, description.compatible_solvers);
    output << ",\n  \"x-parameterConstraints\": [";
    for (std::size_t i = 0; i < description.parameter_constraints.size(); ++i) {
        const auto& constraint = description.parameter_constraints[i];
        output << (i ? ", " : "")
               << "{\"dictionary\": " << quote(dictionaryName(constraint.dictionary))
               << ", \"left\": " << quote(constraint.left)
               << ", \"relation\": " << quote(relationName(constraint.relation))
               << ", \"right\": " << quote(constraint.right)
               << ", \"description\": " << quote(constraint.explanation) << '}';
    }
    output << "],\n  \"x-limitations\": ";
    jsonStrings(output, description.limitations);
    output << ",\n  \"x-caseTemplate\": " << (description.supports_case_template ? "true" : "false")
           << "\n}\n";
    return output.str();
}

void writeModuleDescription(const ModuleDescription& description, ModuleOutputFormat format) {
    std::cout << (format == ModuleOutputFormat::Json ? moduleJson(description)
                                                     : moduleMarkdown(description));
}
void writeModuleList(const std::vector<ModuleListing>& modules) {
    auto sorted = modules;
    std::sort(sorted.begin(), sorted.end(),
              [](const ModuleListing& left, const ModuleListing& right) {
                  return left.name < right.name;
              });
    for (const auto& module : sorted) {
        std::cout << module.name
                  << (module.solver ? " (solver; description " : " (model; description ")
                  << (module.described ? "available" : "not supplied") << ")\n";
    }
}

} // namespace babelsim
