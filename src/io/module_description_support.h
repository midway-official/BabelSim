#pragma once

#include "babelsim/model_description.h"

namespace babelsim::detail::module_description {

// Private vocabulary and value checks shared by the schema validator, text
// renderers and tutorial writer. They do not read a Case or mutate a registry.
const char* dictionaryName(ParameterDictionary dictionary);
const char* typeName(ParameterType type);
const char* fieldRoleName(FieldRole role);
const char* initializationName(FieldInitialization source);
const char* fieldTypeName(DescribedFieldType type);
const char* timeName(TimeMethod method);
const char* relationName(ParameterRelation relation);

std::string quote(const std::string& text);
std::string valueText(const ParameterValue& value, bool json = false);

void checkValue(const ParameterDescription& parameter, const ParameterValue& value,
                const std::string& context);
ParameterValues& valuesFor(ModuleParameters& values, ParameterDictionary dictionary);
const ParameterValues& valuesFor(const ModuleParameters& values, ParameterDictionary dictionary);
void checkRelations(const ModuleDescription& description, const ModuleParameters& values,
                    const std::string& context);

} // namespace babelsim::detail::module_description
