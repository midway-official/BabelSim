#include "babelsim/case.h"

#include "babelsim/config.h"

#include <stdexcept>
#include <string>
#include <utility>

namespace babelsim {
namespace {

[[noreturn]] void invalid(
    const std::filesystem::path& path,
    const ConfigLine& line,
    const std::string& message)
{
    throw std::runtime_error(
        "invalid " + path.string() + ":" + std::to_string(line.number) +
        ": " + message);
}

void setPath(
    const std::filesystem::path& path,
    const ConfigLine& line,
    std::filesystem::path& target,
    const std::filesystem::path& root)
{
    if (line.tokens.size() != 2 || !target.empty()) {
        invalid(path, line, "expected one unique value for " + line.tokens.front());
    }
    target = root / line.tokens[1];
}

}  // 匿名命名空间

CaseDefinition readCase(const std::filesystem::path& case_directory) {
    CaseDefinition result;
    result.root = std::filesystem::absolute(case_directory);
    const auto path = result.root / "case.bs";
    if (!std::filesystem::is_directory(result.root)) {
        throw std::runtime_error("BabelSim case directory does not exist: " +
                                 result.root.string());
    }

    for (const ConfigLine& line : readConfigLines(path)) {
        const std::string& key = line.tokens.front();
        if (key == "solver") {
            if (line.tokens.size() != 2 || !result.solver.empty()) {
                invalid(path, line, "expected one unique solver name");
            }
            result.solver = line.tokens[1];
        } else if (key == "mesh") {
            setPath(path, line, result.mesh_file, result.root);
        } else if (key == "fields") {
            setPath(path, line, result.fields_directory, result.root);
        } else if (key == "physics") {
            setPath(path, line, result.physics_file, result.root);
        } else if (key == "numerics") {
            setPath(path, line, result.numerics_file, result.root);
        } else if (key == "output") {
            setPath(path, line, result.output_file, result.root);
        } else {
            invalid(path, line, "unknown entry " + key);
        }
    }
    if (result.solver.empty() || result.mesh_file.empty() ||
        result.fields_directory.empty() || result.physics_file.empty() ||
        result.numerics_file.empty() || result.output_file.empty()) {
        throw std::runtime_error("case.bs must define solver, mesh, fields, physics, numerics, and output");
    }
    return result;
}

OutputControl readOutputControl(const CaseDefinition& definition) {
    OutputControl result;
    bool has_directory = false;
    bool has_time = false;
    for (const ConfigLine& line : readConfigLines(definition.output_file)) {
        const std::string& key = line.tokens.front();
        if (line.tokens.size() != 2) {
            invalid(definition.output_file, line, "expected one value for " + key);
        }
        if (key == "directory") {
            if (has_directory) invalid(definition.output_file, line, "duplicate directory");
            result.directory = line.tokens[1];
            has_directory = true;
        } else if (key == "timeName" || key == "time_name") {
            if (has_time || line.tokens[1].empty()) {
                invalid(definition.output_file, line, "duplicate or empty time_name");
            }
            result.time_name = line.tokens[1];
            has_time = true;
        } else {
            invalid(definition.output_file, line, "unknown entry " + key);
        }
    }
    if (result.directory.is_absolute()) {
        throw std::runtime_error("output directory must be relative to the case directory");
    }
    return result;
}

}  // babelsim 命名空间
