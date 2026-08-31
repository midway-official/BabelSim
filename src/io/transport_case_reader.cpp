#include "babelsim/transport_io.h"

#include "babelsim/config.h"
#include "babelsim/numerics_io.h"

#include <cmath>
#include <stdexcept>

namespace babelsim {
namespace {

[[noreturn]] void invalid(
    const std::filesystem::path& path, const ConfigLine& line, const std::string& message)
{
    throw std::runtime_error("invalid " + path.string() + ":" +
                             std::to_string(line.number) + ": " + message);
}

double number(const std::filesystem::path& path, const ConfigLine& line) {
    if (line.tokens.size() != 2) invalid(path, line, "expected one value");
    try {
        std::size_t consumed = 0;
        const double value = std::stod(line.tokens[1], &consumed);
        if (consumed == line.tokens[1].size() && std::isfinite(value)) return value;
    } catch (const std::exception&) {
    }
    invalid(path, line, "expected a finite number");
}

}  // 匿名命名空间

TransportCaseControl readTransportCase(const CaseDefinition& definition) {
    TransportCaseControl result;
    bool storage = false;
    bool diffusivity = false;
    bool source = false;
    for (const ConfigLine& line : readConfigLines(definition.physics_file)) {
        const std::string& key = line.tokens.front();
        if (key == "storage" && !storage) {
            result.storage = number(definition.physics_file, line);
            storage = true;
        } else if (key == "diffusivity" && !diffusivity) {
            result.diffusivity = number(definition.physics_file, line);
            diffusivity = true;
        } else if (key == "source" && !source) {
            result.source = number(definition.physics_file, line);
            source = true;
        } else {
            invalid(definition.physics_file, line, "unknown or duplicate transport property " + key);
        }
    }
    if (!storage || !diffusivity || !source || !(result.storage > 0.0) ||
        !(result.diffusivity >= 0.0)) {
        throw std::runtime_error("transport physics needs positive storage, nonnegative diffusivity, and source");
    }

    result.runtime.methods = readMethodsFile(definition.methods_file);
    result.runtime.time = readTimeControlFile(definition.control_file);
    bool scalar_solver = false;
    for (const ConfigLine& line : readConfigLines(definition.solution_file)) {
        if (line.tokens.front() != "scalarSolver" || scalar_solver) {
            invalid(definition.solution_file, line, "transport solution only accepts one scalarSolver");
        }
        readLinearSolverLine(definition.solution_file, line, result.runtime.scalar_solver);
        scalar_solver = true;
    }
    if (!scalar_solver || result.runtime.methods.time == TimeMethod::Steady) {
        throw std::runtime_error("transport solution needs scalarSolver and a transient time method");
    }
    result.runtime.validate();
    return result;
}

}  // babelsim 命名空间
