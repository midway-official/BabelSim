#include "babelsim/thermal_io.h"

#include "babelsim/config.h"
#include "babelsim/numerics_io.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace babelsim {
namespace {

[[noreturn]] void invalid(
    const std::filesystem::path& path,
    const ConfigLine& line,
    const std::string& message)
{
    throw std::runtime_error(
        "invalid " + path.string() + ":" + std::to_string(line.number) + ": " + message);
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

ThermalCaseControl readThermalCase(const CaseDefinition& definition) {
    ThermalCaseControl result;
    bool density = false;
    bool heat_capacity = false;
    bool conductivity = false;
    bool source = false;
    for (const ConfigLine& line : readConfigLines(definition.physics_file)) {
        const std::string& key = line.tokens.front();
        if (key == "density" && !density) {
            result.material.density = number(definition.physics_file, line);
            density = true;
        } else if ((key == "heatCapacity" || key == "heat_capacity") && !heat_capacity) {
            result.material.heat_capacity = number(definition.physics_file, line);
            heat_capacity = true;
        } else if ((key == "conductivity" || key == "thermalConductivity") && !conductivity) {
            result.material.conductivity = number(definition.physics_file, line);
            conductivity = true;
        } else if ((key == "source" || key == "volumetricSource") && !source) {
            result.volumetric_source = number(definition.physics_file, line);
            source = true;
        } else {
            invalid(definition.physics_file, line, "unknown or duplicate thermal property " + key);
        }
    }
    if (!density || !heat_capacity || !conductivity || !source) {
        throw std::runtime_error("thermal physics needs density, heatCapacity, conductivity, and source");
    }

    result.runtime.methods = readMethodsFile(definition.methods_file);
    result.runtime.time = readTimeControlFile(definition.control_file);
    bool scalar_solver = false;
    for (const ConfigLine& line : readConfigLines(definition.solution_file)) {
        if (line.tokens.front() != "scalarSolver" || scalar_solver) {
            invalid(definition.solution_file, line, "thermal solution only accepts one scalarSolver");
        }
        readLinearSolverLine(definition.solution_file, line, result.runtime.scalar_solver);
        scalar_solver = true;
    }
    if (!scalar_solver || result.runtime.methods.time == TimeMethod::Steady) {
        throw std::runtime_error("thermal solution needs scalarSolver and a transient time method");
    }
    result.material.validate();
    result.runtime.validate();
    return result;
}

}  // babelsim 命名空间
