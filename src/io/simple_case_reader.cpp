#include "babelsim/simple_io.h"

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

int integer(const std::filesystem::path& path, const ConfigLine& line) {
    if (line.tokens.size() != 2) invalid(path, line, "expected one value");
    try {
        std::size_t consumed = 0;
        const int value = std::stoi(line.tokens[1], &consumed);
        if (consumed == line.tokens[1].size()) return value;
    } catch (const std::exception&) {
    }
    invalid(path, line, "expected an integer");
}

}  // 匿名命名空间

SimpleCaseControl readSimpleCase(const CaseDefinition& definition) {
    SimpleCaseControl result;
    bool density = false;
    bool viscosity = false;
    for (const ConfigLine& line : readConfigLines(definition.physics_file)) {
        const std::string& key = line.tokens.front();
        if (key == "density" && !density) {
            result.fluid.density = number(definition.physics_file, line);
            density = true;
        } else if ((key == "dynamicViscosity" || key == "dynamic_viscosity") && !viscosity) {
            result.fluid.dynamic_viscosity = number(definition.physics_file, line);
            viscosity = true;
        } else {
            invalid(definition.physics_file, line, "unknown or duplicate incompressible property " + key);
        }
    }
    if (!density || !viscosity) {
        throw std::runtime_error("simple physics needs density and dynamicViscosity");
    }

    result.runtime.methods = readMethodsFile(definition.methods_file);
    result.runtime.time = readTimeControlFile(definition.control_file);
    if (result.runtime.methods.time != TimeMethod::Steady) {
        throw std::runtime_error("simple solver requires time steady in methods dictionary");
    }

    bool iterations = false;
    bool velocity_relaxation = false;
    bool pressure_relaxation = false;
    bool continuity_tolerance = false;
    bool velocity_tolerance = false;
    bool velocity_solver = false;
    bool pressure_solver = false;
    for (const ConfigLine& line : readConfigLines(definition.solution_file)) {
        const std::string& key = line.tokens.front();
        if (key == "maxIterations" || key == "max_iterations") {
            if (iterations) invalid(definition.solution_file, line, "duplicate maxIterations");
            result.simple.max_iterations = integer(definition.solution_file, line);
            iterations = true;
        } else if (key == "nonOrthogonalCorrections" || key == "non_orthogonal_corrections") {
            result.simple.non_orthogonal_corrections = integer(definition.solution_file, line);
        } else if (key == "velocityRelaxation" || key == "velocity_relaxation") {
            if (velocity_relaxation) invalid(definition.solution_file, line, "duplicate velocityRelaxation");
            result.simple.velocity_relaxation = number(definition.solution_file, line);
            velocity_relaxation = true;
        } else if (key == "pressureRelaxation" || key == "pressure_relaxation") {
            if (pressure_relaxation) invalid(definition.solution_file, line, "duplicate pressureRelaxation");
            result.simple.pressure_relaxation = number(definition.solution_file, line);
            pressure_relaxation = true;
        } else if (key == "continuityTolerance" || key == "continuity_tolerance") {
            if (continuity_tolerance) invalid(definition.solution_file, line, "duplicate continuityTolerance");
            result.simple.continuity_tolerance = number(definition.solution_file, line);
            continuity_tolerance = true;
        } else if (key == "velocityTolerance" || key == "velocity_tolerance") {
            if (velocity_tolerance) invalid(definition.solution_file, line, "duplicate velocityTolerance");
            result.simple.velocity_tolerance = number(definition.solution_file, line);
            velocity_tolerance = true;
        } else if (key == "velocitySolver" || key == "velocity_solver") {
            if (velocity_solver) invalid(definition.solution_file, line, "duplicate velocitySolver");
            readLinearSolverLine(definition.solution_file, line, result.simple.velocity_solver);
            velocity_solver = true;
        } else if (key == "pressureSolver" || key == "pressure_solver") {
            if (pressure_solver) invalid(definition.solution_file, line, "duplicate pressureSolver");
            readLinearSolverLine(definition.solution_file, line, result.simple.pressure_solver);
            pressure_solver = true;
        } else {
            invalid(definition.solution_file, line, "unknown SIMPLE solution setting " + key);
        }
    }
    if (!iterations || !velocity_relaxation || !pressure_relaxation || !continuity_tolerance ||
        !velocity_tolerance || !velocity_solver || !pressure_solver) {
        throw std::runtime_error("simple solution dictionary is incomplete");
    }
    result.fluid.validate();
    result.simple.validate();
    result.runtime.validate();
    return result;
}

}  // babelsim 命名空间
