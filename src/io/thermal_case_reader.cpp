#include "babelsim/thermal_io.h"

#include "babelsim/config.h"

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

double number(const std::filesystem::path& path, const ConfigLine& line, std::size_t index) {
    try {
        std::size_t consumed = 0;
        const double value = std::stod(line.tokens.at(index), &consumed);
        if (consumed == line.tokens.at(index).size() && std::isfinite(value)) return value;
    } catch (const std::exception&) {
    }
    invalid(path, line, "expected a finite number");
}

int integer(const std::filesystem::path& path, const ConfigLine& line, std::size_t index) {
    try {
        std::size_t consumed = 0;
        const int value = std::stoi(line.tokens.at(index), &consumed);
        if (consumed == line.tokens.at(index).size()) return value;
    } catch (const std::exception&) {
    }
    invalid(path, line, "expected an integer");
}

void oneValue(const std::filesystem::path& path, const ConfigLine& line) {
    if (line.tokens.size() != 2) invalid(path, line, "expected one value for " + line.tokens.front());
}

void readLinear(
    const std::filesystem::path& path,
    const ConfigLine& line,
    LinearSolverConfig& result)
{
    if (line.tokens.size() != 6) {
        invalid(path, line, "linear solver needs method, preconditioner, tolerances, and iterations");
    }
    if (line.tokens[1] == "cg") result.solver = LinearSolverType::ConjugateGradient;
    else if (line.tokens[1] == "bicgstab") result.solver = LinearSolverType::BiCGSTAB;
    else invalid(path, line, "unknown linear solver " + line.tokens[1]);
    if (line.tokens[2] == "incompleteCholesky") {
        result.preconditioner = PreconditionerType::IncompleteCholesky;
    } else if (line.tokens[2] == "ilut") {
        result.preconditioner = PreconditionerType::ILUT;
    } else {
        invalid(path, line, "unknown preconditioner " + line.tokens[2]);
    }
    result.absolute_tolerance = number(path, line, 3);
    result.relative_tolerance = number(path, line, 4);
    result.max_iterations = integer(path, line, 5);
}

void readMethods(
    const std::filesystem::path& path,
    const ConfigLine& line,
    Methods& methods)
{
    oneValue(path, line);
    const std::string& value = line.tokens[1];
    if (line.tokens[0] == "interpolation") {
        if (value == "linear") methods.interpolation = InterpolationMethod::Linear;
        else if (value == "corrected") methods.interpolation = InterpolationMethod::Corrected;
        else invalid(path, line, "unknown interpolation method");
    } else if (line.tokens[0] == "gradient") {
        if (value == "greenGauss") methods.gradient = GradientMethod::GreenGauss;
        else if (value == "leastSquares") methods.gradient = GradientMethod::LeastSquares;
        else invalid(path, line, "unknown gradient method");
    } else if (line.tokens[0] == "convection") {
        if (value == "upwind") methods.convection = ConvectionMethod::Upwind;
        else if (value == "central") methods.convection = ConvectionMethod::Central;
        else invalid(path, line, "unknown convection method");
    } else if (line.tokens[0] == "diffusion") {
        if (value == "orthogonal") methods.diffusion = DiffusionMethod::Orthogonal;
        else if (value == "corrected") methods.diffusion = DiffusionMethod::Corrected;
        else if (value == "limitedCorrected") methods.diffusion = DiffusionMethod::LimitedCorrected;
        else invalid(path, line, "unknown diffusion method");
    } else if (line.tokens[0] == "time") {
        if (value == "steady") methods.time = TimeMethod::Steady;
        else if (value == "euler") methods.time = TimeMethod::Euler;
        else if (value == "bdf2") methods.time = TimeMethod::BDF2;
        else invalid(path, line, "unknown time method");
    } else {
        invalid(path, line, "unknown numerical method " + line.tokens[0]);
    }
}

}  // 匿名命名空间

ThermalCaseControl readThermalCase(const CaseDefinition& definition) {
    ThermalCaseControl result;
    bool density = false;
    bool heat_capacity = false;
    bool conductivity = false;
    bool source = false;
    for (const ConfigLine& line : readConfigLines(definition.physics_file)) {
        oneValue(definition.physics_file, line);
        const std::string& key = line.tokens[0];
        if (key == "density" && !density) {
            result.material.density = number(definition.physics_file, line, 1);
            density = true;
        } else if ((key == "heatCapacity" || key == "heat_capacity") && !heat_capacity) {
            result.material.heat_capacity = number(definition.physics_file, line, 1);
            heat_capacity = true;
        } else if ((key == "conductivity" || key == "thermalConductivity") && !conductivity) {
            result.material.conductivity = number(definition.physics_file, line, 1);
            conductivity = true;
        } else if ((key == "source" || key == "volumetricSource") && !source) {
            result.volumetric_source = number(definition.physics_file, line, 1);
            source = true;
        } else {
            invalid(definition.physics_file, line, "unknown or duplicate thermal property " + key);
        }
    }
    if (!density || !heat_capacity || !conductivity || !source) {
        throw std::runtime_error("thermal physics needs density, heatCapacity, conductivity, and source");
    }

    bool interpolation = false;
    bool gradient = false;
    bool convection = false;
    bool diffusion = false;
    bool time = false;
    bool start = false;
    bool end = false;
    bool delta = false;
    bool scalar_solver = false;
    for (const ConfigLine& line : readConfigLines(definition.numerics_file)) {
        const std::string& key = line.tokens[0];
        if (key == "interpolation" || key == "gradient" || key == "convection" ||
            key == "diffusion" || key == "time") {
            readMethods(definition.numerics_file, line, result.runtime.methods);
            if (key == "interpolation") interpolation = true;
            else if (key == "gradient") gradient = true;
            else if (key == "convection") convection = true;
            else if (key == "diffusion") diffusion = true;
            else time = true;
        } else if (key == "startTime") {
            oneValue(definition.numerics_file, line);
            if (start) invalid(definition.numerics_file, line, "duplicate startTime");
            result.runtime.time.start_time = number(definition.numerics_file, line, 1);
            start = true;
        } else if (key == "endTime") {
            oneValue(definition.numerics_file, line);
            if (end) invalid(definition.numerics_file, line, "duplicate endTime");
            result.runtime.time.end_time = number(definition.numerics_file, line, 1);
            end = true;
        } else if (key == "deltaT") {
            oneValue(definition.numerics_file, line);
            if (delta) invalid(definition.numerics_file, line, "duplicate deltaT");
            result.runtime.time.delta_t = number(definition.numerics_file, line, 1);
            delta = true;
        } else if (key == "scalarSolver") {
            if (scalar_solver) invalid(definition.numerics_file, line, "duplicate scalarSolver");
            readLinear(definition.numerics_file, line, result.runtime.scalar_solver);
            scalar_solver = true;
        } else {
            invalid(definition.numerics_file, line, "unknown heat numerical setting " + key);
        }
    }
    if (!interpolation || !gradient || !convection || !diffusion || !time || !start ||
        !end || !delta || !scalar_solver || result.runtime.methods.time == TimeMethod::Steady) {
        throw std::runtime_error("thermal numerics dictionary is incomplete or not transient");
    }
    result.material.validate();
    result.runtime.validate();
    return result;
}

}  // babelsim 命名空间
