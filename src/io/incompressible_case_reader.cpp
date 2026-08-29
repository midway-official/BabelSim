#include "babelsim/incompressible_io.h"

#include "babelsim/config.h"

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
        "invalid " + path.string() + ":" + std::to_string(line.number) +
        ": " + message);
}

double number(const std::filesystem::path& path, const ConfigLine& line, std::size_t at) {
    try {
        std::size_t consumed = 0;
        const double value = std::stod(line.tokens.at(at), &consumed);
        if (consumed != line.tokens.at(at).size()) invalid(path, line, "expected a number");
        return value;
    } catch (const std::exception&) {
        invalid(path, line, "expected a number");
    }
}

int integer(const std::filesystem::path& path, const ConfigLine& line, std::size_t at) {
    try {
        std::size_t consumed = 0;
        const int value = std::stoi(line.tokens.at(at), &consumed);
        if (consumed != line.tokens.at(at).size()) invalid(path, line, "expected an integer");
        return value;
    } catch (const std::exception&) {
        invalid(path, line, "expected an integer");
    }
}

void oneValue(const std::filesystem::path& path, const ConfigLine& line) {
    if (line.tokens.size() != 2) invalid(path, line, "expected one value for " + line.tokens.front());
}

LinearSolverType solverType(
    const std::filesystem::path& path,
    const ConfigLine& line,
    const std::string& value)
{
    if (value == "cg") return LinearSolverType::ConjugateGradient;
    if (value == "bicgstab") return LinearSolverType::BiCGSTAB;
    invalid(path, line, "unknown linear solver " + value);
}

PreconditionerType preconditionerType(
    const std::filesystem::path& path,
    const ConfigLine& line,
    const std::string& value)
{
    if (value == "incompleteCholesky" || value == "incomplete_cholesky") {
        return PreconditionerType::IncompleteCholesky;
    }
    if (value == "ilut") return PreconditionerType::ILUT;
    invalid(path, line, "unknown preconditioner " + value);
}

void linearSolver(
    const std::filesystem::path& path,
    const ConfigLine& line,
    LinearSolverConfig& config)
{
    if (line.tokens.size() != 6) {
        invalid(path, line, "linear solver needs method, preconditioner, tolerances, and iterations");
    }
    config.solver = solverType(path, line, line.tokens[1]);
    config.preconditioner = preconditionerType(path, line, line.tokens[2]);
    config.absolute_tolerance = number(path, line, 3);
    config.relative_tolerance = number(path, line, 4);
    config.max_iterations = integer(path, line, 5);
}

}  // 匿名命名空间

IncompressibleCaseControl readIncompressibleCase(const CaseDefinition& definition) {
    IncompressibleCaseControl result;
    bool density = false;
    bool viscosity = false;
    for (const ConfigLine& line : readConfigLines(definition.physics_file)) {
        const std::string& key = line.tokens.front();
        oneValue(definition.physics_file, line);
        if (key == "density" && !density) {
            result.fluid.density = number(definition.physics_file, line, 1);
            density = true;
        } else if ((key == "dynamicViscosity" || key == "dynamic_viscosity") && !viscosity) {
            result.fluid.dynamic_viscosity = number(definition.physics_file, line, 1);
            viscosity = true;
        } else {
            invalid(definition.physics_file, line, "unknown or duplicate entry " + key);
        }
    }
    if (!density || !viscosity) {
        throw std::runtime_error("incompressible physics needs density and dynamic_viscosity");
    }

    bool interpolation = false;
    bool gradient = false;
    bool convection = false;
    bool diffusion = false;
    bool time = false;
    bool iterations = false;
    bool non_orthogonal_corrections = false;
    bool velocity_relaxation = false;
    bool pressure_relaxation = false;
    bool continuity_tolerance = false;
    bool velocity_tolerance = false;
    bool velocity_solver = false;
    bool pressure_solver = false;
    for (const ConfigLine& line : readConfigLines(definition.numerics_file)) {
        const std::string& key = line.tokens.front();
        if (key == "interpolation") {
            oneValue(definition.numerics_file, line);
            if (interpolation) invalid(definition.numerics_file, line, "duplicate interpolation");
            if (line.tokens[1] == "linear") {
                result.methods.interpolation = InterpolationMethod::Linear;
            } else if (line.tokens[1] == "corrected" ||
                       line.tokens[1] == "linearCorrected" ||
                       line.tokens[1] == "linear_corrected") {
                result.methods.interpolation = InterpolationMethod::Corrected;
            } else {
                invalid(definition.numerics_file, line, "unknown interpolation method");
            }
            interpolation = true;
        } else if (key == "gradient") {
            oneValue(definition.numerics_file, line);
            if (gradient) invalid(definition.numerics_file, line, "duplicate gradient");
            if (line.tokens[1] == "greenGauss" || line.tokens[1] == "green_gauss") {
                result.methods.gradient = GradientMethod::GreenGauss;
            } else if (line.tokens[1] == "leastSquares" || line.tokens[1] == "least_squares") {
                result.methods.gradient = GradientMethod::LeastSquares;
            } else {
                invalid(definition.numerics_file, line, "unknown gradient method");
            }
            gradient = true;
        } else if (key == "convection") {
            oneValue(definition.numerics_file, line);
            if (convection) invalid(definition.numerics_file, line, "duplicate convection");
            if (line.tokens[1] == "upwind") result.methods.convection = ConvectionMethod::Upwind;
            else if (line.tokens[1] == "central") result.methods.convection = ConvectionMethod::Central;
            else invalid(definition.numerics_file, line, "unknown convection method");
            convection = true;
        } else if (key == "diffusion") {
            oneValue(definition.numerics_file, line);
            if (diffusion) invalid(definition.numerics_file, line, "duplicate diffusion");
            if (line.tokens[1] == "orthogonal") result.methods.diffusion = DiffusionMethod::Orthogonal;
            else if (line.tokens[1] == "corrected") result.methods.diffusion = DiffusionMethod::Corrected;
            else if (line.tokens[1] == "limitedCorrected" ||
                     line.tokens[1] == "limited_corrected") {
                result.methods.diffusion = DiffusionMethod::LimitedCorrected;
            }
            else invalid(definition.numerics_file, line, "unknown diffusion method");
            diffusion = true;
        } else if (key == "time") {
            oneValue(definition.numerics_file, line);
            if (time) invalid(definition.numerics_file, line, "duplicate time method");
            if (line.tokens[1] == "steady") result.methods.time = TimeMethod::Steady;
            else if (line.tokens[1] == "euler") result.methods.time = TimeMethod::Euler;
            else if (line.tokens[1] == "bdf2") result.methods.time = TimeMethod::BDF2;
            else invalid(definition.numerics_file, line, "unknown time method");
            time = true;
        } else if (key == "maxIterations" || key == "max_iterations") {
            oneValue(definition.numerics_file, line);
            if (iterations) invalid(definition.numerics_file, line, "duplicate max_iterations");
            result.simple.max_iterations = integer(definition.numerics_file, line, 1);
            iterations = true;
        } else if (key == "nonOrthogonalCorrections" ||
                   key == "non_orthogonal_corrections") {
            oneValue(definition.numerics_file, line);
            if (non_orthogonal_corrections) {
                invalid(
                    definition.numerics_file, line,
                    "duplicate non_orthogonal_corrections");
            }
            result.simple.non_orthogonal_corrections =
                integer(definition.numerics_file, line, 1);
            non_orthogonal_corrections = true;
        } else if (key == "velocityRelaxation" || key == "velocity_relaxation") {
            oneValue(definition.numerics_file, line);
            if (velocity_relaxation) invalid(definition.numerics_file, line, "duplicate velocity_relaxation");
            result.simple.velocity_relaxation = number(definition.numerics_file, line, 1);
            velocity_relaxation = true;
        } else if (key == "pressureRelaxation" || key == "pressure_relaxation") {
            oneValue(definition.numerics_file, line);
            if (pressure_relaxation) invalid(definition.numerics_file, line, "duplicate pressure_relaxation");
            result.simple.pressure_relaxation = number(definition.numerics_file, line, 1);
            pressure_relaxation = true;
        } else if (key == "continuityTolerance" || key == "continuity_tolerance") {
            oneValue(definition.numerics_file, line);
            if (continuity_tolerance) invalid(definition.numerics_file, line, "duplicate continuity_tolerance");
            result.simple.continuity_tolerance = number(definition.numerics_file, line, 1);
            continuity_tolerance = true;
        } else if (key == "velocityTolerance" || key == "velocity_tolerance") {
            oneValue(definition.numerics_file, line);
            if (velocity_tolerance) invalid(definition.numerics_file, line, "duplicate velocity_tolerance");
            result.simple.velocity_tolerance = number(definition.numerics_file, line, 1);
            velocity_tolerance = true;
        } else if (key == "velocitySolver" || key == "velocity_solver") {
            if (velocity_solver) invalid(definition.numerics_file, line, "duplicate velocity_solver");
            linearSolver(definition.numerics_file, line, result.simple.velocity_solver);
            velocity_solver = true;
        } else if (key == "pressureSolver" || key == "pressure_solver") {
            if (pressure_solver) invalid(definition.numerics_file, line, "duplicate pressure_solver");
            linearSolver(definition.numerics_file, line, result.simple.pressure_solver);
            pressure_solver = true;
        } else {
            invalid(definition.numerics_file, line, "unknown entry " + key);
        }
    }
    if (!interpolation || !gradient || !convection || !diffusion || !time || !iterations ||
        !velocity_relaxation || !pressure_relaxation || !continuity_tolerance ||
        !velocity_tolerance || !velocity_solver || !pressure_solver) {
        throw std::runtime_error("incompressible numerics dictionary is incomplete");
    }
    result.fluid.validate();
    result.simple.validate();
    return result;
}

}  // babelsim 命名空间
