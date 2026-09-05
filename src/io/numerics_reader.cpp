#include "babelsim/numerics_io.h"

#include <algorithm>
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

InterpolationMethod interpolation(
    const std::filesystem::path& path, const ConfigLine& line, const std::string& value)
{
    if (value == "linear") return InterpolationMethod::Linear;
    if (value == "corrected" || value == "linearCorrected" || value == "linear_corrected") {
        return InterpolationMethod::Corrected;
    }
    invalid(path, line, "unknown interpolation method " + value);
}

GradientMethod gradient(
    const std::filesystem::path& path, const ConfigLine& line, const std::string& value)
{
    if (value == "greenGauss" || value == "green_gauss") return GradientMethod::GreenGauss;
    if (value == "leastSquares" || value == "least_squares") return GradientMethod::LeastSquares;
    invalid(path, line, "unknown gradient method " + value);
}

ConvectionMethod convection(
    const std::filesystem::path& path, const ConfigLine& line, const std::string& value)
{
    if (value == "upwind") return ConvectionMethod::Upwind;
    if (value == "linearUpwind" || value == "linear_upwind" ||
        value == "secondOrderUpwind") return ConvectionMethod::LinearUpwind;
    if (value == "central") return ConvectionMethod::Central;
    invalid(path, line, "unknown convection method " + value);
}

DiffusionMethod diffusion(
    const std::filesystem::path& path, const ConfigLine& line, const std::string& value)
{
    if (value == "orthogonal") return DiffusionMethod::Orthogonal;
    if (value == "corrected") return DiffusionMethod::Corrected;
    if (value == "limitedCorrected" || value == "limited_corrected") {
        return DiffusionMethod::LimitedCorrected;
    }
    invalid(path, line, "unknown diffusion method " + value);
}

TimeMethod timeMethod(
    const std::filesystem::path& path, const ConfigLine& line, const std::string& value)
{
    if (value == "steady") return TimeMethod::Steady;
    if (value == "euler") return TimeMethod::Euler;
    if (value == "bdf2") return TimeMethod::BDF2;
    invalid(path, line, "unknown time method " + value);
}

// 四种空间方法共用相同的“Field 名不能重复”规则；类型仍由各自的 enum 保证。
template <typename Entry, typename Method>
void addOverride(
    const std::filesystem::path& path,
    const ConfigLine& line,
    std::vector<Entry>& entries,
    Method method)
{
    const std::string& field = line.tokens[1];
    if (field.empty() || std::find_if(entries.begin(), entries.end(), [&](const Entry& entry) {
            return entry.field == field;
        }) != entries.end()) {
        invalid(path, line, "duplicate or empty Field method override");
    }
    entries.push_back({field, method});
}

}  // 匿名命名空间

Methods readMethodsFile(const std::filesystem::path& path) {
    Methods result;
    bool has_interpolation = false;
    bool has_gradient = false;
    bool has_convection = false;
    bool has_diffusion = false;
    bool has_time = false;
    for (const ConfigLine& line : readConfigLines(path)) {
        const std::string& key = line.tokens.front();
        if (key == "time") {
            oneValue(path, line);
            if (has_time) invalid(path, line, "duplicate time method");
            result.time = timeMethod(path, line, line.tokens[1]);
            has_time = true;
            continue;
        }
        if (line.tokens.size() != 2 && line.tokens.size() != 3) {
            invalid(path, line, "method needs a default value or a Field name and value");
        }
        const bool override = line.tokens.size() == 3;
        const std::string& value = line.tokens[override ? 2 : 1];
        if (key == "interpolation") {
            const InterpolationMethod method = interpolation(path, line, value);
            if (override) addOverride(path, line, result.interpolation_overrides, method);
            else if (!has_interpolation) { result.interpolation = method; has_interpolation = true; }
            else invalid(path, line, "duplicate interpolation method");
        } else if (key == "gradient") {
            const GradientMethod method = gradient(path, line, value);
            if (override) addOverride(path, line, result.gradient_overrides, method);
            else if (!has_gradient) { result.gradient = method; has_gradient = true; }
            else invalid(path, line, "duplicate gradient method");
        } else if (key == "convection") {
            const ConvectionMethod method = convection(path, line, value);
            if (override) addOverride(path, line, result.convection_overrides, method);
            else if (!has_convection) { result.convection = method; has_convection = true; }
            else invalid(path, line, "duplicate convection method");
        } else if (key == "diffusion") {
            const DiffusionMethod method = diffusion(path, line, value);
            if (override) addOverride(path, line, result.diffusion_overrides, method);
            else if (!has_diffusion) { result.diffusion = method; has_diffusion = true; }
            else invalid(path, line, "duplicate diffusion method");
        } else {
            invalid(path, line, "unknown numerical method " + key);
        }
    }
    if (!has_interpolation || !has_gradient || !has_convection || !has_diffusion || !has_time) {
        throw std::runtime_error("methods dictionary is incomplete");
    }
    return result;
}

TimeControl readTimeControlFile(const std::filesystem::path& path) {
    TimeControl result;
    bool start = false;
    bool end = false;
    bool delta = false;
    for (const ConfigLine& line : readConfigLines(path)) {
        oneValue(path, line);
        const std::string& key = line.tokens.front();
        if (key == "startTime") {
            if (start) invalid(path, line, "duplicate startTime");
            result.start_time = number(path, line, 1);
            start = true;
        } else if (key == "endTime") {
            if (end) invalid(path, line, "duplicate endTime");
            result.end_time = number(path, line, 1);
            end = true;
        } else if (key == "deltaT") {
            if (delta) invalid(path, line, "duplicate deltaT");
            result.delta_t = number(path, line, 1);
            delta = true;
        } else {
            invalid(path, line, "unknown run control " + key);
        }
    }
    if (!start || !end || !delta) throw std::runtime_error("control dictionary is incomplete");
    result.validate();
    return result;
}

void readLinearSolverLine(
    const std::filesystem::path& path,
    const ConfigLine& line,
    LinearSolverConfig& result)
{
    if (line.tokens.size() < 6) {
        invalid(path, line, "linear solver needs method, preconditioner, tolerances, and iterations");
    }
    if (line.tokens[1] == "cg") result.solver = LinearSolverType::ConjugateGradient;
    else if (line.tokens[1] == "bicgstab") result.solver = LinearSolverType::BiCGSTAB;
    else if (line.tokens[1] == "gmres") result.solver = LinearSolverType::GMRES;
    else if (line.tokens[1] == "amg") result.solver = LinearSolverType::AlgebraicMultigrid;
    else invalid(path, line, "unknown linear solver " + line.tokens[1]);
    if (line.tokens[2] == "none") {
        result.preconditioner = PreconditionerType::None;
    } else if (line.tokens[2] == "incompleteCholesky" || line.tokens[2] == "incomplete_cholesky") {
        result.preconditioner = PreconditionerType::IncompleteCholesky;
    } else if (line.tokens[2] == "ilut") {
        result.preconditioner = PreconditionerType::ILUT;
    } else if (line.tokens[2] == "amg") {
        result.preconditioner = PreconditionerType::AlgebraicMultigrid;
    } else {
        invalid(path, line, "unknown preconditioner " + line.tokens[2]);
    }
    result.absolute_tolerance = number(path, line, 3);
    result.relative_tolerance = number(path, line, 4);
    result.max_iterations = integer(path, line, 5);
    for (std::size_t index = 6; index < line.tokens.size(); ++index) {
        const std::string& option = line.tokens[index];
        const std::size_t separator = option.find('=');
        if (separator == std::string::npos || separator == 0 ||
            separator + 1 == option.size()) {
            invalid(path, line, "linear solver option must use name=value");
        }
        ConfigLine value_line = line;
        value_line.tokens = {option.substr(0, separator), option.substr(separator + 1)};
        if (value_line.tokens[0] == "gmresRestart") {
            result.gmres_restart = integer(path, value_line, 1);
        } else if (value_line.tokens[0] == "amgMaxLevels") {
            result.amg_max_levels = integer(path, value_line, 1);
        } else if (value_line.tokens[0] == "amgCoarseSize") {
            result.amg_coarse_size = integer(path, value_line, 1);
        } else if (value_line.tokens[0] == "amgSmoothingSteps") {
            result.amg_smoothing_steps = integer(path, value_line, 1);
        } else {
            invalid(path, line, "unknown linear solver option " + value_line.tokens[0]);
        }
    }
}

}  // babelsim 命名空间
