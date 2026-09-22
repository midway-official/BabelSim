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

double number(const std::filesystem::path& path, const ConfigLine& line, std::size_t index) {
    try {
        std::size_t consumed = 0;
        const double value = std::stod(line.tokens.at(index), &consumed);
        if (consumed == line.tokens.at(index).size() && std::isfinite(value)) return value;
    } catch (const std::exception&) {
    }
    invalid(path, line, "expected a finite number");
}

void oneValue(const std::filesystem::path& path, const ConfigLine& line) {
    if (line.tokens.size() != 2) invalid(path, line, "expected one value for " + line.tokens.front());
}

InterpolationMethod interpolation(
    const std::filesystem::path& path, const ConfigLine& line, const std::string& value)
{
    if (value == "linear") return InterpolationMethod::Linear;
    if (value == "corrected") {
        return InterpolationMethod::Corrected;
    }
    invalid(path, line, "unknown interpolation method " + value);
}

GradientMethod gradient(
    const std::filesystem::path& path, const ConfigLine& line, const std::string& value)
{
    if (value == "greenGauss") return GradientMethod::GreenGauss;
    if (value == "leastSquares") return GradientMethod::LeastSquares;
    invalid(path, line, "unknown gradient method " + value);
}

ConvectionMethod convection(
    const std::filesystem::path& path, const ConfigLine& line, const std::string& value)
{
    if (value == "upwind") return ConvectionMethod::Upwind;
    if (value == "linearUpwind") return ConvectionMethod::LinearUpwind;
    if (value == "central") return ConvectionMethod::Central;
    invalid(path, line, "unknown convection method " + value);
}

DiffusionMethod diffusion(
    const std::filesystem::path& path, const ConfigLine& line, const std::string& value)
{
    if (value == "orthogonal") return DiffusionMethod::Orthogonal;
    if (value == "corrected") return DiffusionMethod::Corrected;
    if (value == "limitedCorrected") {
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

}  // 匿名命名空间

Methods readMethodsFile(const std::filesystem::path& path) {
    Methods result;
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

        if (key.rfind("equation.", 0) == 0 || key.rfind("operation.", 0) == 0) {
            oneValue(path, line);
            const auto split = key.rfind('.');
            const auto scope = key.substr(0, split);
            const auto option = key.substr(split + 1);
            const auto first = scope.find('.');
            const auto rest = scope.substr(first + 1);
            const auto term = rest.find(".term.");
            const bool valid = !rest.empty() &&
                (rest.find('.') == std::string::npos ||
                 (scope.rfind("equation.", 0) == 0 && term != std::string::npos &&
                  term > 0 && term + 6 < rest.size() &&
                  rest.substr(0, term).find('.') == std::string::npos &&
                  rest.substr(term + 6).find('.') == std::string::npos));
            if (!valid || split == first) invalid(path, line, "invalid numerical selector " + key);
            auto& entry = result.named[scope];
            entry.source = path.string() + ":" + std::to_string(line.number);
            auto set = [&](auto& target, auto value) {
                if (target) invalid(path, line, "duplicate numerical option " + key);
                target = value;
            };
            const auto& value = line.tokens[1];
            if (option == "interpolation") set(entry.options.interpolation, interpolation(path,line,value));
            else if (option == "gradient") set(entry.options.gradient, gradient(path,line,value));
            else if (option == "convection") set(entry.options.convection, convection(path,line,value));
            else if (option == "diffusion") set(entry.options.diffusion, diffusion(path,line,value));
            else if (option == "coefficientInterpolation") set(entry.options.coefficientInterpolation, interpolation(path,line,value));
            else if (option == "coefficientGradient") set(entry.options.coefficientGradient, gradient(path,line,value));
            else invalid(path, line, "unknown spatial numerical option " + option);
            continue;
        }
        invalid(path, line, "only time, equation.* and operation.* are valid methods entries");
    }
    if (!has_time) throw std::runtime_error("methods dictionary is incomplete: missing time");
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

}  // babelsim 命名空间
