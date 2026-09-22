#pragma once

#include <string>
#include <optional>
#include <map>
#include <stdexcept>

namespace babelsim {

enum class InterpolationMethod {
    Linear,
    Corrected,
};

enum class GradientMethod {
    GreenGauss,
    LeastSquares,
};

enum class ConvectionMethod {
    Upwind,
    LinearUpwind,
    Central,
};

enum class DiffusionMethod {
    Orthogonal,
    Corrected,
    LimitedCorrected,
};

enum class TimeMethod {
    Steady,
    Euler,
    BDF2,
};

// Per-call spatial choices. A numerical operation must receive a complete
// equation/operation profile; there is no global or field-level fallback.
// Coefficient reconstruction is independent from reconstruction of the unknown.
struct OperatorOptions {
    std::optional<InterpolationMethod> interpolation;
    std::optional<GradientMethod> gradient;
    std::optional<ConvectionMethod> convection;
    std::optional<DiffusionMethod> diffusion;
    std::optional<InterpolationMethod> coefficientInterpolation;
    std::optional<GradientMethod> coefficientGradient;

    void requireComplete(const std::string& context) const {
        if (!interpolation || !gradient || !convection || !diffusion)
            throw std::invalid_argument("incomplete spatial configuration for " + context);
    }

    // Options for an explicitly evaluated diffusion coefficient interpolation.
    OperatorOptions coefficientOptions() const {
        auto result = *this;
        if (coefficientInterpolation) result.interpolation = coefficientInterpolation;
        if (coefficientGradient) result.gradient = coefficientGradient;
        return result;
    }

    void overlay(const OperatorOptions& other) {
        if (other.interpolation) interpolation = other.interpolation;
        if (other.gradient) gradient = other.gradient;
        if (other.convection) convection = other.convection;
        if (other.diffusion) diffusion = other.diffusion;
        if (other.coefficientInterpolation) coefficientInterpolation = other.coefficientInterpolation;
        if (other.coefficientGradient) coefficientGradient = other.coefficientGradient;
    }
};

struct Methods {
    TimeMethod time = TimeMethod::Steady;

    struct NamedOptions {
        OperatorOptions options;
        std::string source;
        mutable bool used = false;
    };
    // Keys are equation.<name>[.term.<name>] or operation.<name>.
    std::map<std::string, NamedOptions> named;

    OperatorOptions equationOptions(const std::string& name) const {
        const auto key = "equation." + name;
        const auto found = named.find(key);
        if (found == named.end())
            throw std::invalid_argument("missing numerical configuration " + key);
        found->second.used = true;
        found->second.options.requireComplete(key);
        return found->second.options;
    }
    OperatorOptions operationOptions(const std::string& name) const {
        const auto key = "operation." + name;
        const auto found = named.find(key);
        if (found == named.end())
            throw std::invalid_argument("missing numerical configuration " + key);
        found->second.used = true;
        found->second.options.requireComplete(key);
        return found->second.options;
    }
    OperatorOptions equationTermOptions(const std::string& equation, const std::string& term) const {
        const auto key = "equation." + equation + ".term." + term;
        const auto found = named.find(key);
        if (found == named.end()) return {};
        found->second.used = true;
        return found->second.options;
    }
    void requireAllUsed() const {
        for (const auto& entry : named)
            if (!entry.second.used)
                throw std::runtime_error(entry.second.source + ": unused numerical configuration " + entry.first);
    }
};

}  // babelsim 命名空间
