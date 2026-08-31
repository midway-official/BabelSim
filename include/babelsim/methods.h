#pragma once

#include <string>
#include <vector>

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

struct Methods {
    InterpolationMethod interpolation = InterpolationMethod::Corrected;
    GradientMethod gradient = GradientMethod::LeastSquares;
    ConvectionMethod convection = ConvectionMethod::Upwind;
    DiffusionMethod diffusion = DiffusionMethod::Corrected;
    TimeMethod time = TimeMethod::Steady;

    // 默认格式适用于未单独指定的 Field。覆盖只在一次算子求值/装配开始时按 Field 名查找，
    // 不进入 cell/face 热循环，因此不会改变核心数值路径的数据布局。
    struct InterpolationOverride {
        std::string field;
        InterpolationMethod method = InterpolationMethod::Corrected;
    };
    struct GradientOverride {
        std::string field;
        GradientMethod method = GradientMethod::LeastSquares;
    };
    struct ConvectionOverride {
        std::string field;
        ConvectionMethod method = ConvectionMethod::Upwind;
    };
    struct DiffusionOverride {
        std::string field;
        DiffusionMethod method = DiffusionMethod::Corrected;
    };

    std::vector<InterpolationOverride> interpolation_overrides;
    std::vector<GradientOverride> gradient_overrides;
    std::vector<ConvectionOverride> convection_overrides;
    std::vector<DiffusionOverride> diffusion_overrides;

    InterpolationMethod interpolationFor(const std::string& field_name) const {
        for (const InterpolationOverride& entry : interpolation_overrides) {
            if (entry.field == field_name) return entry.method;
        }
        return interpolation;
    }
    GradientMethod gradientFor(const std::string& field_name) const {
        for (const GradientOverride& entry : gradient_overrides) {
            if (entry.field == field_name) return entry.method;
        }
        return gradient;
    }
    ConvectionMethod convectionFor(const std::string& field_name) const {
        for (const ConvectionOverride& entry : convection_overrides) {
            if (entry.field == field_name) return entry.method;
        }
        return convection;
    }
    DiffusionMethod diffusionFor(const std::string& field_name) const {
        for (const DiffusionOverride& entry : diffusion_overrides) {
            if (entry.field == field_name) return entry.method;
        }
        return diffusion;
    }
};

}  // babelsim 命名空间
