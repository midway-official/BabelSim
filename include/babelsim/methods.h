#pragma once

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
};

}  // babelsim 命名空间
