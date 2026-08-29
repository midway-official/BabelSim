#pragma once

namespace babelsim {

enum class InterpolationMethod {
    Linear,
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
};

enum class TimeMethod {
    Steady,
    Euler,
    BDF2,
};

struct Methods {
    InterpolationMethod interpolation = InterpolationMethod::Linear;
    GradientMethod gradient = GradientMethod::GreenGauss;
    ConvectionMethod convection = ConvectionMethod::Upwind;
    DiffusionMethod diffusion = DiffusionMethod::Corrected;
    TimeMethod time = TimeMethod::Steady;
};

}  // namespace babelsim

