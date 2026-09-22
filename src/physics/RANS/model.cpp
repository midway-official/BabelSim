#include "model.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>

namespace babelsim::rans {
namespace {

std::string selectedModel(const Parameters& settings) {
    if (!settings.contains("turbulenceModel")) return "none";
    const std::string selected = settings.word("turbulenceModel");
    std::string result;
    for (unsigned char character : selected) {
        if (std::isalnum(character)) result += static_cast<char>(std::tolower(character));
    }
    return result;
}

}  // 匿名命名空间

Model* create(
    Case& problem,
    VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity)
{
    const std::string model = selectedModel(problem.physics());
    if (model == "none" || model == "laminar" || model == "off") return nullptr;
    if (model == "sa" || model == "spalartallmaras") {
        return makeSpalartAllmaras(
            problem, velocity, face_flux, effective_viscosity,
            density, molecular_viscosity);
    }
    if (model == "komega" || model == "wilcox1988") {
        return makeKOmega(
            problem, velocity, face_flux, effective_viscosity,
            density, molecular_viscosity);
    }
    if (model == "kepsilon" || model == "standardkepsilon") {
        return makeKEpsilon(
            problem, velocity, face_flux, effective_viscosity,
            density, molecular_viscosity);
    }
    throw std::invalid_argument(
        "unsupported turbulenceModel; expected none, SA, kOmega or kEpsilon");
}

Turbulence::Turbulence(Case& problem, VectorField& velocity, const ScalarField& phi)
    : viscosity_(nullptr),
      momentumOptions_(readEquationControl(problem, "momentum", velocity,
          {"convection", "diffusion"}).spatial),
      model_(nullptr, destroy)
{
    const double rho = problem.physics().positive("density");
    const double mu = problem.physics().positive("dynamicViscosity");
    viscosity_ = &problem.createScalarField("muEffective", mu);
    model_.reset(create(problem, velocity, phi, *viscosity_, rho, mu));
}

void destroy(Model* model) noexcept { delete model; }

TensorField Turbulence::deviatoricStressRemainder(const VectorField& velocity) const {
    const auto gradU = math::grad(velocity, momentumOptions_);
    return effectiveViscosity() * (math::transpose(gradU)
        - (2.0 / 3.0) * math::isotropic(math::trace(gradU)));
}

}  // babelsim::rans 命名空间
