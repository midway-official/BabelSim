#include "model.h"

#include <stdexcept>
#include <string>

namespace babelsim::rans {
namespace {

std::string selectedModel(const Parameters& settings) {
    if (!settings.contains("turbulenceModel")) return "none";
    return settings.word("turbulenceModel");
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
    if (model == "none") return nullptr;
    if (model == "SA") {
        return makeSpalartAllmaras(
            problem, velocity, face_flux, effective_viscosity,
            density, molecular_viscosity);
    }
    if (model == "kOmega") {
        return makeKOmega(
            problem, velocity, face_flux, effective_viscosity,
            density, molecular_viscosity);
    }
    if (model == "kEpsilon") {
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
