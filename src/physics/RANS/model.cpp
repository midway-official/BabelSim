#include "model.h"

#include <stdexcept>
#include <string>

namespace babelsim::rans {
Model* createModel(
    Case& problem,
    VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    const builtin::FluidConfig& fluid,
    const OperatorOptions& momentumOptions)
{
    const auto& model = fluid.turbulenceModel;
    if (model == "none") return nullptr;
    if (model == "SA") {
        return makeSpalartAllmaras(
            problem, velocity, face_flux, effective_viscosity,
            fluid.density, fluid.dynamicViscosity, momentumOptions);
    }
    if (model == "kOmega") {
        return makeKOmega(
            problem, velocity, face_flux, effective_viscosity,
            fluid.density, fluid.dynamicViscosity, momentumOptions);
    }
    if (model == "kEpsilon") {
        return makeKEpsilon(
            problem, velocity, face_flux, effective_viscosity,
            fluid.density, fluid.dynamicViscosity, momentumOptions);
    }
    throw std::invalid_argument(
        "unsupported turbulenceModel; expected none, SA, kOmega or kEpsilon");
}

Turbulence::Turbulence(Case& problem, VectorField& velocity, const ScalarField& phi,
    const builtin::FluidConfig& fluid, const OperatorOptions& momentumOptions)
    : viscosity_(nullptr), momentumOptions_(momentumOptions), model_(nullptr, destroy)
{
    if (&velocity.mesh() != &problem.mesh() || &phi.mesh() != &problem.mesh() ||
        velocity.location() != FieldLocation::Cell || phi.location() != FieldLocation::Face)
        throw std::invalid_argument("RANS requires Case cell velocity and initialized face flux");
    momentumOptions_.requireComplete("RANS momentum options");
    if (!(fluid.density > 0) || !(fluid.dynamicViscosity > 0) ||
        !std::isfinite(fluid.density) || !std::isfinite(fluid.dynamicViscosity))
        throw std::invalid_argument("RANS requires positive finite fluid properties");
    viscosity_ = &problem.create.scalarField("muEffective", fluid.dynamicViscosity);
    model_.reset(createModel(problem, velocity, phi, *viscosity_, fluid, momentumOptions));
}

void destroy(Model* model) noexcept { delete model; }

TensorField Turbulence::deviatoricStressRemainder(const VectorField& velocity) const {
    const auto gradU = math::grad(velocity, momentumOptions_);
    return effectiveViscosity() * (math::transpose(gradU)
        - (2.0 / 3.0) * math::isotropic(math::trace(gradU)));
}

}  // babelsim::rans 命名空间
