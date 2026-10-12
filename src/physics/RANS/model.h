#pragma once
#include "api.h"

namespace babelsim::rans {
Model* makeSpalartAllmaras(
    Case& problem,
    VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity,
    const OperatorOptions& momentum_options);
Model* makeKOmega(
    Case& problem,
    VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity,
    const OperatorOptions& momentum_options);
Model* makeKEpsilon(
    Case& problem,
    VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity,
    const OperatorOptions& momentum_options);

}  // babelsim::rans 命名空间
