#pragma once
#include "babelsim/model_description.h"

namespace babelsim::builtin {
// Plain values, with no duplicated defaults. Descriptions bind every member.
struct HeatConfig { double density, heatCapacity, conductivity, source; };
struct TransportConfig { double storage, diffusivity, source; };
struct FluidConfig { double density, dynamicViscosity; std::string turbulenceModel; };
struct ConvergenceConfig { double continuity, velocity, momentum, pressureCorrection; };
struct SimpleCoupling { int maxIterations, nonOrthogonalCorrections; double velocityRelaxation, pressureRelaxation; };
struct SimpleConfig { FluidConfig fluid; SimpleCoupling coupling; ConvergenceConfig convergence; };
struct PisoCoupling {
    int maxIterations, correctors, maxCorrectors, nonOrthogonalCorrections, initialNonOrthogonalCorrections;
    double tolerance, velocityRelaxation;
    int initializePotentialFlow;
};
struct TurbulenceCoupling { std::string mode; int maxIterations; std::optional<int> correctors; };
struct PisoConfig { FluidConfig fluid; PisoCoupling coupling; ConvergenceConfig convergence; TurbulenceCoupling turbulence; };
struct RansControls { double relaxation, tolerance; };
struct SAConfig {
    RansControls controls;
    double nuTildaMin, wallDistanceMin, cb1, cb2, sigma, kappa, cw2, cw3, cv1, ct3, ct4, cw1;
};
struct KOmegaConfig {
    RansControls controls;
    double kMin, omegaMin, viscosityRatioMax, betaStar, beta, gamma, sigmaK, sigmaOmega;
};
struct KEpsilonConfig {
    RansControls controls;
    double kMin, epsilonMin, Cmu, C1, C2, sigmaK, sigmaEpsilon;
};
const ConfigurationDescription<HeatConfig>& heat();
const ConfigurationDescription<TransportConfig>& transport();
const ConfigurationDescription<SimpleConfig>& simple();
const ConfigurationDescription<SimpleConfig>& transientSimple();
const ConfigurationDescription<PisoConfig>& piso();
const ConfigurationDescription<SAConfig>& spalartAllmaras();
const ConfigurationDescription<KOmegaConfig>& kOmega();
const ConfigurationDescription<KEpsilonConfig>& kEpsilon();
} // namespace babelsim::builtin
