#include "babelsim/builtin_descriptions.h"

#include <limits>

namespace babelsim::builtin {
namespace {
const std::vector<std::string>& supportedBoundaries() {
    // Descriptions may be requested by a registration in another translation
    // unit before this unit's global objects are initialized.
    static const std::vector<std::string> result{
        "fixedValue", "fixedGradient", "zeroGradient", "symmetry", "inletOutlet"};
    return result;
}

ParameterDescription number(const std::string& name, const std::string& unit,
    const std::string& explanation, std::optional<double> fallback,
    std::optional<double> minimum = std::nullopt, bool exclusive = false,
    ParameterDictionary dictionary = ParameterDictionary::Physics) {
    ParameterDescription result;
    result.name = name;
    result.dictionary = dictionary;
    result.unit = unit;
    result.explanation = explanation;
    result.required = !fallback;
    if (fallback) result.default_value = *fallback;
    result.minimum = minimum;
    result.exclusive_minimum = exclusive;
    return result;
}
ParameterDescription integer(const std::string& name, int fallback, int minimum,
    int maximum, const std::string& explanation) {
    auto result = number(name, "1", explanation, std::nullopt);
    result.dictionary = ParameterDictionary::Solution;
    result.type = ParameterType::Integer;
    result.required = false;
    result.default_value = fallback;
    result.minimum = minimum;
    result.maximum = maximum;
    return result;
}
ParameterDescription word(const std::string& name, const std::string& fallback,
    const std::string& explanation, ParameterDictionary dictionary = ParameterDictionary::Physics) {
    auto result = number(name, "1", explanation, std::nullopt);
    result.dictionary = dictionary;
    result.type = ParameterType::Word;
    result.required = false;
    result.default_value = fallback;
    return result;
}
ParameterDescription fraction(const std::string& name, double fallback, const std::string& explanation) {
    auto result = number(name, "1", explanation, fallback, 0.0, true, ParameterDictionary::Solution);
    result.maximum = 1.0;
    return result;
}
FieldDescription scalar(const std::string& name, const std::string& unit,
    const std::string& explanation, FieldRole role = FieldRole::Input) {
    return {name, DescribedFieldType::Scalar, DescribedFieldLocation::Cell,
        role, unit, explanation, "0"};
}
FieldDescription vector(const std::string& name, const std::string& unit,
    const std::string& explanation, FieldRole role = FieldRole::Input) {
    return {name, DescribedFieldType::Vector, DescribedFieldLocation::Cell,
        role, unit, explanation, "0 0 0"};
}
FieldDescription flux() {
    return {"phi", DescribedFieldType::Scalar, DescribedFieldLocation::Face,
        FieldRole::Derived, "m^3/s", "Oriented volume flux supplied by the flow solver.", "0"};
}
ModuleDescription flow(const std::string& name, const std::string& explanation,
    std::vector<TimeMethod> methods) {
    ModuleDescription result;
    result.name = name;
    result.explanation = explanation;
    result.fields = {vector("U", "m/s", "Velocity."), scalar("p", "Pa", "Pressure."), flux(),
        scalar("muEffective", "Pa*s", "Molecular plus turbulent dynamic viscosity.", FieldRole::Derived)};
    auto density = number("density", "kg/m^3", "Constant fluid density.", std::nullopt, 0.0, true);
    density.template_value = 1.0;
    auto viscosity = number("dynamicViscosity", "Pa*s", "Molecular dynamic viscosity.", std::nullopt, 0.0, true);
    viscosity.template_value = 0.01;
    auto model = word("turbulenceModel", "none", "Optional RANS closure.");
    model.choices = {"none", "SA", "kOmega", "kEpsilon"};
    result.parameters = {density, viscosity, model};
    result.equations = {{"momentum", "U", {"convection", "diffusion"}, "Velocity equation."},
        {"pressureCorrection", "p", {"diffusion"}, "Pressure correction/projection equation."}};
    result.time_methods = std::move(methods);
    result.boundary_types = supportedBoundaries();
    result.limitations = {"Constant-density incompressible flow. Pressure and flux correction remain explicit solver algorithms.",
        "A selected RANS model adds its own description, input fields, and numerical profiles."};
    return result;
}
void addConvergence(ModuleDescription& result, double velocityRelaxation) {
    result.parameters.push_back(integer("maxIterations", 1000, 1, std::numeric_limits<int>::max(),
        "Maximum SIMPLE coupling iterations (per physical step for transientSimple)."));
    result.parameters.push_back(integer("nonOrthogonalCorrections", 1, 0, 20, "Extra pressure nonorthogonal corrections."));
    result.parameters.push_back(fraction("velocityRelaxation", velocityRelaxation, "Momentum underrelaxation."));
    result.parameters.push_back(fraction("pressureRelaxation", 0.3, "Pressure underrelaxation."));
    for (const auto& parameter : {
        number("continuityTolerance", "1", "Normalized continuity acceptance tolerance.", 1e-8, 0.0, true, ParameterDictionary::Solution),
        number("velocityTolerance", "1", "Relative velocity change tolerance.", 1e-7, 0.0, true, ParameterDictionary::Solution),
        number("momentumTolerance", "1", "Normalized momentum residual tolerance.", 1e-6, 0.0, true, ParameterDictionary::Solution),
        number("pressureCorrectionTolerance", "1", "Normalized pressure correction residual tolerance.", 1e-6, 0.0, true, ParameterDictionary::Solution)})
        result.parameters.push_back(parameter);
}
ModuleDescription turbulence(const std::string& name, const std::string& explanation) {
    ModuleDescription result;
    result.name = name;
    result.explanation = explanation;
    result.fields = {vector("U", "m/s", "Velocity supplied by the flow solver.", FieldRole::Derived), flux(),
        scalar("muEffective", "Pa*s", "Effective dynamic viscosity.", FieldRole::Derived),
        scalar("mut", "Pa*s", "Turbulent dynamic viscosity.", FieldRole::Derived)};
    result.parameters = {fraction("turbulenceRelaxation", 0.7, "Turbulence transport underrelaxation."),
        number("turbulenceTolerance", "1", "Model transport residual and change tolerance.", 1e-6, 0.0, true, ParameterDictionary::Solution)};
    result.time_methods = {TimeMethod::Steady, TimeMethod::Euler, TimeMethod::BDF2};
    result.boundary_types = supportedBoundaries();
    result.compatible_solvers = {"simple", "transientSimple", "piso"};
    result.limitations = {"Constant-density incompressible eddy-viscosity closure; density and molecular viscosity come from the host solver.",
        "The host solver owns flow coupling and saves model histories once per physical step."};
    return result;
}
ParameterValue saCw1(const ParameterValues& parameters) {
    const double kappa = parameters.number("saKappa");
    return parameters.number("saCb1") / (kappa * kappa) +
        (1.0 + parameters.number("saCb2")) / parameters.number("saSigma");
}
ParameterValue maxCorrectors(const ParameterValues& parameters) { return parameters.integer("nCorrectors"); }
ParameterValue initialCorrections(const ParameterValues& parameters) { return parameters.integer("nonOrthogonalCorrections"); }
} // namespace

const ModuleDescription& heat() {
    static const ModuleDescription description = [] {
        ModuleDescription result;
        result.name = "heat";
        result.explanation = "Transient heat conduction: rho*cp*dT/dt - div(k*grad(T)) = Q.";
        result.fields = {scalar("T", "K", "Temperature.")};
        auto density = number("density", "kg/m^3", "Constant material density.", std::nullopt, 0.0, true);
        density.template_value = 1.0;
        auto capacity = number("heatCapacity", "J/(kg*K)", "Specific heat capacity.", std::nullopt, 0.0, true);
        capacity.template_value = 1.0;
        auto conductivity = number("conductivity", "W/(m*K)", "Thermal conductivity.", std::nullopt, 0.0);
        conductivity.template_value = 0.1;
        auto source = number("source", "W/m^3", "Volumetric heat source; integration is performed by equ::source.", std::nullopt);
        source.template_value = 1.0;
        result.parameters = {density, capacity, conductivity, source};
        result.equations = {{"temperature", "T", {"diffusion"}, "Nonconservative capacity derivative rho*cp*dT/dt with constant capacity."}};
        result.time_methods = {TimeMethod::Euler, TimeMethod::BDF2};
        result.boundary_types = supportedBoundaries();
        result.limitations = {"Capacity is constant. This module does not implement d(rho*cp*T)/dt for time-varying coefficients."};
        result.supports_case_template = true;
        return result;
    }();
    return description;
}
const ModuleDescription& transport() {
    static const ModuleDescription description = [] {
        ModuleDescription result;
        result.name = "transport";
        result.explanation = "Scalar transport: storage*dC/dt + div(phi*C) - div(D*grad(C)) = S.";
        result.fields = {scalar("C", "1", "Transported concentration."), vector("U", "m/s", "Prescribed velocity."), flux()};
        auto storage = number("storage", "1", "Constant scalar capacity.", std::nullopt, 0.0, true);
        storage.template_value = 1.0;
        auto diffusivity = number("diffusivity", "m^2/s", "Scalar diffusivity.", std::nullopt, 0.0);
        diffusivity.template_value = 0.1;
        auto source = number("source", "1/s", "Volumetric scalar source; integration is performed by equ::source.", std::nullopt);
        source.template_value = 1.0;
        result.parameters = {storage, diffusivity, source};
        result.equations = {{"transport", "C", {"convection", "diffusion"}, "Prescribed velocity; no flow coupling."}};
        result.time_methods = {TimeMethod::Euler, TimeMethod::BDF2};
        result.boundary_types = supportedBoundaries();
        result.limitations = {"Input U remains fixed. phi is a volume flux; storage is a constant capacity.",
            "The generic tutorial uses U=0 to give the analytic source accumulation C=S*t/storage. Configure velocity and inflow/outflow boundaries together when adding convection."};
        result.supports_case_template = true;
        return result;
    }();
    return description;
}
const ModuleDescription& simple() {
    static const ModuleDescription description = [] {
        auto result = flow("simple", "Steady incompressible SIMPLE coupling.", {TimeMethod::Steady});
        result.fields.push_back(scalar("pPrime", "Pa", "Homogeneous pressure correction.", FieldRole::Derived));
        result.equations[1].unknown = "pPrime";
        addConvergence(result, 0.7);
        return result;
    }();
    return description;
}
const ModuleDescription& transientSimple() {
    static const ModuleDescription description = [] {
        auto result = flow("transientSimple", "Transient incompressible SIMPLE coupling.", {TimeMethod::Euler, TimeMethod::BDF2});
        result.fields.push_back(scalar("pPrime", "Pa", "Homogeneous pressure correction.", FieldRole::Derived));
        result.equations[1].unknown = "pPrime";
        addConvergence(result, 0.7);
        return result;
    }();
    return description;
}
const ModuleDescription& piso() {
    static const ModuleDescription description = [] {
        auto result = flow("piso", "Transient incompressible PISO pressure-velocity coupling.", {TimeMethod::Euler, TimeMethod::BDF2});
        result.parameters.push_back(integer("maxIterations", 1, 1, std::numeric_limits<int>::max(), "Momentum passes per physical step."));
        result.parameters.push_back(integer("nCorrectors", 2, 1, 100, "Minimum pressure-velocity corrections per pass."));
        auto maximum = integer("maxCorrectors", 2, 1, 100, "Maximum pressure-velocity corrections per pass.");
        maximum.default_value.reset(); maximum.default_factory = maxCorrectors;
        maximum.default_explanation = "nCorrectors";
        result.parameters.push_back(maximum);
        result.parameters.push_back(number("couplingTolerance", "1", "Optional corrected frozen momentum residual gate.",
            std::numeric_limits<double>::max(), 0.0, true, ParameterDictionary::Solution));
        result.parameters.push_back(integer("nonOrthogonalCorrections", 1, 0, 20, "Extra nonorthogonal pressure corrections."));
        auto initial = integer("initialNonOrthogonalCorrections", 1, 0, 20, "Nonorthogonal corrections for the initial projection.");
        initial.default_value.reset(); initial.default_factory = initialCorrections;
        initial.default_explanation = "nonOrthogonalCorrections";
        result.parameters.push_back(initial);
        result.parameters.push_back(fraction("velocityRelaxation", 1.0, "Momentum underrelaxation."));
        result.parameters.push_back(integer("initializePotentialFlow", 0, 0, 1, "Enable initial potential-flow projection."));
        result.parameters.push_back(number("continuityTolerance", "1", "Normalized continuity tolerance.", 1e-8, 0.0, true, ParameterDictionary::Solution));
        result.parameters.push_back(number("velocityTolerance", "1", "Relative velocity change tolerance.", 1e-7, 0.0, true, ParameterDictionary::Solution));
        result.parameters.push_back(number("momentumTolerance", "1", "Normalized momentum residual tolerance.", 1e-6, 0.0, true, ParameterDictionary::Solution));
        result.parameters.push_back(number("pressureCorrectionTolerance", "1", "Pressure residual tolerance.", 1e-6, 0.0, true, ParameterDictionary::Solution));
        auto coupling = word("turbulenceCoupling", "iterated", "Turbulence coupling strategy.", ParameterDictionary::Solution);
        coupling.choices = {"iterated", "segregated"}; result.parameters.push_back(coupling);
        result.parameters.push_back(integer("turbulenceMaxIterations", 1000, 1, std::numeric_limits<int>::max(), "Maximum turbulence Picard corrections."));
        auto corrections = integer("turbulenceCorrectors", 1, 1, std::numeric_limits<int>::max(), "Model/flow-dependent correction count; host solver computes the effective default.");
        corrections.default_value.reset(); result.parameters.push_back(corrections);
        result.parameter_constraints = {
            {ParameterDictionary::Solution, "maxCorrectors", ParameterRelation::GreaterEqual, "nCorrectors", "The maximum must cover the minimum correction count."},
            {ParameterDictionary::Solution, "turbulenceCorrectors", ParameterRelation::LessEqual, "turbulenceMaxIterations", "Turbulence corrections cannot exceed their limit."}};
        result.limitations.push_back("Segregated turbulence requires maxIterations=1, turbulenceRelaxation=1 when a model is active, and one transport correction; these conditional rules remain checked by the host algorithm.");
        return result;
    }();
    return description;
}
const ModuleDescription& spalartAllmaras() {
    static const ModuleDescription description = [] {
        auto result = turbulence("SA", "Positive-variable Spalart-Allmaras with ft2 and supplied geometric wall distance.");
        result.fields.push_back(scalar("nuTilda", "m^2/s", "SA working variable."));
        result.fields.push_back(scalar("wallDistance", "m", "Distance to the nearest wall."));
        for (const auto& parameter : {
            number("saNuTildaMin", "m^2/s", "Positive lower bound for nuTilda.", 1e-14, 0.0, true),
            number("saWallDistanceMin", "m", "Positive lower bound for wall distance.", 1e-12, 0.0, true),
            number("saCb1", "1", "SA cb1 constant.", 0.1355, 0.0, true),
            number("saCb2", "1", "SA cb2 constant.", 0.622, 0.0, true),
            number("saSigma", "1", "SA sigma constant.", 2.0/3.0, 0.0, true),
            number("saKappa", "1", "von Karman constant.", 0.41, 0.0, true),
            number("saCw2", "1", "SA cw2 constant.", 0.3, 0.0, true),
            number("saCw3", "1", "SA cw3 constant.", 2.0, 0.0, true),
            number("saCv1", "1", "SA cv1 constant.", 7.1, 0.0, true),
            number("saCt3", "1", "SA ct3 constant.", 1.2, 0.0, true),
            number("saCt4", "1", "SA ct4 constant.", 0.5, 0.0, true)}) result.parameters.push_back(parameter);
        auto cw1 = number("saCw1", "1", "SA cw1 constant.", 1.0, 0.0, true);
        cw1.default_value.reset(); cw1.default_factory = saCw1;
        cw1.default_explanation = "saCb1/(saKappa*saKappa) + (1+saCb2)/saSigma";
        result.parameters.push_back(cw1);
        result.equations = {{"nuTildaTransport", "nuTilda", {"convection", "diffusion"}, "SA transport equation."}};
        result.limitations.push_back("Requires supplied wallDistance; no primary trip term and no negative-variable SA formulation.");
        return result;
    }();
    return description;
}
const ModuleDescription& kOmega() {
    static const ModuleDescription description = [] {
        auto result = turbulence("kOmega", "Wilcox 1988 k-omega eddy-viscosity transport.");
        result.fields.push_back(scalar("k", "m^2/s^2", "Turbulent kinetic energy."));
        result.fields.push_back(scalar("omega", "1/s", "Specific dissipation rate."));
        for (const auto& parameter : {
            number("kMin", "m^2/s^2", "Positive lower bound for k.", 1e-12, 0.0, true),
            number("omegaMin", "1/s", "Positive lower bound for omega.", 1e-12, 0.0, true),
            number("maxTurbulentViscosityRatio", "1", "Upper bound on turbulent/molecular viscosity ratio.", 1e5, 0.0, true),
            number("kOmegaBetaStar", "1", "k destruction constant.", 0.09, 0.0, true),
            number("kOmegaBeta", "1", "omega destruction constant.", 3.0/40.0, 0.0, true),
            number("kOmegaGamma", "1", "omega production constant.", 5.0/9.0, 0.0, true),
            number("kOmegaSigmaK", "1", "k diffusion constant.", 0.5, 0.0, true),
            number("kOmegaSigmaOmega", "1", "omega diffusion constant.", 0.5, 0.0, true)}) result.parameters.push_back(parameter);
        result.equations = {{"kTransport", "k", {"convection", "diffusion"}, "Turbulent kinetic energy."},
            {"omegaTransport", "omega", {"convection", "diffusion"}, "Specific dissipation rate."}};
        return result;
    }();
    return description;
}
const ModuleDescription& kEpsilon() {
    static const ModuleDescription description = [] {
        auto result = turbulence("kEpsilon", "Standard k-epsilon eddy-viscosity transport.");
        result.fields.push_back(scalar("k", "m^2/s^2", "Turbulent kinetic energy."));
        result.fields.push_back(scalar("epsilon", "m^2/s^3", "Turbulent dissipation rate."));
        for (const auto& parameter : {
            number("kMin", "m^2/s^2", "Positive lower bound for k.", 1e-12, 0.0, true),
            number("epsilonMin", "m^2/s^3", "Positive lower bound for epsilon.", 1e-12, 0.0, true),
            number("kEpsilonCmu", "1", "Eddy-viscosity constant.", 0.09, 0.0, true),
            number("kEpsilonC1", "1", "epsilon production constant.", 1.44, 0.0, true),
            number("kEpsilonC2", "1", "epsilon destruction constant.", 1.92, 0.0, true),
            number("kEpsilonSigmaK", "1", "k diffusion denominator.", 1.0, 0.0, true),
            number("kEpsilonSigmaEpsilon", "1", "epsilon diffusion denominator.", 1.3, 0.0, true)}) result.parameters.push_back(parameter);
        result.equations = {{"kTransport", "k", {"convection", "diffusion"}, "Turbulent kinetic energy."},
            {"epsilonTransport", "epsilon", {"convection", "diffusion"}, "Turbulent dissipation rate."}};
        return result;
    }();
    return description;
}
namespace {
const ModuleRegistration saDescription(spalartAllmaras());
const ModuleRegistration omegaDescription(kOmega());
const ModuleRegistration epsilonDescription(kEpsilon());
}
} // namespace babelsim::builtin
