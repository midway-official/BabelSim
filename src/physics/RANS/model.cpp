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
    const ConfigLine& entry = settings.entry("turbulenceModel");
    if (entry.tokens.size() != 2) {
        throw std::invalid_argument("turbulenceModel expects exactly one model name");
    }
    std::string result;
    for (unsigned char character : entry.tokens[1]) {
        if (std::isalnum(character)) result += static_cast<char>(std::tolower(character));
    }
    return result;
}

}  // 匿名命名空间

Model::Model(
    Case& problem,
    const VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity)
    : m_problem(problem),
      m_velocity(velocity),
      m_face_flux(face_flux),
      m_effective_viscosity(effective_viscosity),
      m_density(density),
      m_molecular_viscosity(molecular_viscosity),
      m_relaxation(positiveSetting(problem.physics(), "turbulenceRelaxation", 0.7)),
      m_tolerance(positiveSetting(problem.physics(), "turbulenceTolerance", 1e-6)),
      m_velocity_gradient(problem.tensorField("ransGradU", Tensor3{})),
      m_strain_measure(problem.scalarField("ransStrain2", 0.0))
{
    if (m_relaxation > 1.0) {
        throw std::invalid_argument("turbulenceRelaxation must not exceed one");
    }
}

SolveResult Model::solveTransport(
    ScalarField& variable,
    const ScalarField& diffusivity,
    const ScalarField& source) const
{
    if (numericalMethods().time == TimeMethod::Steady) {
        return solve(
            eqn::div(m_density, m_face_flux, variable) ==
                eqn::laplacian(diffusivity, variable) + eqn::source(source),
            relaxed(m_relaxation));
    }
    return solve(
        eqn::ddt(m_density, variable) +
            eqn::div(m_density, m_face_flux, variable) ==
            eqn::laplacian(diffusivity, variable) + eqn::source(source),
        relaxed(m_relaxation));
}

void Model::updateKinematics() {
    math::evaluate(math::grad(m_velocity), m_velocity_gradient);
    m_strain_measure.evaluate(m_velocity_gradient, strainMeasure);
}

void Model::setEddyViscosity(const ScalarField& turbulent_viscosity) {
    m_effective_viscosity.fill(m_molecular_viscosity);
    m_effective_viscosity.addScaled(1.0, turbulent_viscosity);
}

SolveResult Model::combine(const SolveResult& first, const SolveResult& second) {
    SolveResult result;
    result.status = first.status == SolveStatus::NumericalFailure ||
            second.status == SolveStatus::NumericalFailure
        ? SolveStatus::NumericalFailure
        : first.converged() && second.converged()
            ? SolveStatus::Converged : SolveStatus::MaxIterations;
    result.iterations = first.iterations + second.iterations;
    result.initial_residual = std::hypot(first.initial_residual, second.initial_residual);
    result.final_residual = std::hypot(first.final_residual, second.final_residual);
    result.relative_residual = std::max(first.relative_residual, second.relative_residual);
    return result;
}

double positiveSetting(const Parameters& settings, const char* key, double fallback) {
    const double value = settings.number(key, fallback);
    if (!(value > 0.0) || !std::isfinite(value)) {
        throw std::invalid_argument(std::string(key) + " must be positive");
    }
    return value;
}

double strainMeasure(const Tensor3& gradient) {
    const double divergence = gradient[0][0] + gradient[1][1] + gradient[2][2];
    double squared = 0.0;
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            double value = 0.5 * (gradient[row][column] + gradient[column][row]);
            if (row == column) value -= divergence / 3.0;
            squared += value * value;
        }
    }
    return 2.0 * squared;
}

double vorticityMagnitude(const Tensor3& gradient) {
    double squared = 0.0;
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            const double value = 0.5 *
                (gradient[row][column] - gradient[column][row]);
            squared += value * value;
        }
    }
    return std::sqrt(2.0 * squared);
}

Model* create(
    Case& problem,
    const VectorField& velocity,
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

void destroy(Model* model) noexcept { delete model; }
SolveResult correct(Model& model) { return model.correct(); }
double relativeChange(const Model& model) { return model.relativeChange(); }
double tolerance(const Model& model) { return model.tolerance(); }
const char* name(const Model& model) { return model.modelName(); }

}  // babelsim::rans 命名空间
