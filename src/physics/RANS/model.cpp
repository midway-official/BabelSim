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

void saveOld(Model& model,double dt) { model.saveOld(dt); }
void destroy(Model* model) noexcept { delete model; }
SolveResult correct(Model& model) { return model.correct(); }
double relativeChange(const Model& model) { return model.relativeChange(); }
double tolerance(const Model& model) { return model.tolerance(); }
const char* name(const Model& model) { return model.modelName(); }

double relativeResidual(const Model& model) { return model.relativeResidual(); }

}  // babelsim::rans 命名空间
