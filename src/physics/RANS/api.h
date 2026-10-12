#pragma once

#include "babelsim/case.h"
#include "babelsim/solver.h"
#include "babelsim/history.h"
#include "babelsim/builtin_descriptions.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace babelsim::rans {

// Per-equation diagnostics retain names and units. Never combine absolute
// residuals of different transported quantities such as k and omega.
struct TransportEquationResult {
    std::string field;
    SolveResult linearSolve;
    double initialResidual; // normalized, before relaxation and solve
    double relativeChange; // after solving and bounding, vs previous outer iterate
};

struct TransportResult {
    std::vector<TransportEquationResult> equations;

    bool healthy() const {
        for (const auto& equation : equations)
            if (!equation.linearSolve.healthy() || !std::isfinite(equation.initialResidual)
                || !std::isfinite(equation.relativeChange)) return false;
        return true;
    }
    bool linearConverged() const {
        for (const auto& equation : equations)
            if (!equation.linearSolve.converged()) return false;
        return true;
    }
    double initialResidual() const {
        double result = 0.0;
        for (const auto& equation : equations) {
            if (!std::isfinite(equation.initialResidual))
                return std::numeric_limits<double>::infinity();
            result = std::max(result, equation.initialResidual);
        }
        return result;
    }
    double relativeChange() const {
        double result = 0.0;
        for (const auto& equation : equations) {
            if (!std::isfinite(equation.relativeChange))
                return std::numeric_limits<double>::infinity();
            result = std::max(result, equation.relativeChange);
        }
        return result;
    }
};

// Interface only: every concrete model owns its fields, closure and transport.
class Model {
public:
    virtual ~Model() = default;
    virtual void registerTimeState(time::HistoryGroup&) = 0;
    virtual void beginStep() {}
    virtual const char* modelName() const = 0;
    virtual TransportResult solveTransport() = 0;
    virtual double tolerance() const = 0;
    virtual double relaxationFactor() const = 0;
};

Model* createModel(Case&, VectorField&, const ScalarField&, ScalarField&,
    const builtin::FluidConfig&, const OperatorOptions&);
void destroy(Model*) noexcept;
using Handle = std::unique_ptr<Model, void(*)(Model*)>;

// SIMPLE-facing coupling: configuration and field ownership, no numerical model.
class Turbulence {
public:
    Turbulence(Case& problem, VectorField& velocity, const ScalarField& phi,
        const builtin::FluidConfig& fluid, const OperatorOptions& momentumOptions);
    explicit operator bool() const { return bool(model_); }
    const ScalarField& effectiveViscosity() const { return *viscosity_; }
    // Explicit remainder of the eddy-viscosity deviatoric stress after the
    // implicit div(muEff*grad(U)) term. Do not add the full stress a second time.
    TensorField deviatoricStressRemainder(const VectorField& velocity) const;
    void registerTimeState(time::HistoryGroup& histories) {
        histories.trackForRollback(*viscosity_);
        if (model_) model_->registerTimeState(histories);
    }
    void beginStep() { if (model_) model_->beginStep(); }
    TransportResult solveTransport() {
        return model_ ? model_->solveTransport() : TransportResult{};
    }
    double relaxation() const { return model_ ? model_->relaxationFactor() : 1.0; }
    double tolerance() const { return model_ ? model_->tolerance() : 0.0; }

private:
    ScalarField* viscosity_;
    OperatorOptions momentumOptions_;
    Handle model_;
};

inline Turbulence create(Case& problem, VectorField& velocity, const ScalarField& phi,
    const builtin::FluidConfig& fluid, const OperatorOptions& momentumOptions) {
    return Turbulence(problem, velocity, phi, fluid, momentumOptions);
}

} // namespace babelsim::rans
