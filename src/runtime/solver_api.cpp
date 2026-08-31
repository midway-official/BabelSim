#include "babelsim/runtime.h"
#include "internal/fvm_execution.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace babelsim {

SolveResult solve(const ScalarEquationDefinition& equation, EquationControl control) {
    RunTime& time = RunTime::current();
    const SolveResult result = detail::execution().solve(equation, control);
    if (!result.converged() && time.primary()) {
        std::cerr << "linear solve failed at time=" << time.time()
                  << " iterations=" << result.iterations
                  << " initial=" << result.initial_residual
                  << " final=" << result.final_residual
                  << " relative=" << result.relative_residual
                  << " healthy=" << result.healthy() << '\n';
    }
    return result;
}

namespace {
SolveResult aggregate(const std::array<SolveResult, 3>& components) {
    SolveResult result;
    result.status = SolveStatus::Converged;
    for (const SolveResult& component : components) {
        if (!component.healthy()) result.status = SolveStatus::NumericalFailure;
        else if (!component.converged() && result.status != SolveStatus::NumericalFailure)
            result.status = SolveStatus::MaxIterations;
        result.iterations += component.iterations;
        result.initial_residual = std::hypot(result.initial_residual, component.initial_residual);
        result.final_residual = std::hypot(result.final_residual, component.final_residual);
        result.relative_residual = std::max(result.relative_residual, component.relative_residual);
    }
    if (!result.converged() && RunTime::current().primary())
        std::cerr << "vector equation failed at time=" << RunTime::current().time()
                  << " worst relative residual=" << result.relative_residual << '\n';
    return result;
}

}  // 匿名命名空间

SolveResult solve(const VectorEquationDefinition& equation, EquationControl control) {
    if (control.fix_reference) throw std::invalid_argument("referenceValue requires a scalar equation");
    return aggregate(detail::execution().solve(equation, {control.relaxation, nullptr}));
}

SolveResult solveWithResponse(const VectorEquationDefinition& equation, ScalarField& response,
                              EquationControl control) {
    if (control.fix_reference) throw std::invalid_argument("referenceValue requires a scalar equation");
    return aggregate(detail::execution().solve(equation, {control.relaxation, &response}));
}

namespace fvc {

void evaluate(NormalGradient operation, ScalarField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(ScalarDiffusionFlux operation, ScalarField& result) {
    detail::execution().evaluate(operation, result);
}

void evaluate(ScalarGradient operation, VectorField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(VectorGradient operation, TensorField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(FaceFlux operation, ScalarField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(FaceDivergence operation, ScalarField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(VectorDivergence operation, ScalarField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(ScalarConvection operation, ScalarField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(VectorConvection operation, VectorField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(ScalarInterpolation operation, ScalarField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(VectorInterpolation operation, VectorField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(ScalarReconstruction operation, ScalarField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(VectorReconstruction operation, VectorField& result) {
    detail::execution().evaluate(operation, result);
}
void evaluate(ScalarLaplacian operation, ScalarField& result) {
    detail::execution().evaluate(operation, result);
}

void subtract(
    const ScalarField& coefficient,
    ScalarGradient operation,
    VectorField& target)
{
    detail::execution().subtract(coefficient, operation, target);
}

void subtract(ScalarDiffusionFlux operation, ScalarField& target) {
    detail::execution().subtract(operation, target);
}

}  // fvc 命名空间

namespace diagnostics {

double relativeChange(const VectorField& current, const VectorField& previous) {
    return detail::execution().relativeChange(current, previous);
}

double relativeChange(const ScalarField& current, const ScalarField& previous) {
    return detail::execution().relativeChange(current, previous);
}

double relativeMagnitude(const ScalarField& value, const ScalarField& reference) {
    return detail::execution().relativeMagnitude(value, reference);
}

FluxBalance fluxBalance(const ScalarField& face_flux) {
    return detail::execution().fluxBalance(face_flux);
}

bool all(bool local_condition) {
    return detail::execution().all(local_condition);
}

}  // diagnostics 命名空间

}  // babelsim 命名空间
