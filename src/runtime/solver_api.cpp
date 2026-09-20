#include "babelsim/runtime.h"
#include "internal/fvm_execution.h"

namespace babelsim {

const Methods& numericalMethods() { return RunTime::current().methods(); }
bool primaryProcess() { return RunTime::current().primary(); }

namespace math {

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
void evaluate(TensorDivergence operation, VectorField& result) {
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

void add(FaceFlux operation, ScalarField& target, FaceRegion region) {
    detail::execution().add(operation, target, region);
}

void subtract(ScalarDiffusionFlux operation, ScalarField& target, FaceRegion region) {
    detail::execution().subtract(operation, target, region);
}

}  // math 命名空间

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