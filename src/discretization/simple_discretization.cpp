#include "babelsim/simple.h"

#include "internal/scalar_equation_control.h"
#include "internal/simple_discretization.h"
#include "internal/vector_equation_control.h"

namespace babelsim {

namespace simple {

void initializeFaceFlux(const VectorField& velocity, ScalarField& face_flux) {
    fvc::evaluate(fvc::flux(velocity), face_flux);
}

std::array<SolveResult, 3> solveMomentumEquation(
    const VectorEquationDefinition& equation,
    double relaxation,
    ScalarField& rAU)
{
    VectorEquationControl control;
    control.relaxation = relaxation;
    control.mobility = &rAU;
    return detail::solve(equation, control);
}

SolveResult solvePressureCorrectionEquation(
    const ScalarEquationDefinition& equation,
    bool fix_reference)
{
    ScalarEquationControl control;
    control.fix_reference = fix_reference;
    return detail::solve(equation, control);
}

}  // simple 命名空间

void SimpleSolver::NumericalWorkspace::predictMomentumFlux(
    const Methods& methods,
    const VectorField& U,
    const ScalarField& p,
    const ScalarField& rAU,
    ScalarField& phiHbyA)
{
    fvc::evaluate(fvc::grad(p), grad_p);
    fvc::evaluate(fvc::flux(U), phiHbyA);
    rAU_grad_p.assignProduct(rAU, grad_p);
    fvc::evaluate(fvc::interpolate(rAU_grad_p), rAU_grad_p_face);
    fvc::evaluate(fvc::interpolate(rAU), rAU_face);
    detail::applyMomentumInterpolation(
        p, grad_p, rAU_face, rAU_grad_p_face, phiHbyA,
        methods.diffusionFor(p.name()));
}

void SimpleSolver::NumericalWorkspace::preparePressureEquation(
    const ScalarField& phiHbyA)
{
    fvc::evaluate(fvc::div(phiHbyA), div_phiHbyA);
}

}  // babelsim 命名空间
