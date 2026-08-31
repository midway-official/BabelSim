#include "babelsim/incompressible.h"

#include "internal/incompressible_discretization.h"

namespace babelsim {

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
        methods.diffusion);
}

void SimpleSolver::NumericalWorkspace::preparePressureEquation(
    const ScalarField& phiHbyA)
{
    fvc::evaluate(fvc::div(phiHbyA), div_phiHbyA);
}

}  // babelsim 命名空间
