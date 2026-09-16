#include "babelsim/equ.h"
#include "babelsim/geometry.h"
#include "babelsim/runtime.h"
#include "internal/field_access.h"
#include "test_util.h"
#include <iostream>
using namespace babelsim;
int main() {
    // Non-unit cell volume detects accidental normalization or double integration.
    auto mesh=makeHexBox({2,1,1},{0,0,0},{4,1,1});
    auto runtime=RunTime::forMesh(mesh);
    ScalarField x(mesh,FieldLocation::Cell,"x",1.0);
    ScalarField old(x), capacity(mesh,FieldLocation::Cell,"capacity",3.0);
    auto a=equ::createEquation(x);
    equ::ddt(a,capacity,old,0.5);
    equ::source(a,6.0);
    capacity.fill(99); old.fill(99); // Assembled coefficients must be frozen.
    auto d=a.diagonal(); auto b=a.rhs();
    require(near(detail::fieldData(d)[0],12),"integrated time diagonal");
    require(near(detail::fieldData(b)[0],24),"source sign/volume or frozen history");
    require(equ::solve(a,x).converged(),"procedural solve");
    require(near(detail::fieldData(x)[0],2),"frozen system changed during solve");
    auto r=equ::residual(a,x);
    require(near(detail::fieldData(r)[0],0),"b-Ax residual");
    auto saved=a.copy();
    a.reset();
    equ::add(a,saved,2); equ::scale(a,0.5);
    const auto volume = geometry::cellVolumes(mesh);
    auto mobility=volume / a.diagonal();
    require(near(detail::fieldData(mobility)[0],1.0/6),"V/aP response");
    equ::relax(a,x,0.5);
    auto relaxed=volume / a.diagonal();
    require(near(detail::fieldData(relaxed)[0],1.0/12),"relaxed response");
    require(equ::solve(a,x).converged() && near(detail::fieldData(x)[0],2),"relax fixed point");

    // Variable-step BDF2 differentiates a quadratic exactly.
    old.fill(1); ScalarField older(x); older.fill(0);
    a.reset(); equ::ddt(a,1.0,old,0.5,TimeMethod::BDF2,&older,1.0);
    equ::source(a,3.0);
    require(equ::solve(a,x).converged() && near(detail::fieldData(x)[0],2.25),"variable-step BDF2");

    // Known linear diffusion solution verifies boundary RHS and laplacian sign.
    x.fill(0); x.boundary("minus_x")=fixedValue(0.0); x.boundary("plus_x")=fixedValue(4.0);
    a.reset(); equ::laplacian(a,1.0, -1);
    require(equ::solve(a,x).converged(),"diffusion solve");
    require(near(detail::fieldData(x)[0],1) && near(detail::fieldData(x)[1],3),"diffusion boundary/sign");

    // Public sign means the mathematical sign, including boundary RHS terms.
    auto negative = a.copy();
    a.reset(); equ::laplacian(a, 1.0, +1);
    require(near(math::normL2(a.diagonal() + negative.diagonal()), 0), "laplacian LHS sign");
    require(near(math::normL2(a.rhs() + negative.rhs()), 0), "laplacian boundary RHS sign");
    const auto before = x;
    const auto diagonalBefore = a.diagonal();
    const auto rhsBefore = a.rhs();
    require(diagnostics::relativeResidual(a, x) < 1e-12, "assembled residual");
    require(near(math::normL2(x - before), 0), "diagnostic changed unknown");
    require(near(math::normL2(a.diagonal() - diagonalBefore), 0)
        && near(math::normL2(a.rhs() - rhsBefore), 0), "diagnostic changed equation");
    auto independent = x; independent.fill(99);
    require(near(math::normL2(x - before), 0), "field copy aliases original");

    // Correction boundaries inherit homogeneous constraints. An already anchored
    // pressure matrix must not acquire an additional interior reference.
    auto correction = field::homogeneousLike(x, "pPrime");
    require(correction.name() == "pPrime", "homogeneous field name is explicit");
    auto pressure = equ::createEquation(correction);
    equ::laplacian(pressure, 1.0, -1);
    const auto anchored = pressure.diagonal();
    pressure.referenceIfUnanchored(7.0);
    const auto unchanged = pressure.diagonal();
    require(near(detail::fieldData(anchored)[0], detail::fieldData(unchanged)[0]),
            "anchored pressure received a redundant reference");
    correction.fill(2.0);
    require(equ::solve(pressure, correction).converged(), "correction solve");
    require(near(math::normL2(correction), 0.0), "correction boundary is not homogeneous");

    ScalarField freePressure(mesh, FieldLocation::Cell, "freePressure");
    auto freeSystem = equ::createEquation(freePressure);
    equ::laplacian(freeSystem, 1.0, -1);
    freeSystem.referenceIfUnanchored(3.0);
    require(equ::solve(freeSystem, freePressure).converged(), "unanchored pressure solve");
    require(near(detail::fieldData(freePressure)[0], 3.0) &&
            near(detail::fieldData(freePressure)[1], 3.0), "pressure reference value");

    VectorField u(mesh,FieldLocation::Cell,"U");
    auto v=equ::createEquation(u); equ::reaction(v,2.0); equ::source(v,Vec3{2,4,6});
    require(equ::solve(v).converged(),"bound vector solve");
    require(near(detail::fieldData(u)[0],Vec3{1,2,3}),"vector source/diagonal");
    std::cout << "procedural_equation_test: assembly, frozen inputs, algebra, BDF2, diffusion and vector solve passed\n";
}
