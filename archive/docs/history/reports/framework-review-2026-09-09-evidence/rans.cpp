#include "physics/simple_common.h"
#include "physics/RANS/model.h"
#include "internal/field_access.h"
#include <iostream>
using namespace babelsim;
int main(){
 Case c("/tmp/babelsim-review-rans-case");
 auto& U=c.vectorField("U",Vec3{}); auto& phi=c.faceFlux("phi",U);
 auto& mu=c.scalarField("muEffective",0.0);
 const double rho=c.physics().positive("density"), molecular=c.physics().positive("dynamicViscosity");
 auto* m=rans::create(c,U,phi,mu,rho,molecular);
 c.start();
 for(int i=0;i<5;++i){ auto r=rans::correct(*m); std::cout << "iteration="<<i<<" converged="<<r.converged()<<" dTurb="<<rans::relativeChange(*m)<<" k="<<detail::fieldData(c.scalarField("k"))[0]<<" epsilon="<<detail::fieldData(c.scalarField("epsilon"))[0]<<'\n'; }
 auto& k=c.scalarField("k");
 ScalarField src(c.mesh(),FieldLocation::Cell,"source",-0.001);
 ScalarField gamma(c.mesh(),FieldLocation::Cell,"gamma",0.00109);
 ScalarField lap(c.mesh(),FieldLocation::Cell,"lap");
 math::evaluate(math::laplacian(gamma,k),lap);
 std::cout << "post-clipping k PDE residual per volume=" << detail::fieldData(lap)[0] -0.001 << '\n';
 rans::destroy(m);
}
