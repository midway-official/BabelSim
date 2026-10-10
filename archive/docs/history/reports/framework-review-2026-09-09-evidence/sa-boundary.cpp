#include "physics/RANS/model.h"
#include "physics/simple_common.h"
#include "internal/field_access.h"
#include <iostream>
using namespace babelsim;
int main(){
 Case c("/tmp/babelsim-review-sa-case");
 auto& U=c.vectorField("U",Vec3{}); auto& phi=c.faceFlux("phi",U); auto& mu=c.scalarField("muEffective",0.0);
 auto* m=rans::create(c,U,phi,mu,c.physics().positive("density"),c.physics().positive("dynamicViscosity"));
 ScalarField face(c.mesh(),FieldLocation::Face,"face"); math::evaluate(math::interpolate(mu),face);
 for(Index f=0;f<c.mesh().faceCount();++f) if(c.mesh().boundaryFace(f)) {std::cout<<"muEffectiveBoundary="<<detail::fieldData(face)[f]<<" expected=0.001 because nuTilda_wall=0\n";break;}
 rans::destroy(m);
}
