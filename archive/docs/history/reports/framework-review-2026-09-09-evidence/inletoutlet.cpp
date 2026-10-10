#include "babelsim/runtime.h"
#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "test_util.h"
#include <iostream>
using namespace babelsim;
int main(){
 for(bool mixed : {false,true}) {
  const Mesh mesh=makeHexBox({1,1,1},{0,0,0},{1,1,1});
  ScalarField C(mesh,FieldLocation::Cell,"C",0.0), phi(mesh,FieldLocation::Face,"phi",0.0);
  C.setBoundary(0,mixed ? BoundaryCondition<double>::inletOutlet(1.0) : fixedValue(1.0));
  C.setBoundary(1,fixedValue(0.0));
  detail::fieldData(phi)[detail::meshData(mesh).patches[0].faces.front()]=-1;
  detail::fieldData(phi)[detail::meshData(mesh).patches[1].faces.front()]=1;
  RuntimeControl rc;rc.methods.diffusion=DiffusionMethod::Orthogonal;rc.methods.convection=ConvectionMethod::Upwind;
  RunTime time=RunTime::forMesh(mesh,rc);
  auto r=solve(eqn::div(phi,C)==eqn::laplacian(1.0,C));
  std::cout<<(mixed?"inletOutlet":"fixedValue")<<" converged="<<r.converged()<<" C="<<detail::fieldData(C)[0]<<" expectedInflowSolution=0.6\n";
 }
}
