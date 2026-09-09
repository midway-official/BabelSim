#include "physics/transient_simple/algorithm.h"
#include "babelsim/runtime.h"
#include "internal/field_access.h"
#include "test_util.h"
#include <iostream>
using namespace babelsim;
int main(){
 for(auto method : {TimeMethod::Euler,TimeMethod::BDF2}) {
  const Mesh mesh=makeHexBox({1,1,1},{0,0,0},{1,1,1});
  IncompressibleFields f(mesh); f.velocity.fill({1,0,0});
  for(Index p=0;p<mesh.patchCount();++p) { f.velocity.setBoundary(p,fixedValue(Vec3{})); f.pressure.setBoundary(p,fixedValue(0.0)); }
  RuntimeControl rc; rc.methods.time=method; rc.methods.diffusion=DiffusionMethod::Orthogonal;
  rc.time={0,0.3,0.1}; rc.vector_solver.absolute_tolerance=1e-14; rc.vector_solver.relative_tolerance=1e-12;
  RunTime time=RunTime::forMesh(mesh,rc);
  SimpleControl control; control.velocity_tolerance=1e-11; control.max_iterations=1000;
  TransientSimpleAlgorithm solver(f,{1,0.1},control);
  double prev=1,older=1,error=0;
  while(time.loop()) {
   solver.beginTimeStep(); while(solver.loop()) solver.iterate();
   require(solver.converged(),"transient SIMPLE failed");
   const double exact=method==TimeMethod::Euler || time.step()==1 ? prev/1.12 : (2*prev-.5*older)/1.62;
   error=std::max(error,std::abs(detail::fieldData(f.velocity)[0].x-exact)); older=prev; prev=exact;
  }
  std::cout<<"method="<<static_cast<int>(method)<<" steps="<<time.step()<<" maxDiscreteODEError="<<error<<'\n';
 }
}
