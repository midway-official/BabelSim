#include "babelsim/equ.h"
#include "babelsim/runtime.h"
#include "internal/field_access.h"
#include "test_util.h"
#include <iostream>
#include <iomanip>
using namespace babelsim;
int main(){
 std::cout<<std::setprecision(16);
 {
 std::vector<Vec3> vertices;
 for(int k=0;k<=1;++k)for(int j=0;j<=3;++j)for(int i=0;i<=3;++i)vertices.push_back({i+0.7*j,double(j),double(k)});
 auto mesh=makeHexFromVertices({3,3,1},vertices);
 auto runtime=RunTime::forMesh(mesh);
 for(auto diffusion:{DiffusionMethod::Orthogonal,DiffusionMethod::Corrected}) {
  std::vector<ScalarField> answers;
  for(int steps:{20,40,80,160}){
   double dt=0.2/steps;
   ScalarField x(mesh,FieldLocation::Cell,"T");
   for(Index patch=0;patch<mesh.patchCount();++patch)x.setBoundary(patch,fixedValue(0.0));
   x.evaluate([](Vec3 r){return 1+0.1*r.x*r.y;});
   auto h=time::history(x);auto control=testEquationControl("T");control.spatial.diffusion=diffusion;
   control.linear.relative_tolerance=1e-13;control.linear.absolute_tolerance=1e-14;
   auto eq=equ::createEquation(x,control);
   for(int n=0;n<steps;++n){
    h.save(x,dt);eq.reset();
    equ::ddt(eq,1.,h.previous(),dt,n==0?TimeMethod::Euler:TimeMethod::BDF2,&h.older(),dt);
    equ::laplacian(eq,1.,-1.);require(equ::solve(eq).converged(),"time-order solve");
   }
   answers.push_back(x);
  }
  double e0=math::normL2(answers[0]-answers[1]),e1=math::normL2(answers[1]-answers[2]),e2=math::normL2(answers[2]-answers[3]);
  std::cout<<"BDF2 diffusion="<<(diffusion==DiffusionMethod::Orthogonal?"orthogonal":"corrected")<<" successive_differences="<<e0<<","<<e1<<","<<e2<<" ratios="<<e0/e1<<","<<e1/e2<<"\n";
 }
 }
 {
 auto advMesh=makeHexBox({8,1,1},{0,0,0},{1,1,1});
 auto advRuntime=RunTime::forMesh(advMesh);
 for(auto method:{ConvectionMethod::Upwind,ConvectionMethod::LinearUpwind}) {
  std::vector<ScalarField> answers;
  auto control=testEquationControl("C");control.spatial.convection=method;
  control.linear.relative_tolerance=1e-13;control.linear.absolute_tolerance=1e-14;
  VectorField U(advMesh,FieldLocation::Cell,"U",Vec3{1,0,0});
  auto phi=math::flux(U,control.spatial);
  for(int steps:{20,40,80,160}) {
   double dt=0.2/steps;ScalarField x(advMesh,FieldLocation::Cell,"C");
   x.boundary("minus_x")=fixedValue(0.0);
   x.evaluate([](Vec3 r){return std::sin(3.141592653589793*r.x);});
   auto h=time::history(x);auto eq=equ::createEquation(x,control);
   for(int n=0;n<steps;++n) {
    h.save(x,dt);eq.reset();
    equ::ddt(eq,1.,h.previous(),dt,n==0?TimeMethod::Euler:TimeMethod::BDF2,&h.older(),dt);
    equ::div(eq,phi);require(equ::solve(eq).converged(),"advection time-order solve");
   }
   answers.push_back(x);
  }
  double e0=math::normL2(answers[0]-answers[1]),e1=math::normL2(answers[1]-answers[2]),e2=math::normL2(answers[2]-answers[3]);
  std::cout<<"BDF2 convection="<<(method==ConvectionMethod::Upwind?"upwind":"linearUpwind")<<" successive_differences="<<e0<<","<<e1<<","<<e2<<" ratios="<<e0/e1<<","<<e1/e2<<"\n";
 }
 }

}
