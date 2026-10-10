#include "babelsim/equ.h"
#include "babelsim/geometry.h"
#include "babelsim/runtime.h"
#include "internal/field_access.h"
#include "test_util.h"
#include <iostream>
#include <iomanip>
using namespace babelsim;
int main() {
 std::cout << std::setprecision(16);
 {
 auto mesh=makeHexBox({2,1,1},{0,0,0},{2,1,1});
 auto runtime=RunTime::forMesh(mesh);
 auto options=testEquationControl("probe").spatial;
 VectorField u(mesh,FieldLocation::Cell,"U",Vec3{1,0,0});
 u.boundary("plus_x")=BoundaryCondition<Vec3>::inletOutlet(Vec3{-2,0,0});
 ScalarField previousFlux(mesh,FieldLocation::Face,"phi",-1);
 u.setBoundaryFlux(previousFlux);
 auto before=math::interpolate(u,options);
 auto generated=math::flux(u,options);
 auto after=math::interpolate(u,options);
 for(Index f=0;f<mesh.faceCount();++f) if(mesh.boundaryFace(f)&&mesh.patchName(mesh.boundaryPatch(f))=="plus_x")
 std::cout << "inletOutlet: trace_before="<<detail::fieldData(before)[f].x<<" returned_flux="<<detail::fieldData(generated)[f]<<" trace_after="<<detail::fieldData(after)[f].x<<" expected_flux=-2\n";
 }
 {
 std::vector<Vec3> vertices;
 for(int k=0;k<=1;++k) for(int j=0;j<=3;++j) for(int i=0;i<=3;++i) vertices.push_back({i+0.7*j,double(j),double(k)});
 auto mesh=makeHexFromVertices({3,3,1},vertices);
 auto runtime=RunTime::forMesh(mesh);
 for(auto method:{DiffusionMethod::Corrected,DiffusionMethod::LimitedCorrected}) {
 ScalarField p(mesh,FieldLocation::Cell,"p");
 for(Index patch=0;patch<mesh.patchCount();++patch)p.setBoundary(patch,fixedValue(0.0));
 p.evaluate([](Vec3 x){return 1+x.x*x.y;});
 auto control=testEquationControl("p");control.spatial.diffusion=method;control.linear.relative_tolerance=1e-13;control.linear.absolute_tolerance=1e-14;
 auto a=equ::createEquation(p,control); equ::laplacian(a,1.0,-1);
 auto V=geometry::cellVolumes(mesh);
 auto check=[&](const char* stage) {
 auto flux=equ::faceFlux(a,p);
 double defect=math::normL2(V*math::div(flux)-(equ::apply(a,p)-a.rhs()));
 std::cout<<"diffusion_flux: method="<<(method==DiffusionMethod::Corrected?"corrected":"limitedCorrected")<<" stage="<<stage<<" identity_defect="<<defect<<" matrix_residual="<<math::normL2(equ::residual(a,p))<<" flux_divergence="<<math::normL2(V*math::div(flux))<<"\n";
 };
 check("assembly");equ::solve(a);check("solved");
 }
 }
}
