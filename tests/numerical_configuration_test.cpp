#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/geometry.h"
#include "babelsim/numerics_io.h"
#include "babelsim/operators.h"
#include "babelsim/runtime.h"
#include "internal/field_access.h"
#include "test_util.h"
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace babelsim;
namespace {
template<class F> void rejects(F f, const std::string& message) {
    bool rejected = false;
    try { f(); } catch (const std::exception&) { rejected = true; }
    require(rejected, message);
}
void writeFile(const std::filesystem::path& path, const std::string& text) {
    std::ofstream file(path); file << text;
}
Mesh warpedMesh() {
    std::vector<Vec3> points;
    for (int z=0; z<=1; ++z) for (int y=0; y<=4; ++y) for (int x=0; x<=4; ++x) {
        const double a=x/4.0, b=y/4.0;
        points.push_back({a+0.12*std::sin(3.141592653589793*a)*std::sin(3.141592653589793*b),
            b+0.06*std::sin(6.283185307179586*a)*std::sin(3.141592653589793*b), double(z)});
    }
    return makeHexFromVertices({4,4,1}, std::move(points));
}
}
int main() {
    const auto root=std::filesystem::temp_directory_path()/"babelsim-numerical-configuration-test";
    std::filesystem::remove_all(root);
    std::filesystem::copy("cases/heat",root,std::filesystem::copy_options::recursive);
    const std::string defaults="time euler\n"
        "equation.first.interpolation linear\n"
        "equation.first.gradient greenGauss\n"
        "equation.first.convection linearUpwind\n"
        "equation.first.diffusion orthogonal\n"
        "equation.first.term.advection.convection upwind\n"
        "equation.first.coefficientInterpolation corrected\n"
        "equation.first.coefficientGradient leastSquares\n"
        "equation.second.interpolation linear\n"
        "equation.second.gradient leastSquares\n"
        "equation.second.convection central\n"
        "equation.second.diffusion orthogonal\n"
        "operation.pressureGradient.interpolation linear\n"
        "operation.pressureGradient.gradient leastSquares\n"
        "operation.pressureGradient.convection upwind\n"
        "operation.pressureGradient.diffusion orthogonal\n";
    writeFile(root/"numerics/methods.bs", defaults);
    writeFile(root/"numerics/solution.bs",
        "equation.first.solver bicgstab\n"
        "equation.first.preconditioner ilut\n"
        "equation.first.absoluteTolerance 1e-13\n"
        "equation.first.relativeTolerance 1e-7\n"
        "equation.first.maxIterations 300\n"
        "equation.first.ilutFillFactor 3\n"
        "equation.second.solver cg\n"
        "equation.second.preconditioner incompleteCholesky\n"
        "equation.second.absoluteTolerance 1e-14\n"
        "equation.second.relativeTolerance 1e-10\n"
        "equation.second.maxIterations 1000\n");
    {
        Case problem(root);
        auto& x=problem.scalarField("T");
        const auto first=readEquationControl(problem,"first",x,{"advection","diffusion"});
        const auto second=readEquationControl(problem,"second",x);
        require(first.spatial.convection==ConvectionMethod::LinearUpwind &&
            first.options("advection").convection==ConvectionMethod::Upwind &&
            second.spatial.convection==ConvectionMethod::Central,"field/equation/term precedence");
        require(first.linear.max_iterations==300 && first.linear.ilut_fill_factor==3 &&
            first.linear.absolute_tolerance==1e-13 && first.linear.relative_tolerance==1e-7,
            "partial solver inheritance lost a parent member");
        require(first.spatial.coefficientGradient==GradientMethod::LeastSquares &&
            first.spatial.gradient==GradientMethod::GreenGauss,"coefficient and unknown gradients were coupled");
        auto a=equ::createEquation(x,first), b=equ::createEquation(x,second);
        const auto op=problem.methods().operationOptions("pressureGradient");
        require(op.gradient==GradientMethod::LeastSquares,"named math options");
        problem.methods().requireAllUsed(); problem.solution().requireAllUsed();
        ScalarField phi(x.mesh(),FieldLocation::Face,"phi",1.0);
        equ::div(a,phi,1.0,"advection"); equ::div(b,phi);
        require(math::normL2(a.diagonal()-b.diagonal())>1e-6,"equation choices did not change actual convection assembly");
        b.reset();
        OperatorOptions explicitOptions; explicitOptions.convection=ConvectionMethod::Upwind;
        equ::div(b,phi,1.0,{},explicitOptions);
        require(math::normL2(a.diagonal()-b.diagonal())<1e-12,"explicit operator override did not win");
        auto copy=a.copy(); a.reset();
        require(copy.name()=="first" && a.options("advection").convection==ConvectionMethod::Upwind,
            "copy/reset discarded equation identity or options");
        rejects([&]{equ::div(a,phi,1.0,"typo");},"unknown term silently fell back");
    }
    for (const auto& bad : {"interpolation linear\n", "equation.first.time bdf2\n", "equation.first.gradent leastSquares\n",
                           "equation.first.term..gradient leastSquares\n",
                           "equation.first.gradient greenGauss\nequation.first.gradient leastSquares\n",
                           "equation.first.gradient green_gauss\n",
                           "equation.first.convection secondOrderUpwind\n",
                           "equation.first.diffusion limited_corrected\n"}) {
        writeFile(root/"numerics/methods.bs",defaults+bad);
        rejects([&]{readMethodsFile(root/"numerics/methods.bs");},"invalid named numerical configuration accepted");
    }
    writeFile(root/"numerics/methods.bs",defaults+"equation.typo.gradient leastSquares\n");
    const auto unused=readMethodsFile(root/"numerics/methods.bs");
    rejects([&]{unused.requireAllUsed();},"unused equation selector accepted");
    // A fully named scalar profile needs neither scalarSolver nor vectorSolver.
    writeFile(root/"numerics/methods.bs",
        "time euler\n"
        "equation.only.interpolation linear\n"
        "equation.only.gradient greenGauss\n"
        "equation.only.convection upwind\n"
        "equation.only.diffusion orthogonal\n");
    writeFile(root/"numerics/solution.bs","equation.only.solver bicgstab\nequation.only.preconditioner ilut\n"
        "equation.only.absoluteTolerance 1e-14\nequation.only.relativeTolerance 1e-10\nequation.only.maxIterations 100\n");
    {
        Case problem(root); auto& x=problem.scalarField("T");
        auto a=equ::createEquation(problem,"only",x);
        equ::reaction(a,1.0); equ::source(a,2.0);
        require(equ::solve(a).converged(),"fully named scalar solve without unused type defaults");
        problem.solution().requireAllUsed();
    }
    std::filesystem::remove_all(root);
    {
        auto mesh=warpedMesh(); auto runtime=RunTime::forMesh(mesh);
        ScalarField x(mesh,FieldLocation::Cell,"x"), k(mesh,FieldLocation::Cell,"k");
        for (Index cell=0; cell<mesh.cellCount(); ++cell) {
            const auto c=mesh.cellCentre(cell);
            detail::fieldData(x)[cell]=c.x*c.x+0.7*c.y*c.y+c.x*c.y;
            detail::fieldData(k)[cell]=1+c.x*c.x+2*c.y*c.y;
        }
        OperatorOptions linear;
        linear.interpolation = InterpolationMethod::Linear;
        linear.gradient=GradientMethod::LeastSquares;
        linear.convection = ConvectionMethod::Upwind;
        linear.diffusion=DiffusionMethod::Corrected;
        linear.coefficientInterpolation=InterpolationMethod::Linear;
        linear.coefficientGradient=GradientMethod::GreenGauss;
        auto corrected=linear; corrected.coefficientInterpolation=InterpolationMethod::Corrected;
        corrected.coefficientGradient=GradientMethod::LeastSquares;
        const auto defaultGradient=math::grad(x, linear);
        OperatorOptions gg = linear; gg.gradient=GradientMethod::GreenGauss;
        const auto differentGradient=math::grad(x,gg);
        require(math::normL2(defaultGradient-differentGradient)>1e-5,"gradient test did not exercise distinct methods");
        rejects([&]{ (void)math::grad(x); }, "math operation accepted an incomplete implicit configuration");
        auto a=testEquation(x, "coefficientLinear");
        equ::laplacian(a,k,-1,{},linear);
        const auto firstDiagonal=a.diagonal();
        auto b=testEquation(x, "coefficientCorrected");
        equ::laplacian(b,k,-1,{},corrected);
        require(math::normL2(firstDiagonal-b.diagonal())>1e-5,"coefficient interpolation did not change matrix");
        const auto volume=geometry::cellVolumes(mesh);
        require(math::normL2(equ::residual(a,x)-volume*math::laplacian(k,x,linear))<1e-10,
            "explicit laplacian and implicit assembly selected different coefficient/gradient options");
        require(math::normL2(equ::faceFlux(a,x)+math::flux(k,x,linear))<1e-10,
            "assembled and explicit diffusion flux disagree");
        equ::add(a,b,0.3);
        const auto frozenFlux=equ::faceFlux(a,x);
        require(math::normL2(equ::apply(a,x)-a.rhs()-volume*math::div(frozenFlux))<1e-10,
            "mixed diffusion contributions violate matrix/face-flux conservation");
        k.fill(99);
        require(math::normL2(frozenFlux-equ::faceFlux(a,x))<1e-13,"face flux reread live coefficient instead of snapshot");
        const auto copied=a.copy(); a.reset();
        require(math::normL2(frozenFlux-equ::faceFlux(copied,x))<1e-13,"copy lost frozen diffusion settings");
    }
    std::cout << "numerical_configuration_test: hierarchy, named solvers, explicit math, independent coefficient reconstruction and frozen flux passed\n";
}
