#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/equ.h"
#include "babelsim/geometry.h"
#include "babelsim/runtime.h"
#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    const Mesh mesh = makeHexBox({2, 1, 1}, {0, 0, 0}, {2, 1, 1});
    RuntimeControl control;
    control.time.delta_t = 0.1;
    control.time.end_time = 0.1;
    control.methods.time = TimeMethod::Euler;
    RunTime run_time = RunTime::forMesh(mesh, control);
    VectorField U(mesh, FieldLocation::Cell, "U");
    VectorField force(mesh, FieldLocation::Cell, "force");
    detail::fieldData(force)[0] = {1, 2, 3};
    detail::fieldData(force)[1] = {3, 2, 1};
    ScalarField response(mesh, FieldLocation::Cell, "response");
    require(run_time.loop(), "time step did not begin");
    time::History<Vec3> velocity_history = time::history(U);
    velocity_history.save(U, run_time.deltaT());
    equ::Equation<Vec3> momentum = equ::createEquation(U);
    equ::ddt(momentum, 1.0, velocity_history);
    equ::source(momentum, force, 2.0);
    require(equ::solve(momentum).converged(), "vector Field source did not converge");
    // 对角响应 V/aP（SIMPLE 的 rAU）由已装配方程显式给出，不需要专门的求解入口。
    response = geometry::cellVolumes(mesh) / momentum.diagonal();
    for (Index cell : detail::meshData(mesh).owned_cells) {
        require(near(detail::fieldData(U)[cell], 0.2 * detail::fieldData(force)[cell], 1e-12), "vector source components are incorrect");
        require(near(detail::fieldData(response)[cell], 0.1, 1e-12), "diagonal response is incorrect");
    }

    // 具有常数零空间的 Neumann 方程由参考值确定唯一解。
    ScalarField p(mesh, FieldLocation::Cell, "p");
    equ::Equation<double> poisson = equ::createEquation(p);
    equ::laplacian(poisson, 1.0, -1.0);  // 左端 -div(grad p)
    poisson.reference(0, 3.0);
    require(equ::solve(poisson).converged(), "reference-constrained Poisson equation did not converge");
    for (Index cell : detail::meshData(mesh).owned_cells) require(near(detail::fieldData(p)[cell], 3.0, 1e-10), "reference value was ignored");
    ScalarField normal(mesh, FieldLocation::Face, "normal");
    ScalarField flux(mesh, FieldLocation::Face, "flux");
    ScalarField coefficient(mesh, FieldLocation::Cell, "k", 2.0);
    math::evaluate(math::normalGradient(p), normal);
    const auto pGradient = math::grad(p);
    ScalarField normalWithGradient(mesh, FieldLocation::Face, "normalWithGradient");
    math::evaluate(math::normalGradient(p, pGradient), normalWithGradient);
    math::evaluate(math::flux(coefficient, p), flux);
    for (Index face : detail::meshData(mesh).owned_faces) {
        require(near(detail::fieldData(normalWithGradient)[face], detail::fieldData(normal)[face], 1e-12),
                "supplied gradient changed normal-gradient semantics");
        require(near(detail::fieldData(flux)[face], 2.0 * mesh.faceArea(face) * detail::fieldData(normal)[face], 1e-10),
                "normal gradient and diffusion flux disagree");
    }
    bool rejected = false;
    try { momentum.reference(0, 0.0); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "vector equation silently ignored a scalar reference constraint");
    const auto expectedLaplacian = math::laplacian(p);
    math::evaluate(math::laplacian(p), p);
    require(near(math::normL2(p-expectedLaplacian),0), "eager in-place laplacian changed result");
    const auto expectedFlux = math::flux(flux,p);
    math::evaluate(math::flux(flux,p),flux);
    require(near(math::normL2(flux-expectedFlux),0), "eager flux coefficient alias changed result");
    rejected = false;
    try { equ::relax(poisson, p, 0.0); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "invalid relaxation was accepted");
    // 一个方程只能有一个边界通量上下文：第二个不同通量必须报错，而不是让边界迹
    // 静默取用最后一个通量。
    rejected = false;
    try {
        ScalarField transported(mesh, FieldLocation::Cell, "transported");
        ScalarField other_flux(mesh, FieldLocation::Face, "otherFlux");
        equ::Equation<double> equation = equ::createEquation(transported);
        equ::div(equation, flux);
        equ::div(equation, other_flux);
    } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "a second boundary flux was silently accepted");
    // 欠松弛只改变迭代路径；同一时间步的旧时间层不随重复求解更新。
    ScalarField T(mesh, FieldLocation::Cell, "T");
    time::History<double> temperature_history = time::history(T);
    temperature_history.save(T, run_time.deltaT());
    equ::Equation<double> energy = equ::createEquation(T);
    equ::ddt(energy, 1.0, temperature_history);
    equ::source(energy, 2.0);
    equ::relax(energy, T, 0.5);
    require(equ::solve(energy).converged(), "scalar relaxation failed");
    for (Index cell : detail::meshData(mesh).owned_cells)
        require(near(detail::fieldData(T)[cell], 0.1, 1e-12), "scalar relaxation formula changed");
    std::cout << "public_equation_test: vector source, response, reference and face operators passed\n";
}