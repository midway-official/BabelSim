#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/runtime.h"
#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    const Mesh mesh = Mesh::cartesian({2, 1, 1}, {0, 0, 0}, {2, 1, 1});
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
    require(solveWithResponse(eqn::ddt(U) == eqn::source(2.0, force), response).converged(),
            "vector Field source did not converge");
    for (Index cell : detail::meshData(mesh).owned_cells) {
        require(near(detail::fieldData(U)[cell], 0.2 * detail::fieldData(force)[cell], 1e-12), "vector source components are incorrect");
        require(near(detail::fieldData(response)[cell], 0.1, 1e-12), "diagonal response is incorrect");
    }

    // 具有常数零空间的 Neumann 方程由参考值确定唯一解。
    ScalarField p(mesh, FieldLocation::Cell, "p");
    require(solve(-eqn::laplacian(1.0, p) == eqn::source(0.0), referenceValue(3.0)).converged(),
            "reference-constrained Poisson equation did not converge");
    for (Index cell : detail::meshData(mesh).owned_cells) require(near(detail::fieldData(p)[cell], 3.0, 1e-10), "reference value was ignored");
    ScalarField normal(mesh, FieldLocation::Face, "normal");
    ScalarField flux(mesh, FieldLocation::Face, "flux");
    ScalarField coefficient(mesh, FieldLocation::Cell, "k", 2.0);
    math::evaluate(math::normalGradient(p), normal);
    math::evaluate(math::flux(coefficient, p), flux);
    for (Index face : detail::meshData(mesh).owned_faces)
        require(near(detail::fieldData(flux)[face], 2.0 * mesh.faceArea(face) * detail::fieldData(normal)[face], 1e-10),
                "normal gradient and diffusion flux disagree");
    bool rejected = false;
    try { (void)solve(eqn::ddt(U) == eqn::source(force), referenceValue(0.0)); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "vector equation silently ignored a scalar reference constraint");
    rejected = false;
    try { math::evaluate(math::laplacian(p), p); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "in-place explicit laplacian corrupted its input");
    rejected = false;
    try { math::evaluate(math::flux(flux, p), flux); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "diffusion flux overwrote its coefficient");
    rejected = false;
    try { (void)solve(eqn::ddt(p) == 0.0, relaxed(0.0)); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "invalid relaxation was accepted");
    rejected = false;
    try { (void)solveWithResponse(eqn::ddt(response, U) == eqn::source(force), response); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "momentum response overwrote an equation coefficient");
    // 欠松弛只改变迭代路径；同一时间步的旧时间层不随重复求解更新。
    ScalarField T(mesh, FieldLocation::Cell, "T");
    require(solve(eqn::ddt(T) == 2.0, relaxed(0.5)).converged(), "scalar relaxation failed");
    for (Index cell : detail::meshData(mesh).owned_cells)
        require(near(detail::fieldData(T)[cell], 0.1, 1e-12), "scalar relaxation formula changed");
    std::cout << "public_equation_test: vector source, response, reference and face operators passed\n";
}
