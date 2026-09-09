#include "babelsim/runtime.h"
#include "babelsim/mpi_support.h"
#include "babelsim/parallel.h"
#include "physics/simple_common.h"
#include "physics/RANS/model.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"
#include "test_util.h"
#include <iostream>

using namespace babelsim;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    try {
        require(argc == 2, "expected generated model case");
        Case problem(argv[1]);
        auto& U = problem.vectorField("U");
        auto& phi = problem.faceFlux("phi", U);
        const auto& mesh = problem.mesh();
        U.evaluate([](Vec3 x) { return Vec3{2*x.y,0,0}; });
        for (Index p = 0; p < mesh.patchCount(); ++p) U.setBoundary(p, BoundaryCondition<Vec3>::zeroGradient());
        U.setBoundary(2, fixedValue(Vec3{})); U.setBoundary(3, fixedValue(Vec3{2,0,0}));
        math::evaluate(math::flux(U), phi);
        auto& effective = problem.scalarField("effective", 0.0);
        const double rho = problem.physics().positive("density");
        const double mu = problem.physics().positive("dynamicViscosity");
        const auto model_name = problem.physics().entry("turbulenceModel").tokens[1];
        if (model_name == "SA")
            problem.scalarField("nuTilda").setBoundary(0, fixedValue(0.0));
        (void)readSimpleControl(problem.solution());
        std::unique_ptr<rans::Model, void(*)(rans::Model*)> model(
            rans::create(problem, U, phi, effective, rho, mu), rans::destroy);
        require(model != nullptr, "expected turbulence model");
        require(problem.loop(), "expected transient coefficient test step");
        require(model->correct().healthy(), "model coefficient update failed");
        if (model_name == "SA") {
            ScalarField face(mesh, FieldLocation::Face);
            math::evaluate(math::interpolate(effective), face);
            for (Index f : detail::meshData(mesh).owned_faces)
                if (mesh.boundaryFace(f) && mesh.boundaryPatch(f) == 0)
                    require(std::abs(detail::fieldData(face)[f] - mu) < 1e-13,
                            "SA zero boundary trace did not propagate to effective coefficient");
        }
        const auto cell = [&](const char* name, Index i) {
            return detail::fieldData(static_cast<const ScalarField&>(problem.scalarField(name)))[i];
        };
        const auto check = [](double a, double b) {
            require(std::isfinite(a) && std::isfinite(b) &&
                std::abs(a-b) <= 1e-10*std::max({1.0,std::abs(a),std::abs(b)}),
                "model term differs from independently evaluated published equation");
        };
        for (Index i : detail::meshData(mesh).owned_cells) {
            if (model_name == "kEpsilon" || model_name == "kOmega") {
                const double k = cell("k", i);
                const double second = cell(model_name == "kEpsilon" ? "epsilon" : "omega", i);
                const double mut = model_name == "kEpsilon" ? 0.09*rho*k*k/second : rho*k/second;
                const double production = 4*mut;
                check(cell("mut", i), mut);
                check(cell("ransProduction", i), production);
                check(cell("ransSourceK", i) - cell("ransSinkK", i)*k,
                      production - (model_name == "kEpsilon" ? rho*second : 0.09*rho*second*k));
                if (model_name == "kEpsilon") {
                    check(cell("ransSourceEpsilon", i) - cell("ransSinkSecond", i)*second,
                          1.44*production*second/k - 1.92*rho*second*second/k);
                    check(cell("ransDiffusivityK", i), mu+mut);
                    check(cell("ransDiffusivityEpsilon", i), mu+mut/1.3);
                } else {
                    check(cell("ransSourceOmega", i) - cell("ransSinkSecond", i)*second,
                          (5.0/9.0)*production*second/k - 0.075*rho*second*second);
                    check(cell("ransDiffusivityK", i), mu+0.5*mut);
                    check(cell("ransDiffusivityOmega", i), mu+0.5*mut);
                }
            } else {
                const double nu = cell("nuTilda", i), d = cell("wallDistance", i);
                const double chi = nu/(mu/rho), kappa = 0.41, sigma = 2.0/3.0;
                const double fv1 = std::pow(chi,3)/(std::pow(chi,3)+std::pow(7.1,3));
                const double fv2 = 1-chi/(1+chi*fv1), ft2 = 1.2*std::exp(-0.5*chi*chi);
                const double st = std::max(2 + nu*fv2/(kappa*kappa*d*d), 1e-30);
                const double r = std::clamp(nu/(st*kappa*kappa*d*d),0.0,10.0);
                const double g = r+0.3*(std::pow(r,6)-r);
                const double fw = g*std::pow(65.0/(std::pow(g,6)+64),1.0/6.0);
                const double cw1 = 0.1355/(kappa*kappa)+(1+0.622)/sigma;
                const double source = rho*(0.1355*(1-ft2)*st*nu -
                    (cw1*fw-0.1355*ft2/(kappa*kappa))*nu*nu/(d*d) +
                    (0.622/sigma)*cell("ransGradNuTilda2",i));
                check(cell("mut",i), rho*nu*fv1);
                check(cell("ransDiffusivityNuTilda",i), (mu+rho*nu)/sigma);
                check(cell("ransSourceNuTilda",i)-cell("ransSinkNuTilda",i)*nu,source);
            }
        }

        // F1: div(mu*grad(U)^T) = (0,2,0), although div(mu*grad(U)) = 0.
        ScalarField coefficient(mesh, FieldLocation::Cell, "coefficient");
        coefficient.useCalculatedBoundary(); coefficient.evaluate([](Vec3 x) { return 1+x.x; });
        TensorField gradient(mesh, FieldLocation::Cell), stress(mesh, FieldLocation::Cell);
        VectorField correction(mesh, FieldLocation::Cell);
        evaluateStressCorrection(U, phi, coefficient, gradient, stress, correction);
        for (Index i : detail::meshData(mesh).owned_cells)
            require(norm(detail::fieldData(correction)[i]-Vec3{0,2,0}) < 1e-10,
                    "RANS transposed stress term is missing or incorrectly indexed");
        diagnostics::report("rans_equations_test: published source/diffusion/eddy viscosity and variable-coefficient stress passed");
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
