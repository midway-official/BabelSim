#include "babelsim/runtime.h"
#include "babelsim/equ.h"
#include "babelsim/mpi_support.h"
#include "babelsim/parallel.h"
#include "physics/RANS/api.h"
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
        auto& phi = problem.createFaceField("phi");
        phi = math::flux(U);
        const auto& mesh = problem.mesh();
        U.evaluate([](Vec3 x) { return Vec3{2*x.y,0,0}; });
        for (Index p = 0; p < mesh.patchCount(); ++p) U.setBoundary(p, BoundaryCondition<Vec3>::zeroGradient());
        U.setBoundary(2, fixedValue(Vec3{})); U.setBoundary(3, fixedValue(Vec3{2,0,0}));
        math::evaluate(math::flux(U), phi);
        auto& effective = problem.createScalarField("effective", 0.0);
        const double rho = problem.physics().positive("density");
        const double mu = problem.physics().positive("dynamicViscosity");
        const auto model_name = problem.physics().word("turbulenceModel");
        if (model_name == "SA")
            problem.scalarField("nuTilda").setBoundary(0, fixedValue(0.0));
        for(const char* key:{"maxIterations","nonOrthogonalCorrections","velocityRelaxation",
            "pressureRelaxation","continuityTolerance","velocityTolerance","momentumTolerance",
            "pressureCorrectionTolerance"})
            if(problem.solution().contains(key)) (void)problem.solution().number(key);
        std::unique_ptr<rans::Model, void(*)(rans::Model*)> model(
            rans::create(problem, U, phi, effective, rho, mu), rans::destroy);
        require(model != nullptr, "expected turbulence model");
        auto time=time::start(problem);
        time.advance();
        model->saveOld(time.dt());
        // Independently build one frozen-coefficient Euler step from the published
        // scalar formulas. Test observable solutions rather than private scratch fields.
        const double relaxation = problem.solution().fraction("turbulenceRelaxation", 0.7);
        std::vector<ScalarField> expected;
        std::vector<double> initialResiduals;
        const std::vector<std::string> names = model_name == "SA"
            ? std::vector<std::string>{"nuTilda"}
            : std::vector<std::string>{"k", model_name == "kOmega" ? "omega" : "epsilon"};
        for (const auto& name : names) {
            expected.push_back(problem.scalarField(name));
            auto& value = expected.back();
            const auto old = value;
            ScalarField diffusion(mesh, FieldLocation::Cell);
            ScalarField source(mesh, FieldLocation::Cell);
            ScalarField sink(mesh, FieldLocation::Cell);
            diffusion.useCalculatedBoundary();
            source.useCalculatedBoundary();
            sink.useCalculatedBoundary();
            if (model_name != "SA") {
                const auto& k = problem.scalarField("k");
                const auto& second = problem.scalarField(model_name == "kOmega" ? "omega" : "epsilon");
                diffusion.evaluate(k, second, [&](double kv, double sv) {
                    const double eddy = model_name == "kOmega" ? rho*kv/sv : 0.09*rho*kv*kv/sv;
                    return mu + eddy * (model_name == "kOmega" ? 0.5 : name == "k" ? 1.0 : 1.0/1.3);
                });
                source.evaluate(k, second, [&](double kv, double sv) {
                    const double eddy = model_name == "kOmega" ? rho*kv/sv : 0.09*rho*kv*kv/sv;
                    const double production = 4*eddy; // U=(2y,0,0): 2*S:S=4
                    return name == "k" ? production : production*sv/kv*(model_name == "kOmega" ? 5.0/9.0 : 1.44);
                });
                sink.evaluate(k, second, [&](double kv, double sv) {
                    return model_name == "kOmega" ? rho*sv*(name == "k" ? 0.09 : 0.075)
                        : rho*sv/kv*(name == "k" ? 1.0 : 1.92);
                });
            } else {
                const auto& distance = problem.scalarField("wallDistance");
                ScalarField reaction(mesh, FieldLocation::Cell);
                reaction.useCalculatedBoundary();
                reaction.evaluate(old, distance, [&](double nu, double d) {
                    const double chi = nu/(mu/rho), kappa = 0.41, sigma = 2.0/3.0;
                    const double fv1 = std::pow(chi,3)/(std::pow(chi,3)+std::pow(7.1,3));
                    const double fv2 = 1-chi/(1+chi*fv1), ft2 = 1.2*std::exp(-0.5*chi*chi);
                    const double st = std::max(2 + nu*fv2/(kappa*kappa*d*d), 1e-30);
                    const double r = std::clamp(nu/(st*kappa*kappa*d*d),0.0,10.0);
                    const double g = r+0.3*(std::pow(r,6)-r);
                    const double fw = g*std::pow(65.0/(std::pow(g,6)+64),1.0/6.0);
                    const double cw1 = 0.1355/(kappa*kappa)+(1+0.622)/sigma;
                    return 0.1355*(1-ft2)*st-(cw1*fw-0.1355*ft2/(kappa*kappa))*nu/(d*d);
                });
                diffusion.evaluate(old, [&](double nu) { return (mu+rho*nu)/(2.0/3.0); });
                sink.evaluate(reaction, [&](double rate) { return rho*std::max(-rate,0.0); });
                source.evaluate(reaction, old, [&](double rate, double nu) { return rho*std::max(rate,0.0)*nu; });
                auto gradient = math::grad(old);
                ScalarField gradientSource(mesh, FieldLocation::Cell);
                gradientSource.useCalculatedBoundary();
                gradientSource.evaluate(gradient, [&](Vec3 g) { return rho*0.622/(2.0/3.0)*squaredNorm(g); });
                source.addScaled(1.0, gradientSource);
            }
            auto equation = equ::createEquation(value);
            equ::ddt(equation, rho, old, time.dt());
            equ::div(equation, phi, rho);
            equ::laplacian(equation, diffusion, -1);
            equ::reaction(equation, sink);
            equ::source(equation, source);
            initialResiduals.push_back(diagnostics::relativeResidual(equation, value));
            equ::relax(equation, old, relaxation);
            require(equ::solve(equation, value, readLinearControl(problem,value)).healthy(),
                    "reference transport solve failed");
        }
        const auto report = model->solveTransport();
        require(report.healthy(), "model transport update failed");
        require(report.equations.size() == names.size(), "missing per-equation diagnostics");
        for (std::size_t n=0; n<names.size(); ++n) {
            require(report.equations[n].field == names[n], "incorrect diagnostic field name");
            require(near(report.equations[n].initialResidual, initialResiduals[n], 1e-10),
                    "transport residual was not measured before relaxation");
            const auto& actual = problem.scalarField(names[n]);
            for (Index i : detail::meshData(mesh).owned_cells)
                require(near(detail::fieldData(actual)[i], detail::fieldData(expected[n])[i], 1e-10),
                        "transport step differs from independent source/diffusion equation");
        }
        const auto cell = [&](const char* name, Index i) {
            return detail::fieldData(static_cast<const ScalarField&>(problem.existingScalarField(name)))[i];
        };
        for (Index i : detail::meshData(mesh).owned_cells) {
            double eddy = 0;
            if (model_name == "kOmega") eddy = rho*cell("k",i)/cell("omega",i);
            else if (model_name == "kEpsilon") eddy = 0.09*rho*std::pow(cell("k",i),2)/cell("epsilon",i);
            else {
                const double nu=cell("nuTilda",i), chi=nu/(mu/rho);
                eddy=rho*nu*std::pow(chi,3)/(std::pow(chi,3)+std::pow(7.1,3));
            }
            require(near(cell("mut",i),eddy,1e-10), "eddy viscosity closure differs from published formula");
            require(near(detail::fieldData(effective)[i],mu+eddy,1e-10), "effective viscosity was not updated");
        }
        if (model_name == "SA") {
            const auto face = math::interpolate(effective);
            for (Index f : detail::meshData(mesh).owned_faces)
                if (mesh.boundaryFace(f) && mesh.boundaryPatch(f) == 0)
                    require(std::abs(detail::fieldData(face)[f] - mu) < 1e-13,
                            "SA zero boundary trace did not propagate to effective coefficient");
        }

        // F1: div(mu*grad(U)^T) = (0,2,0), although div(mu*grad(U)) = 0.
        ScalarField coefficient(mesh, FieldLocation::Cell, "coefficient");
        coefficient.useCalculatedBoundary(); coefficient.evaluate([](Vec3 x) { return 1+x.x; });
        TensorField gradient(mesh, FieldLocation::Cell), stress(mesh, FieldLocation::Cell);
        VectorField correction(mesh, FieldLocation::Cell);
        U.setBoundaryFlux(phi);
        gradient.useCalculatedBoundary(); stress.useCalculatedBoundary();
        gradient=math::grad(U);
        stress=coefficient*(math::transpose(gradient)-(2.0/3.0)*math::isotropic(math::trace(gradient)));
        correction=math::div(stress);
        for (Index i : detail::meshData(mesh).owned_cells)
            require(norm(detail::fieldData(correction)[i]-Vec3{0,2,0}) < 1e-10,
                    "RANS transposed stress term is missing or incorrectly indexed");
        if(primaryProcess()) std::cout << "rans_equations_test: published source/diffusion/eddy viscosity and variable-coefficient stress passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
