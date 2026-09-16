#include "api.h"
#include "babelsim/equ.h"

#include <algorithm>
#include <cmath>

namespace babelsim::rans {
namespace {

// 2 dev(symm(grad U)):dev(symm(grad U)); kept local to this model.
double strainSquared(const Tensor3& gradient) {
    const double divergence = trace(gradient);
    double sum = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const double strain = 0.5 * (gradient[i][j] + gradient[j][i])
                - (i == j ? divergence / 3.0 : 0.0);
            sum += strain * strain;
        }
    }
    return 2.0 * sum;
}

// standard k-epsilon; constant-density incompressible transport.
class KEpsilon final : public Model {
public:
    KEpsilon(Case& problem, const VectorField& velocity, const ScalarField& flux,
             ScalarField& effectiveViscosity, double density, double viscosity)
        : U(velocity), phi(flux), muEff(effectiveViscosity), rho(density), mu(viscosity),
          k(problem.scalarField("k")), epsilon(problem.scalarField("epsilon")),
          mut(problem.scalarField("mut", 0.0)),
          relaxation(problem.physics().fraction("turbulenceRelaxation", 0.7)),
          residualTolerance(problem.physics().positive("turbulenceTolerance", 1e-6)),
          kMin(problem.physics().positive("kMin", 1e-12)),
          epsilonMin(problem.physics().positive("epsilonMin", 1e-12)),
          Cmu(problem.physics().positive("kEpsilonCmu", 0.09)),
          C1(problem.physics().positive("kEpsilonC1", 1.44)),
          C2(problem.physics().positive("kEpsilonC2", 1.92)),
          sigmaK(problem.physics().positive("kEpsilonSigmaK", 1.0)),
          sigmaEpsilon(problem.physics().positive("kEpsilonSigmaEpsilon", 1.3)),
          kSolver(readLinearControl(problem, k)),
          epsilonSolver(readLinearControl(problem, epsilon)),
          kHistory(math::history(k)), epsilonHistory(math::history(epsilon))
    {
        mut.useCalculatedBoundary();
        muEff.useCalculatedBoundary();
        boundTransportFields();
        updateViscosity();
        problem.output(mut);
    }

    const char* modelName() const override { return "standard k-epsilon"; }
    double tolerance() const override { return residualTolerance; }

    void saveOld(double dt) override {
        math::saveOld(kHistory, k, dt);
        math::saveOld(epsilonHistory, epsilon, dt);
    }

    TransportResult solveTransport() override {
        const auto previousK = math::copy(k);
        const auto previousEpsilon = math::copy(epsilon);
        const_cast<VectorField&>(U).setBoundaryFlux(phi);

        // Freeze closure coefficients at the current nonlinear iterate.
        updateViscosity();
        const auto production = mut * math::map(math::grad(U), strainSquared);
        const auto kDiffusivity = mu + mut / sigmaK;
        const auto epsilonDiffusivity = mu + mut / sigmaEpsilon;
        const auto epsilonOverK = epsilon / math::max(k, kMin);
        const auto kDestructionRate = rho * epsilonOverK;
        const auto epsilonDestructionRate = C2 * rho * epsilonOverK;
        const auto epsilonProduction = C1 * production * epsilonOverK;

        // k: production on RHS; destruction linearized as rate * k on LHS.
        auto kEquation = equ::createEquation(k);
        equ::ddt(kEquation, rho, kHistory);
        equ::div(kEquation, phi, rho);
        equ::laplacian(kEquation, kDiffusivity, -1.0);
        equ::reaction(kEquation, kDestructionRate);
        equ::source(kEquation, production);
        const double kResidual = equ::relativeResidual(kEquation, k);
        equ::relax(kEquation, previousK, relaxation);
        const auto kSolve = equ::solve(kEquation, k, kSolver);
        if (!diagnostics::all(kSolve.healthy()))
            return {{{"k", kSolve, kResidual, 0.0}}};

        // epsilon: use the same frozen closure state as the k equation.
        auto epsilonEquation = equ::createEquation(epsilon);
        equ::ddt(epsilonEquation, rho, epsilonHistory);
        equ::div(epsilonEquation, phi, rho);
        equ::laplacian(epsilonEquation, epsilonDiffusivity, -1.0);
        equ::reaction(epsilonEquation, epsilonDestructionRate);
        equ::source(epsilonEquation, epsilonProduction);
        const double epsilonResidual = equ::relativeResidual(epsilonEquation, epsilon);
        equ::relax(epsilonEquation, previousEpsilon, relaxation);
        const auto epsilonSolve = equ::solve(epsilonEquation, epsilon, epsilonSolver);
        if (!diagnostics::all(epsilonSolve.healthy()))
            return {{{"k", kSolve, kResidual, 0.0}, {"epsilon", epsilonSolve, epsilonResidual, 0.0}}};

        // Bound unknowns, then publish the viscosity for the next momentum solve.
        boundTransportFields();
        updateViscosity();
        return {{
            {"k", kSolve, kResidual, diagnostics::relativeChange(k, previousK)},
            {"epsilon", epsilonSolve, epsilonResidual,
                diagnostics::relativeChange(epsilon, previousEpsilon)}
        }};
    }

private:
    void boundTransportFields() {
        k.setBoundaryFlux(phi);
        epsilon.setBoundaryFlux(phi);
        k = math::max(k, kMin);
        epsilon = math::max(epsilon, epsilonMin);
    }

    void updateViscosity() {
        mut = Cmu * rho * k * k / math::max(epsilon, epsilonMin);
        muEff = mu + mut;
    }

    const VectorField& U;
    const ScalarField& phi;
    ScalarField& muEff;
    const double rho, mu;
    ScalarField& k;
    ScalarField& epsilon;
    ScalarField& mut;
    const double relaxation, residualTolerance, kMin, epsilonMin;
    const double Cmu, C1, C2, sigmaK, sigmaEpsilon;
    const LinearSolverConfig kSolver, epsilonSolver;
    math::History<double> kHistory, epsilonHistory;
};

} // namespace

Model* makeKEpsilon(Case& problem, const VectorField& velocity, const ScalarField& flux,
                   ScalarField& effectiveViscosity, double density, double viscosity) {
    return new KEpsilon(problem, velocity, flux, effectiveViscosity, density, viscosity);
}

} // namespace babelsim::rans
