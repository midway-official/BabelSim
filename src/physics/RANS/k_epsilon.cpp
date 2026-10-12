#include "api.h"
#include "babelsim/equ.h"
#include "babelsim/builtin_descriptions.h"

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
    KEpsilon(Case& problem, VectorField& velocity, const ScalarField& flux,
             ScalarField& effectiveViscosity, double density, double viscosity)
        : KEpsilon(problem, velocity, flux, effectiveViscosity, density, viscosity,
            problem.parameters(builtin::kEpsilon())) {}

    KEpsilon(Case& problem, VectorField& velocity, const ScalarField& flux,
             ScalarField& effectiveViscosity, double density, double viscosity,
             const ModuleParameters& parameters)
        : U(velocity), phi(flux), muEff(effectiveViscosity), rho(density), mu(viscosity),
          k(problem.scalarField("k")), epsilon(problem.scalarField("epsilon")),
          mut(problem.createScalarField("mut", 0.0)),
          relaxation(parameters.solution.number("turbulenceRelaxation")),
          residualTolerance(parameters.solution.number("turbulenceTolerance")),
          kMin(parameters.physics.number("kMin")),
          epsilonMin(parameters.physics.number("epsilonMin")),
          Cmu(parameters.physics.number("kEpsilonCmu")),
          C1(parameters.physics.number("kEpsilonC1")),
          C2(parameters.physics.number("kEpsilonC2")),
          sigmaK(parameters.physics.number("kEpsilonSigmaK")),
          sigmaEpsilon(parameters.physics.number("kEpsilonSigmaEpsilon")),
          momentumOptions(readEquationControl(problem, "momentum", velocity,
              {"convection", "diffusion"}).spatial),
          kControl(readEquationControl(problem, "kTransport", k, {"convection", "diffusion"})),
          epsilonControl(readEquationControl(problem, "epsilonTransport", epsilon, {"convection", "diffusion"})),
          kHistory(time::history(k)), epsilonHistory(time::history(epsilon))
    {
        mut.useCalculatedBoundary();
        muEff.useCalculatedBoundary();
        boundTransportFields();
        updateViscosity();
        problem.output(mut);
    }

    const char* modelName() const override { return "standard k-epsilon"; }
    double tolerance() const override { return residualTolerance; }

    void addHistories(time::HistoryGroup& histories) override {
        histories.add(k, kHistory);
        histories.add(epsilon, epsilonHistory);
        histories.track(mut);
    }

    TransportResult solveTransport() override {
        const ScalarField previousK = k;
        const ScalarField previousEpsilon = epsilon;
        U.setBoundaryFlux(phi);

        // Freeze closure coefficients at the current nonlinear iterate.
        updateViscosity();
        const auto production = mut * math::map(math::grad(U, momentumOptions), strainSquared);
        const auto kDiffusivity = mu + mut / sigmaK;
        const auto epsilonDiffusivity = mu + mut / sigmaEpsilon;
        const auto epsilonOverK = epsilon / math::max(k, kMin);
        const auto kDestructionRate = rho * epsilonOverK;
        const auto epsilonDestructionRate = C2 * rho * epsilonOverK;
        const auto epsilonProduction = C1 * production * epsilonOverK;

        // k: production on RHS; destruction linearized as rate * k on LHS.
        auto kEquation = equ::createEquation(k, kControl);
        equ::ddt(kEquation, rho, kHistory);
        equ::div(kEquation, phi, rho, "convection");
        equ::laplacian(kEquation, kDiffusivity, -1, "diffusion");
        equ::reaction(kEquation, kDestructionRate);
        equ::source(kEquation, production);
        const double kResidual = diagnostics::relativeResidual(kEquation, k);
        equ::relax(kEquation, previousK, relaxation);
        const auto kSolve = equ::solve(kEquation);
        if (!diagnostics::all(kSolve.healthy()))
            return {{{"k", kSolve, kResidual, 0.0}}};

        // epsilon: use the same frozen closure state as the k equation.
        auto epsilonEquation = equ::createEquation(epsilon, epsilonControl);
        equ::ddt(epsilonEquation, rho, epsilonHistory);
        equ::div(epsilonEquation, phi, rho, "convection");
        equ::laplacian(epsilonEquation, epsilonDiffusivity, -1, "diffusion");
        equ::reaction(epsilonEquation, epsilonDestructionRate);
        equ::source(epsilonEquation, epsilonProduction);
        const double epsilonResidual = diagnostics::relativeResidual(epsilonEquation, epsilon);
        equ::relax(epsilonEquation, previousEpsilon, relaxation);
        const auto epsilonSolve = equ::solve(epsilonEquation);
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

    VectorField& U;
    const ScalarField& phi;
    ScalarField& muEff;
    const double rho, mu;
    ScalarField& k;
    ScalarField& epsilon;
    ScalarField& mut;
    const double relaxation, residualTolerance, kMin, epsilonMin;
    const double Cmu, C1, C2, sigmaK, sigmaEpsilon;
    const OperatorOptions momentumOptions;
    const EquationControl kControl, epsilonControl;
    time::History<double> kHistory, epsilonHistory;
};

} // namespace

Model* makeKEpsilon(Case& problem, VectorField& velocity, const ScalarField& flux,
                   ScalarField& effectiveViscosity, double density, double viscosity) {
    return new KEpsilon(problem, velocity, flux, effectiveViscosity, density, viscosity);
}

} // namespace babelsim::rans
