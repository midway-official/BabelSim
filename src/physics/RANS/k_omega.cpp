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

// Wilcox1988m k-omega; constant-density incompressible transport.
class KOmega final : public Model {
public:
    KOmega(Case& problem, VectorField& velocity, const ScalarField& flux,
             ScalarField& effectiveViscosity, double density, double viscosity)
        : U(velocity), phi(flux), muEff(effectiveViscosity), rho(density), mu(viscosity),
          k(problem.scalarField("k")), omega(problem.scalarField("omega")),
          mut(problem.createScalarField("mut", 0.0)),
          relaxation(problem.solution().fraction("turbulenceRelaxation", 0.7)),
          residualTolerance(problem.solution().positive("turbulenceTolerance", 1e-6)),
          kMin(problem.physics().positive("kMin", 1e-12)),
          omegaMin(problem.physics().positive("omegaMin", 1e-12)),
          betaStar(problem.physics().positive("kOmegaBetaStar", 0.09)),
          beta(problem.physics().positive("kOmegaBeta", 3.0 / 40.0)),
          gamma(problem.physics().positive("kOmegaGamma", 5.0 / 9.0)),
          sigmaK(problem.physics().positive("kOmegaSigmaK", 0.5)),
          sigmaOmega(problem.physics().positive("kOmegaSigmaOmega", 0.5)),
          momentumOptions(readEquationControl(problem, "momentum", velocity,
              {"convection", "diffusion"}).spatial),
          kControl(readEquationControl(problem, "kTransport", k, {"convection", "diffusion"})),
          omegaControl(readEquationControl(problem, "omegaTransport", omega, {"convection", "diffusion"})),
          kHistory(time::history(k)), omegaHistory(time::history(omega))
    {
        mut.useCalculatedBoundary();
        muEff.useCalculatedBoundary();
        boundTransportFields();
        updateViscosity();
        problem.output(mut);
    }

    const char* modelName() const override { return "Wilcox1988m k-omega"; }
    double tolerance() const override { return residualTolerance; }

    void saveOld(double dt) override {
        kHistory.save(k, dt);
        omegaHistory.save(omega, dt);
    }

    TransportResult solveTransport() override {
        const ScalarField previousK = k;
        const ScalarField previousOmega = omega;
        U.setBoundaryFlux(phi);

        // Freeze closure coefficients at the current nonlinear iterate.
        updateViscosity();
        const auto production = mut * math::map(math::grad(U, momentumOptions), strainSquared);
        const auto kDiffusivity = mu + sigmaK * mut;
        const auto omegaDiffusivity = mu + sigmaOmega * mut;
        const auto kDestructionRate = betaStar * rho * omega;
        const auto omegaDestructionRate = beta * rho * omega;
        const auto omegaProduction = gamma * production * omega / math::max(k, kMin);

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

        // omega: use the same frozen closure state as the k equation.
        auto omegaEquation = equ::createEquation(omega, omegaControl);
        equ::ddt(omegaEquation, rho, omegaHistory);
        equ::div(omegaEquation, phi, rho, "convection");
        equ::laplacian(omegaEquation, omegaDiffusivity, -1, "diffusion");
        equ::reaction(omegaEquation, omegaDestructionRate);
        equ::source(omegaEquation, omegaProduction);
        const double omegaResidual = diagnostics::relativeResidual(omegaEquation, omega);
        equ::relax(omegaEquation, previousOmega, relaxation);
        const auto omegaSolve = equ::solve(omegaEquation);
        if (!diagnostics::all(omegaSolve.healthy()))
            return {{{"k", kSolve, kResidual, 0.0}, {"omega", omegaSolve, omegaResidual, 0.0}}};

        // Bound unknowns, then publish the viscosity for the next momentum solve.
        boundTransportFields();
        updateViscosity();
        return {{
            {"k", kSolve, kResidual, diagnostics::relativeChange(k, previousK)},
            {"omega", omegaSolve, omegaResidual,
                diagnostics::relativeChange(omega, previousOmega)}
        }};
    }

private:
    void boundTransportFields() {
        k.setBoundaryFlux(phi);
        omega.setBoundaryFlux(phi);
        k = math::max(k, kMin);
        omega = math::max(omega, omegaMin);
    }

    void updateViscosity() {
        mut = rho * k / math::max(omega, omegaMin);
        muEff = mu + mut;
    }

    VectorField& U;
    const ScalarField& phi;
    ScalarField& muEff;
    const double rho, mu;
    ScalarField& k;
    ScalarField& omega;
    ScalarField& mut;
    const double relaxation, residualTolerance, kMin, omegaMin;
    const double betaStar, beta, gamma, sigmaK, sigmaOmega;
    const OperatorOptions momentumOptions;
    const EquationControl kControl, omegaControl;
    time::History<double> kHistory, omegaHistory;
};

} // namespace

Model* makeKOmega(Case& problem, VectorField& velocity, const ScalarField& flux,
                   ScalarField& effectiveViscosity, double density, double viscosity) {
    return new KOmega(problem, velocity, flux, effectiveViscosity, density, viscosity);
}

} // namespace babelsim::rans
