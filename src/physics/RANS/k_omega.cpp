#include "api.h"
#include "babelsim/equ.h"
#include "babelsim/monitor.h"

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
          viscosityRatioMax(problem.physics().positive("maxTurbulentViscosityRatio", 1e5)),
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
        correction = 0;
        kHistory.save(k, dt);
        omegaHistory.save(omega, dt);
    }

    TransportResult solveTransport() override {
        boundedOmegaCells = 0.0;
        const ScalarField previousK = k;
        const ScalarField previousOmega = omega;
        U.setBoundaryFlux(phi);

        // Freeze closure coefficients at the current nonlinear iterate.
        updateViscosity();
        const auto production = mut * math::map(math::grad(U, momentumOptions), strainSquared);
        const auto kDiffusivity = mu + sigmaK * mut;
        const auto omegaDiffusivity = mu + sigmaOmega * mut;
        const auto omegaDestructionRate = beta * rho * omega;
        const auto omegaProduction = gamma * production * omega / math::max(k, kMin);

        // Solve omega first, then use its new value in k's implicit sink.
        // Production and diffusivity stay frozen for this segregated update.
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
            return {{{"omega", omegaSolve, omegaResidual, 0.0}}};
        const double rawOmegaMin = -math::max(-omega);
        boundOmega();

        auto kEquation = equ::createEquation(k, kControl);
        equ::ddt(kEquation, rho, kHistory);
        equ::div(kEquation, phi, rho, "convection");
        equ::laplacian(kEquation, kDiffusivity, -1, "diffusion");
        equ::reaction(kEquation, betaStar * rho * omega);
        equ::source(kEquation, production);
        const double kResidual = diagnostics::relativeResidual(kEquation, k);
        equ::relax(kEquation, previousK, relaxation);
        const auto kSolve = equ::solve(kEquation);
        if (!diagnostics::all(kSolve.healthy()))
            return {{{"omega", omegaSolve, omegaResidual, 0.0}, {"k", kSolve, kResidual, 0.0}}};

        const double rawKMin = -math::max(-k);
        // Bound unknowns, then publish the viscosity for the next momentum solve.
        boundTransportFields();
        updateViscosity();
        ++correction;
        if (correction == 1 || correction % 100 == 0 || boundedOmegaCells > 0.0) {
            monitor::Reporter("kOmega").record({{"iteration", correction},
                {"rK", kResidual}, {"rOmega", omegaResidual},
                {"kMinRaw", rawKMin}, {"omegaMinRaw", rawOmegaMin},
                {"kMax", math::max(k)}, {"omegaMax", math::max(omega)},
                {"mutMax", math::max(mut)}, {"productionMax", math::max(production)},
                {"omegaBoundedCells", boundedOmegaCells}});
        }
        return {{
            {"omega", omegaSolve, omegaResidual,
                diagnostics::relativeChange(omega, previousOmega)},
            {"k", kSolve, kResidual, diagnostics::relativeChange(k, previousK)}
        }};
    }

private:
    void boundTransportFields() {
        k.setBoundaryFlux(phi);
        omega.setBoundaryFlux(phi);
        k = math::max(k, kMin);
        boundOmega();
    }

    void boundOmega() {
        // A constant tiny omega floor alone permits k/omega to become enormous
        // after a negative BDF2 overshoot. Bound the viscosity ratio as well.
        boundedOmegaCells = std::max(boundedOmegaCells,
            math::sum(fieldBinary(omega, k, [&](double w, double energy) {
                return w < std::max(omegaMin, rho * energy / (viscosityRatioMax * mu)) ? 1.0 : 0.0;
            })));
        omega = fieldBinary(omega, k, [&](double w, double energy) {
            return std::max({w, omegaMin, rho * energy / (viscosityRatioMax * mu)});
        });
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
    const double relaxation, residualTolerance, kMin, omegaMin, viscosityRatioMax;
    const double betaStar, beta, gamma, sigmaK, sigmaOmega;
    const OperatorOptions momentumOptions;
    const EquationControl kControl, omegaControl;
    time::History<double> kHistory, omegaHistory;
    int correction = 0;
    double boundedOmegaCells = 0.0;
};

} // namespace

Model* makeKOmega(Case& problem, VectorField& velocity, const ScalarField& flux,
                   ScalarField& effectiveViscosity, double density, double viscosity) {
    return new KOmega(problem, velocity, flux, effectiveViscosity, density, viscosity);
}

} // namespace babelsim::rans
