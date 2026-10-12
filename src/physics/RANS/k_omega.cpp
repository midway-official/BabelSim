#include "api.h"
#include "babelsim/equ.h"
#include "babelsim/builtin_descriptions.h"
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
             ScalarField& effectiveViscosity, double density, double viscosity,
             const OperatorOptions& hostMomentumOptions)
        : KOmega(problem, velocity, flux, effectiveViscosity, density, viscosity, hostMomentumOptions,
            problem.read.config(builtin::kOmega())) {}

    KOmega(Case& problem, VectorField& velocity, const ScalarField& flux,
             ScalarField& effectiveViscosity, double density, double viscosity,
             const OperatorOptions& hostMomentumOptions, const builtin::KOmegaConfig& parameters)
        : U(velocity), phi(flux), muEff(effectiveViscosity), rho(density), mu(viscosity),
          k(problem.load.scalarField("k")), omega(problem.load.scalarField("omega")),
          mut(problem.create.scalarField("mut", 0.0)),
          relaxation(parameters.controls.relaxation),
          residualTolerance(parameters.controls.tolerance),
          kMin(parameters.kMin),
          omegaMin(parameters.omegaMin),
          viscosityRatioMax(parameters.viscosityRatioMax),
          betaStar(parameters.betaStar),
          beta(parameters.beta),
          gamma(parameters.gamma),
          sigmaK(parameters.sigmaK),
          sigmaOmega(parameters.sigmaOmega),
          momentumOptions(hostMomentumOptions),
          kControl(readEquationControl(problem, builtin::kOmega().equation("kTransport"), k)),
          omegaControl(readEquationControl(problem, builtin::kOmega().equation("omegaTransport"), omega)),
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
    double relaxationFactor() const override { return relaxation; }

    void registerTimeState(time::HistoryGroup& histories) override {
        histories.add(k, kHistory);
        histories.add(omega, omegaHistory);
        histories.trackForRollback(mut);
        histories.watch(correction);
        histories.watch(boundedOmegaCells);
    }
    void beginStep() override { correction = 0; }

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
                   ScalarField& effectiveViscosity, double density, double viscosity,
             const OperatorOptions& hostMomentumOptions) {
    return new KOmega(problem, velocity, flux, effectiveViscosity, density, viscosity, hostMomentumOptions);
}

} // namespace babelsim::rans
