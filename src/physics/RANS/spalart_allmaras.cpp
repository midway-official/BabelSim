#include "api.h"
#include "babelsim/equ.h"
#include "babelsim/builtin_descriptions.h"

#include <algorithm>
#include <cmath>

namespace babelsim::rans {
namespace {

double vorticityMagnitude(const Tensor3& gradient) {
    double sum = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const double rotation = 0.5 * (gradient[i][j] - gradient[j][i]);
            sum += rotation * rotation;
        }
    }
    return std::sqrt(2.0 * sum);
}

// Standard positive-variable SA, including ft2, without the primary trip term.
// wallDistance is the geometric distance to the nearest wall, supplied by Case.
class SpalartAllmaras final : public Model {
public:
    SpalartAllmaras(Case& problem, VectorField& velocity, const ScalarField& flux,
                   ScalarField& effectiveViscosity, double density, double viscosity,
             const OperatorOptions& hostMomentumOptions)
        : SpalartAllmaras(problem, velocity, flux, effectiveViscosity, density, viscosity, hostMomentumOptions,
            problem.read.config(builtin::spalartAllmaras())) {}

    SpalartAllmaras(Case& problem, VectorField& velocity, const ScalarField& flux,
                   ScalarField& effectiveViscosity, double density, double viscosity,
             const OperatorOptions& hostMomentumOptions, const builtin::SAConfig& parameters)
        : U(velocity), phi(flux), muEff(effectiveViscosity), rho(density), mu(viscosity),
          nuTilda(problem.load.scalarField("nuTilda")),
          wallDistance(problem.load.scalarField("wallDistance")),
          mut(problem.create.scalarField("mut", 0.0)),
          relaxation(parameters.controls.relaxation),
          residualTolerance(parameters.controls.tolerance),
          nuTildaMin(parameters.nuTildaMin),
          wallDistanceMin(parameters.wallDistanceMin),
          cb1(parameters.cb1),
          cb2(parameters.cb2),
          sigma(parameters.sigma),
          kappa(parameters.kappa),
          cw2(parameters.cw2),
          cw3(parameters.cw3),
          cv1(parameters.cv1),
          ct3(parameters.ct3),
          ct4(parameters.ct4),
          cw1(parameters.cw1),
          momentumOptions(hostMomentumOptions),
          nuTildaControl(readEquationControl(problem, builtin::spalartAllmaras().equation("nuTildaTransport"), nuTilda)),
          nuTildaHistory(time::history(nuTilda))
    {
        mut.useCalculatedBoundary();
        muEff.useCalculatedBoundary();
        boundTransportField();
        updateViscosity();
        problem.output(mut);
    }

    const char* modelName() const override { return "Spalart-Allmaras"; }
    double tolerance() const override { return residualTolerance; }
    double relaxationFactor() const override { return relaxation; }

    void registerTimeState(time::HistoryGroup& histories) override {
        histories.add(nuTilda, nuTildaHistory);
        histories.trackForRollback(mut);
    }

    TransportResult solveTransport() override {
        const ScalarField previousNuTilda = nuTilda;
        U.setBoundaryFlux(phi);

        // Closure functions, evaluated at the current nonlinear iterate.
        const double nu = mu / rho;
        const auto chi = nuTilda / nu;
        const auto fv1 = viscosityDamping(chi);
        const auto fv2 = 1.0 - chi / (1.0 + chi * fv1);
        const auto ft2 = math::map(chi, [this](double x) {
            return ct3 * std::exp(-ct4 * x * x);
        });
        const auto distance = math::max(wallDistance, wallDistanceMin);
        const auto inverseDistanceSquared = 1.0 / (distance * distance);
        const auto rotation = math::map(math::grad(U, momentumOptions), vorticityMagnitude);
        // Existing positive S-tilde regularization; this is not SA-neg.
        const auto modifiedRotation = math::max(
            rotation + nuTilda * fv2 * inverseDistanceSquared / (kappa * kappa), 1e-30);
        const auto r = math::map(
            nuTilda * inverseDistanceSquared / (kappa * kappa * modifiedRotation),
            [](double value) { return std::clamp(value, 0.0, 10.0); });
        const auto fw = math::map(r, [this](double value) {
            const double g = value + cw2 * (std::pow(value, 6.0) - value);
            const double cw3ToSix = std::pow(cw3, 6.0);
            return g * std::pow((1.0 + cw3ToSix) /
                (std::pow(g, 6.0) + cw3ToSix), 1.0 / 6.0);
        });

        // Net reaction a*nuTilda: positive part explicit, negative part implicit.
        const auto productionRate = cb1 * (1.0 - ft2) * modifiedRotation;
        const auto destructionRate = (cw1 * fw - cb1 / (kappa * kappa) * ft2)
            * nuTilda * inverseDistanceSquared;
        const auto netReactionRate = productionRate - destructionRate;
        const auto implicitDestruction = rho * math::max(-netReactionRate, 0.0);
        const auto gradNuTilda = math::grad(nuTilda, nuTildaControl.spatial);
        const auto gradientSource = (rho * cb2 / sigma) * math::dot(gradNuTilda, gradNuTilda);
        const auto explicitSource = rho * math::max(netReactionRate, 0.0) * nuTilda
            + gradientSource;
        const auto diffusivity = (mu + rho * nuTilda) / sigma;

        // rho D(nuTilda)/Dt - div(diffusivity grad(nuTilda)) + sink = source.
        auto nuTildaEquation = equ::createEquation(nuTilda, nuTildaControl);
        equ::ddt(nuTildaEquation, rho, nuTildaHistory);
        equ::div(nuTildaEquation, phi, rho, "convection");
        equ::laplacian(nuTildaEquation, diffusivity, -1, "diffusion");
        equ::reaction(nuTildaEquation, implicitDestruction);
        equ::source(nuTildaEquation, explicitSource);
        const double transportResidual = diagnostics::relativeResidual(nuTildaEquation, nuTilda);
        equ::relax(nuTildaEquation, previousNuTilda, relaxation);
        const auto nuTildaSolve = equ::solve(nuTildaEquation);
        if (!diagnostics::all(nuTildaSolve.healthy()))
            return {{{"nuTilda", nuTildaSolve, transportResidual, 0.0}}};

        boundTransportField();
        updateViscosity();
        return {{{"nuTilda", nuTildaSolve, transportResidual,
            diagnostics::relativeChange(nuTilda, previousNuTilda)}}};
    }

private:
    ScalarField viscosityDamping(const ScalarField& chi) const {
        return math::map(chi, [this](double value) {
            const double chiCubed = value * value * value;
            return chiCubed / (chiCubed + cv1 * cv1 * cv1);
        });
    }

    void boundTransportField() {
        nuTilda.setBoundaryFlux(phi);
        nuTilda = math::max(nuTilda, nuTildaMin);
    }

    void updateViscosity() {
        mut = rho * nuTilda * viscosityDamping(nuTilda / (mu / rho));
        muEff = mu + mut;
    }

    VectorField& U;
    const ScalarField& phi;
    ScalarField& muEff;
    const double rho, mu;
    ScalarField& nuTilda;
    const ScalarField& wallDistance;
    ScalarField& mut;
    const double relaxation, residualTolerance, nuTildaMin, wallDistanceMin;
    const double cb1, cb2, sigma, kappa, cw2, cw3, cv1, ct3, ct4, cw1;
    const OperatorOptions momentumOptions;
    const EquationControl nuTildaControl;
    time::History<double> nuTildaHistory;
};

} // namespace

Model* makeSpalartAllmaras(Case& problem, VectorField& velocity, const ScalarField& flux,
                          ScalarField& effectiveViscosity, double density, double viscosity,
             const OperatorOptions& hostMomentumOptions) {
    return new SpalartAllmaras(problem, velocity, flux, effectiveViscosity, density, viscosity, hostMomentumOptions);
}

} // namespace babelsim::rans
