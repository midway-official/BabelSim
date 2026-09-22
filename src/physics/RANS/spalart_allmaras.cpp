#include "api.h"
#include "babelsim/equ.h"

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
                   ScalarField& effectiveViscosity, double density, double viscosity)
        : U(velocity), phi(flux), muEff(effectiveViscosity), rho(density), mu(viscosity),
          nuTilda(problem.scalarField("nuTilda")),
          wallDistance(problem.scalarField("wallDistance")),
          mut(problem.createScalarField("mut", 0.0)),
          relaxation(problem.solution().fraction("turbulenceRelaxation", 0.7)),
          residualTolerance(problem.solution().positive("turbulenceTolerance", 1e-6)),
          nuTildaMin(problem.physics().positive("saNuTildaMin", 1e-14)),
          wallDistanceMin(problem.physics().positive("saWallDistanceMin", 1e-12)),
          cb1(problem.physics().positive("saCb1", 0.1355)),
          cb2(problem.physics().positive("saCb2", 0.622)),
          sigma(problem.physics().positive("saSigma", 2.0 / 3.0)),
          kappa(problem.physics().positive("saKappa", 0.41)),
          cw2(problem.physics().positive("saCw2", 0.3)),
          cw3(problem.physics().positive("saCw3", 2.0)),
          cv1(problem.physics().positive("saCv1", 7.1)),
          ct3(problem.physics().positive("saCt3", 1.2)),
          ct4(problem.physics().positive("saCt4", 0.5)),
          cw1(problem.physics().positive("saCw1",
              cb1 / (kappa * kappa) + (1.0 + cb2) / sigma)),
          momentumOptions(readEquationControl(problem, "momentum", velocity,
              {"convection", "diffusion"}).spatial),
          nuTildaControl(readEquationControl(problem, "nuTildaTransport", nuTilda, {"convection", "diffusion"})),
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

    void saveOld(double dt) override {
        nuTildaHistory.save(nuTilda, dt);
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
                          ScalarField& effectiveViscosity, double density, double viscosity) {
    return new SpalartAllmaras(problem, velocity, flux, effectiveViscosity, density, viscosity);
}

} // namespace babelsim::rans
