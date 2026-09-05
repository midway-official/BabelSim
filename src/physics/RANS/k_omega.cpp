#include "model.h"

#include <algorithm>

namespace babelsim::rans {
namespace {

// Wilcox 1988m 高雷诺数 k-omega。该明确命名的版本不含 2006 版交叉扩散，
// P 使用不可压缩线性涡黏性近似 mu_t * 2*S':S'。
class KOmega final : public Model {
public:
    KOmega(
        Case& problem,
        const VectorField& velocity,
        const ScalarField& face_flux,
        ScalarField& effective_viscosity,
        double density,
        double molecular_viscosity)
        : Model(problem, velocity, face_flux, effective_viscosity,
                density, molecular_viscosity),
          m_beta_star(positiveSetting(problem.physics(), "kOmegaBetaStar", 0.09)),
          m_beta(positiveSetting(problem.physics(), "kOmegaBeta", 3.0 / 40.0)),
          m_gamma(positiveSetting(problem.physics(), "kOmegaGamma", 5.0 / 9.0)),
          m_sigma_k(positiveSetting(problem.physics(), "kOmegaSigmaK", 0.5)),
          m_sigma_omega(positiveSetting(problem.physics(), "kOmegaSigmaOmega", 0.5)),
          m_k_min(positiveSetting(problem.physics(), "kMin", 1e-12)),
          m_omega_min(positiveSetting(problem.physics(), "omegaMin", 1e-12)),
          m_k(problem.scalarField("k")),
          m_omega(problem.scalarField("omega")),
          m_mut(problem.scalarField("mut", 0.0)),
          m_previous_k(problem.scalarField("ransPreviousK", 0.0)),
          m_previous_omega(problem.scalarField("ransPreviousOmega", 0.0)),
          m_inverse_k(problem.scalarField("ransInverseK", 0.0)),
          m_inverse_omega(problem.scalarField("ransInverseOmega", 0.0)),
          m_production(problem.scalarField("ransProduction", 0.0)),
          m_diffusivity_k(problem.scalarField("ransDiffusivityK", 0.0)),
          m_diffusivity_omega(problem.scalarField("ransDiffusivityOmega", 0.0)),
          m_source_k(problem.scalarField("ransSourceK", 0.0)),
          m_source_omega(problem.scalarField("ransSourceOmega", 0.0)),
          m_work(problem.scalarField("ransWork", 0.0))
    {
        boundFields();
        updateKinematics();
        updateViscosity();
        m_problem.output(m_mut);
    }

    const char* modelName() const override { return "Wilcox1988m k-omega"; }

    SolveResult correct() override {
        m_previous_k.assign(m_k);
        m_previous_omega.assign(m_omega);
        updateKinematics();
        updateSourcesAndDiffusivities();

        const SolveResult k_result = solveTransport(m_k, m_diffusivity_k, m_source_k);
        const SolveResult omega_result = solveTransport(
            m_omega, m_diffusivity_omega, m_source_omega);
        boundFields();
        updateViscosity();
        m_relative_change = std::max(
            diagnostics::relativeChange(m_k, m_previous_k),
            diagnostics::relativeChange(m_omega, m_previous_omega));
        return combine(k_result, omega_result);
    }

    double relativeChange() const override { return m_relative_change; }

private:
    void boundFields() {
        m_k.evaluate(m_k, [this](double value) { return std::max(value, m_k_min); });
        m_omega.evaluate(m_omega, [this](double value) {
            return std::max(value, m_omega_min);
        });
    }

    void updateViscosity() {
        m_inverse_omega.evaluate(m_omega, [this](double value) {
            return 1.0 / std::max(value, m_omega_min);
        });
        m_mut.assignProduct(m_k, m_inverse_omega);
        m_mut.assignScaled(m_density, m_mut);
        setEddyViscosity(m_mut);
    }

    void updateSourcesAndDiffusivities() {
        updateViscosity();
        m_production.assignProduct(m_mut, m_strain_measure);

        m_work.assignProduct(m_k, m_omega);
        m_source_k.assign(m_production);
        m_source_k.addScaled(-m_beta_star * m_density, m_work);

        m_inverse_k.evaluate(m_k, [this](double value) {
            return 1.0 / std::max(value, m_k_min);
        });
        m_work.assignProduct(m_production, m_omega);
        m_work.assignProduct(m_inverse_k, m_work);
        m_source_omega.assignScaled(m_gamma, m_work);
        m_work.assignProduct(m_omega, m_omega);
        m_source_omega.addScaled(-m_beta * m_density, m_work);

        m_diffusivity_k.fill(m_molecular_viscosity);
        m_diffusivity_k.addScaled(m_sigma_k, m_mut);
        m_diffusivity_omega.fill(m_molecular_viscosity);
        m_diffusivity_omega.addScaled(m_sigma_omega, m_mut);
    }

    double m_beta_star;
    double m_beta;
    double m_gamma;
    double m_sigma_k;
    double m_sigma_omega;
    double m_k_min;
    double m_omega_min;
    ScalarField& m_k;
    ScalarField& m_omega;
    ScalarField& m_mut;
    ScalarField& m_previous_k;
    ScalarField& m_previous_omega;
    ScalarField& m_inverse_k;
    ScalarField& m_inverse_omega;
    ScalarField& m_production;
    ScalarField& m_diffusivity_k;
    ScalarField& m_diffusivity_omega;
    ScalarField& m_source_k;
    ScalarField& m_source_omega;
    ScalarField& m_work;
    double m_relative_change = 0.0;
};

}  // 匿名命名空间

Model* makeKOmega(
    Case& problem,
    const VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity)
{
    return new KOmega(
        problem, velocity, face_flux, effective_viscosity,
        density, molecular_viscosity);
}

}  // babelsim::rans 命名空间
