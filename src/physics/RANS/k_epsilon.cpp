#include "model.h"

#include <algorithm>

namespace babelsim::rans {
namespace {

// 标准高雷诺数 k-epsilon（线性涡黏性、不可压缩形式）。破坏项采用上一轮
// 湍流场显式 Picard 线性化，避免把模型专用反应项加入通用 Equation API。
class KEpsilon final : public Model {
public:
    KEpsilon(
        Case& problem,
        const VectorField& velocity,
        const ScalarField& face_flux,
        ScalarField& effective_viscosity,
        double density,
        double molecular_viscosity)
        : Model(problem, velocity, face_flux, effective_viscosity,
                density, molecular_viscosity),
          m_c_mu(positiveSetting(problem.physics(), "kEpsilonCmu", 0.09)),
          m_c1(positiveSetting(problem.physics(), "kEpsilonC1", 1.44)),
          m_c2(positiveSetting(problem.physics(), "kEpsilonC2", 1.92)),
          m_sigma_k(positiveSetting(problem.physics(), "kEpsilonSigmaK", 1.0)),
          m_sigma_epsilon(positiveSetting(
              problem.physics(), "kEpsilonSigmaEpsilon", 1.3)),
          m_k_min(positiveSetting(problem.physics(), "kMin", 1e-12)),
          m_epsilon_min(positiveSetting(problem.physics(), "epsilonMin", 1e-12)),
          m_k(problem.scalarField("k")),
          m_epsilon(problem.scalarField("epsilon")),
          m_mut(problem.scalarField("mut", 0.0)),
          m_previous_k(problem.scalarField("ransPreviousK", 0.0)),
          m_previous_epsilon(problem.scalarField("ransPreviousEpsilon", 0.0)),
          m_inverse_k(problem.scalarField("ransInverseK", 0.0)),
          m_inverse_epsilon(problem.scalarField("ransInverseEpsilon", 0.0)),
          m_k_squared(problem.scalarField("ransKSquared", 0.0)),
          m_production(problem.scalarField("ransProduction", 0.0)),
          m_diffusivity_k(problem.scalarField("ransDiffusivityK", 0.0)),
          m_diffusivity_epsilon(problem.scalarField("ransDiffusivityEpsilon", 0.0)),
          m_source_k(problem.scalarField("ransSourceK", 0.0)),
          m_source_epsilon(problem.scalarField("ransSourceEpsilon", 0.0)),
          m_work(problem.scalarField("ransWork", 0.0))
    {
        boundFields();
        updateKinematics();
        updateViscosity();
        m_problem.output(m_mut);
    }

    const char* modelName() const override { return "standard k-epsilon"; }

    SolveResult correct() override {
        m_previous_k.assign(m_k);
        m_previous_epsilon.assign(m_epsilon);
        updateKinematics();
        updateSourcesAndDiffusivities();

        const SolveResult k_result = solveTransport(m_k, m_diffusivity_k, m_source_k);
        const SolveResult epsilon_result = solveTransport(
            m_epsilon, m_diffusivity_epsilon, m_source_epsilon);
        boundFields();
        updateViscosity();
        m_relative_change = std::max(
            diagnostics::relativeChange(m_k, m_previous_k),
            diagnostics::relativeChange(m_epsilon, m_previous_epsilon));
        return combine(k_result, epsilon_result);
    }

    double relativeChange() const override { return m_relative_change; }

private:
    void boundFields() {
        m_k.evaluate(m_k, [this](double value) { return std::max(value, m_k_min); });
        m_epsilon.evaluate(m_epsilon, [this](double value) {
            return std::max(value, m_epsilon_min);
        });
    }

    void updateViscosity() {
        m_inverse_epsilon.evaluate(m_epsilon, [this](double value) {
            return 1.0 / std::max(value, m_epsilon_min);
        });
        m_k_squared.assignProduct(m_k, m_k);
        m_mut.assignProduct(m_k_squared, m_inverse_epsilon);
        m_mut.assignScaled(m_c_mu * m_density, m_mut);
        setEddyViscosity(m_mut);
    }

    void updateSourcesAndDiffusivities() {
        updateViscosity();
        m_production.assignProduct(m_mut, m_strain_measure);

        m_source_k.assign(m_production);
        m_source_k.addScaled(-m_density, m_epsilon);

        m_inverse_k.evaluate(m_k, [this](double value) {
            return 1.0 / std::max(value, m_k_min);
        });
        m_work.assignProduct(m_production, m_epsilon);
        m_work.assignProduct(m_inverse_k, m_work);
        m_source_epsilon.assignScaled(m_c1, m_work);
        m_work.assignProduct(m_epsilon, m_epsilon);
        m_work.assignProduct(m_inverse_k, m_work);
        m_source_epsilon.addScaled(-m_c2 * m_density, m_work);

        m_diffusivity_k.fill(m_molecular_viscosity);
        m_diffusivity_k.addScaled(1.0 / m_sigma_k, m_mut);
        m_diffusivity_epsilon.fill(m_molecular_viscosity);
        m_diffusivity_epsilon.addScaled(1.0 / m_sigma_epsilon, m_mut);
    }

    double m_c_mu;
    double m_c1;
    double m_c2;
    double m_sigma_k;
    double m_sigma_epsilon;
    double m_k_min;
    double m_epsilon_min;
    ScalarField& m_k;
    ScalarField& m_epsilon;
    ScalarField& m_mut;
    ScalarField& m_previous_k;
    ScalarField& m_previous_epsilon;
    ScalarField& m_inverse_k;
    ScalarField& m_inverse_epsilon;
    ScalarField& m_k_squared;
    ScalarField& m_production;
    ScalarField& m_diffusivity_k;
    ScalarField& m_diffusivity_epsilon;
    ScalarField& m_source_k;
    ScalarField& m_source_epsilon;
    ScalarField& m_work;
    double m_relative_change = 0.0;
};

}  // 匿名命名空间

Model* makeKEpsilon(
    Case& problem,
    const VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity)
{
    return new KEpsilon(
        problem, velocity, face_flux, effective_viscosity,
        density, molecular_viscosity);
}

}  // babelsim::rans 命名空间
