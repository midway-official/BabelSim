#include "model.h"

#include <algorithm>
#include <cmath>

namespace babelsim::rans {
namespace {

// 标准正变量 Spalart-Allmaras。wallDistance 必须是 Case 提供的“到最近壁面
// 的真实几何距离”场，不能用网格线方向距离或最近单元中心距离替代。
class SpalartAllmaras final : public Model {
public:
    SpalartAllmaras(
        Case& problem,
        const VectorField& velocity,
        const ScalarField& face_flux,
        ScalarField& effective_viscosity,
        double density,
        double molecular_viscosity)
        : Model(problem, velocity, face_flux, effective_viscosity,
                density, molecular_viscosity),
          m_cb1(positiveSetting(problem.physics(), "saCb1", 0.1355)),
          m_cb2(positiveSetting(problem.physics(), "saCb2", 0.622)),
          m_sigma(positiveSetting(problem.physics(), "saSigma", 2.0 / 3.0)),
          m_kappa(positiveSetting(problem.physics(), "saKappa", 0.41)),
          m_cw2(positiveSetting(problem.physics(), "saCw2", 0.3)),
          m_cw3(positiveSetting(problem.physics(), "saCw3", 2.0)),
          m_cv1(positiveSetting(problem.physics(), "saCv1", 7.1)),
          m_ct3(positiveSetting(problem.physics(), "saCt3", 1.2)),
          m_ct4(positiveSetting(problem.physics(), "saCt4", 0.5)),
          m_cw1(positiveSetting(
              problem.physics(), "saCw1",
              m_cb1 / (m_kappa * m_kappa) + (1.0 + m_cb2) / m_sigma)),
          m_nu_tilda_min(positiveSetting(problem.physics(), "saNuTildaMin", 1e-14)),
          m_wall_distance_min(positiveSetting(
              problem.physics(), "saWallDistanceMin", 1e-12)),
          m_nu_tilda(problem.scalarField("nuTilda")),
          m_wall_distance(problem.scalarField("wallDistance")),
          m_mut(problem.scalarField("mut", 0.0)),
          m_previous_nu_tilda(problem.scalarField("ransPreviousNuTilda", 0.0)),
          m_nu_tilda_gradient(problem.vectorField("ransGradNuTilda", Vec3{})),
          m_vorticity_magnitude(problem.scalarField("ransVorticity", 0.0)),
          m_gradient_squared(problem.scalarField("ransGradNuTilda2", 0.0)),
          m_chi(problem.scalarField("ransChi", 0.0)),
          m_fv1(problem.scalarField("ransFv1", 0.0)),
          m_fv2(problem.scalarField("ransFv2", 0.0)),
          m_ft2(problem.scalarField("ransFt2", 0.0)),
          m_inverse_distance_squared(problem.scalarField("ransInvD2", 0.0)),
          m_stilda(problem.scalarField("ransSTilda", 0.0)),
          m_r(problem.scalarField("ransR", 0.0)),
          m_fw(problem.scalarField("ransFw", 0.0)),
          m_diffusivity(problem.scalarField("ransDiffusivityNuTilda", 0.0)),
          m_source(problem.scalarField("ransSourceNuTilda", 0.0)),
          m_work1(problem.scalarField("ransWork1", 0.0)),
          m_work2(problem.scalarField("ransWork2", 0.0))
    {
        boundField();
        updateKinematics();
        updateFunctions();
        updateViscosity();
        m_problem.output(m_mut);
    }

    const char* modelName() const override { return "Spalart-Allmaras"; }

    SolveResult correct() override {
        m_previous_nu_tilda.assign(m_nu_tilda);
        updateKinematics();
        updateEquationFields();
        const SolveResult result = solveTransport(
            m_nu_tilda, m_diffusivity, m_source);
        boundField();
        updateFunctions();
        updateViscosity();
        m_relative_change = diagnostics::relativeChange(
            m_nu_tilda, m_previous_nu_tilda);
        return result;
    }

    double relativeChange() const override { return m_relative_change; }

private:
    void boundField() {
        m_nu_tilda.evaluate(m_nu_tilda, [this](double value) {
            return std::max(value, m_nu_tilda_min);
        });
    }

    void updateFunctions() {
        const double molecular_nu = m_molecular_viscosity / m_density;
        const double cv1_cubed = m_cv1 * m_cv1 * m_cv1;
        m_vorticity_magnitude.evaluate(m_velocity_gradient, vorticityMagnitude);
        m_chi.evaluate(m_nu_tilda, [molecular_nu](double value) {
            return value / molecular_nu;
        });
        m_fv1.evaluate(m_chi, [cv1_cubed](double chi) {
            const double chi_cubed = chi * chi * chi;
            return chi_cubed / (chi_cubed + cv1_cubed);
        });
        m_work1.assignProduct(m_chi, m_fv1);
        m_work1.evaluate(m_work1, [](double value) { return 1.0 / (1.0 + value); });
        m_fv2.assignProduct(m_chi, m_work1);
        m_fv2.evaluate(m_fv2, [](double value) { return 1.0 - value; });
        m_ft2.evaluate(m_chi, [this](double chi) {
            return m_ct3 * std::exp(-m_ct4 * chi * chi);
        });
        m_inverse_distance_squared.evaluate(m_wall_distance, [this](double distance) {
            const double bounded = std::max(distance, m_wall_distance_min);
            return 1.0 / (bounded * bounded);
        });

        m_work1.assignProduct(m_nu_tilda, m_inverse_distance_squared);
        m_work1.assignProduct(m_fv2, m_work1);
        m_stilda.assign(m_vorticity_magnitude);
        m_stilda.addScaled(1.0 / (m_kappa * m_kappa), m_work1);
        m_stilda.evaluate(m_stilda, [](double value) {
            return std::max(value, 1e-30);
        });

        m_r.evaluate(m_stilda, [](double value) { return 1.0 / value; });
        m_r.assignProduct(m_inverse_distance_squared, m_r);
        m_r.assignProduct(m_nu_tilda, m_r);
        m_r.assignScaled(1.0 / (m_kappa * m_kappa), m_r);
        m_r.evaluate(m_r, [](double value) { return std::clamp(value, 0.0, 10.0); });
        m_fw.evaluate(m_r, [this](double r) {
            const double r6 = std::pow(r, 6.0);
            const double g = r + m_cw2 * (r6 - r);
            const double cw3_6 = std::pow(m_cw3, 6.0);
            return g * std::pow((1.0 + cw3_6) /
                                (std::pow(g, 6.0) + cw3_6), 1.0 / 6.0);
        });
    }

    void updateViscosity() {
        m_mut.assignProduct(m_nu_tilda, m_fv1);
        m_mut.assignScaled(m_density, m_mut);
        setEddyViscosity(m_mut);
    }

    void updateEquationFields() {
        updateFunctions();
        math::evaluate(math::grad(m_nu_tilda), m_nu_tilda_gradient);
        m_gradient_squared.evaluate(m_nu_tilda_gradient, [](const Vec3& gradient) {
            return squaredNorm(gradient);
        });

        // cb1(1-ft2) S~ nu~
        m_work1.evaluate(m_ft2, [](double value) { return 1.0 - value; });
        m_source.assignProduct(m_stilda, m_nu_tilda);
        m_source.assignProduct(m_work1, m_source);
        m_source.assignScaled(m_cb1, m_source);

        // -[cw1 fw - cb1 ft2/kappa^2] (nu~/d)^2
        m_work1.assignProduct(m_nu_tilda, m_nu_tilda);
        m_work1.assignProduct(m_inverse_distance_squared, m_work1);
        m_work2.assignScaled(m_cw1, m_fw);
        m_work2.addScaled(-m_cb1 / (m_kappa * m_kappa), m_ft2);
        m_source.addProduct(-1.0, m_work2, m_work1);
        m_source.addScaled(m_cb2 / m_sigma, m_gradient_squared);
        m_source.assignScaled(m_density, m_source);

        // 守恒形式中的扩散系数 rho*(nu+nu~)/sigma。
        m_diffusivity.fill(m_molecular_viscosity);
        m_diffusivity.addScaled(m_density, m_nu_tilda);
        m_diffusivity.assignScaled(1.0 / m_sigma, m_diffusivity);
    }

    double m_cb1;
    double m_cb2;
    double m_sigma;
    double m_kappa;
    double m_cw2;
    double m_cw3;
    double m_cv1;
    double m_ct3;
    double m_ct4;
    double m_cw1;
    double m_nu_tilda_min;
    double m_wall_distance_min;
    ScalarField& m_nu_tilda;
    ScalarField& m_wall_distance;
    ScalarField& m_mut;
    ScalarField& m_previous_nu_tilda;
    VectorField& m_nu_tilda_gradient;
    ScalarField& m_vorticity_magnitude;
    ScalarField& m_gradient_squared;
    ScalarField& m_chi;
    ScalarField& m_fv1;
    ScalarField& m_fv2;
    ScalarField& m_ft2;
    ScalarField& m_inverse_distance_squared;
    ScalarField& m_stilda;
    ScalarField& m_r;
    ScalarField& m_fw;
    ScalarField& m_diffusivity;
    ScalarField& m_source;
    ScalarField& m_work1;
    ScalarField& m_work2;
    double m_relative_change = 0.0;
};

}  // 匿名命名空间

Model* makeSpalartAllmaras(
    Case& problem,
    const VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity)
{
    return new SpalartAllmaras(
        problem, velocity, face_flux, effective_viscosity,
        density, molecular_viscosity);
}

}  // babelsim::rans 命名空间
