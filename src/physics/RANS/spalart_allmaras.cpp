#include "api.h"
#include "babelsim/equ.h"

#include <algorithm>
#include <cmath>

namespace babelsim::rans {
namespace {
double strainMeasure(const Tensor3& gradient) {
    const double divergence = gradient[0][0] + gradient[1][1] + gradient[2][2];
    double squared = 0.0;
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            double value = 0.5 * (gradient[row][column] + gradient[column][row]);
            if (row == column) value -= divergence / 3.0;
            squared += value * value;
        }
    }
    return 2.0 * squared;
}

double vorticityMagnitude(const Tensor3& gradient) {
    double squared = 0.0;
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            const double value = 0.5 *
                (gradient[row][column] - gradient[column][row]);
            squared += value * value;
        }
    }
    return std::sqrt(2.0 * squared);
}



// 标准正变量 Spalart-Allmaras。wallDistance 必须是 Case 提供的“到最近壁面
// 的真实几何距离”场，不能用网格线方向距离或最近单元中心距离替代。
class SpalartAllmaras final : public Model {
    // This model owns its numerical state and transport implementation.
    Case& m_problem;
    const VectorField& m_velocity;
    const ScalarField& m_face_flux;
    ScalarField& m_effective_viscosity;
    double m_density, m_molecular_viscosity, m_relaxation, m_tolerance;
    TensorField& m_velocity_gradient;
    ScalarField& m_strain_measure;
    LinearSolverConfig m_linear_options;
    double m_relative_residual=std::numeric_limits<double>::infinity();
public:
    SpalartAllmaras(
        Case& problem,
        const VectorField& velocity,
        const ScalarField& face_flux,
        ScalarField& effective_viscosity,
        double density,
        double molecular_viscosity)
        : m_problem(problem), m_velocity(velocity), m_face_flux(face_flux),
          m_effective_viscosity(effective_viscosity), m_density(density),
          m_molecular_viscosity(molecular_viscosity),
          m_relaxation(problem.physics().fraction("turbulenceRelaxation",0.7)),
          m_tolerance(problem.physics().positive("turbulenceTolerance",1e-6)),
          m_velocity_gradient(problem.tensorField("ransGradU",Tensor3{})),
          m_strain_measure(problem.scalarField("ransStrain2",0.0)),
          m_cb1(problem.physics().positive("saCb1", 0.1355)),
          m_cb2(problem.physics().positive("saCb2", 0.622)),
          m_sigma(problem.physics().positive("saSigma", 2.0 / 3.0)),
          m_kappa(problem.physics().positive("saKappa", 0.41)),
          m_cw2(problem.physics().positive("saCw2", 0.3)),
          m_cw3(problem.physics().positive("saCw3", 2.0)),
          m_cv1(problem.physics().positive("saCv1", 7.1)),
          m_ct3(problem.physics().positive("saCt3", 1.2)),
          m_ct4(problem.physics().positive("saCt4", 0.5)),
          m_cw1(problem.physics().positive(
              "saCw1",
              m_cb1 / (m_kappa * m_kappa) + (1.0 + m_cb2) / m_sigma)),
          m_nu_tilda_min(problem.physics().positive("saNuTildaMin", 1e-14)),
          m_wall_distance_min(problem.physics().positive(
              "saWallDistanceMin", 1e-12)),
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
          m_sink(problem.scalarField("ransSinkNuTilda", 0.0)),
          m_work1(problem.scalarField("ransWork1", 0.0)),
          m_work2(problem.scalarField("ransWork2", 0.0))
    {
        m_effective_viscosity.useCalculatedBoundary();
        m_velocity_gradient.useCalculatedBoundary();
        m_strain_measure.useCalculatedBoundary();
        m_linear_options=readLinearControl(problem, m_nu_tilda);
        for (ScalarField* field : {&m_mut, &m_vorticity_magnitude, &m_gradient_squared, &m_chi, &m_fv1, &m_fv2, &m_ft2, &m_inverse_distance_squared, &m_stilda, &m_r, &m_fw, &m_diffusivity, &m_source, &m_sink, &m_work1, &m_work2})
            field->useCalculatedBoundary();
        m_nu_tilda_gradient.useCalculatedBoundary();
        boundField();
        updateKinematics();
        updateFunctions();
        updateViscosity();
        m_problem.output(m_mut);

    }

    const char* modelName() const override { return "Spalart-Allmaras"; }

    double tolerance() const override {return m_tolerance;}
    double relativeResidual() const override {return m_relative_residual;}
    void saveOld(double dt) override {
        math::saveOld(m_nu_tilda_old, m_nu_tilda, dt);
    }

    SolveResult correct() override {
        m_previous_nu_tilda.assign(m_nu_tilda);
        updateKinematics();
        updateEquationFields();
        auto nu_tildaEquation=assemble_nu_tilda();
        equ::relax(nu_tildaEquation,m_nu_tilda,m_relaxation);
        const auto result=equ::solve(nu_tildaEquation,m_nu_tilda,m_linear_options);
        boundField();
        updateEquationFields();
        updateViscosity();
        m_relative_residual = equ::relativeResidual(assemble_nu_tilda(),m_nu_tilda);
        m_relative_change = diagnostics::relativeChange(
            m_nu_tilda, m_previous_nu_tilda);
        return result;
    }

    double relativeChange() const override { return m_relative_change; }

private:
    equ::Matrix<double> assemble_nu_tilda() {
        auto A=equ::createEquation(m_nu_tilda);
        equ::ddt(A,m_density,m_nu_tilda_old);
        equ::div(A,m_face_flux,m_density);
        equ::laplacian(A,m_diffusivity,-1.0);
        equ::reaction(A,m_sink);
        equ::source(A,m_source);
        return A;
    }
    void updateKinematics() {
        const_cast<VectorField&>(m_velocity).setBoundaryFlux(m_face_flux);
        m_velocity_gradient=math::grad(m_velocity);
        m_strain_measure.evaluate(m_velocity_gradient,strainMeasure);
    }
    void setEddyViscosity(const ScalarField& turbulent) {
        m_effective_viscosity.fill(m_molecular_viscosity);
        m_effective_viscosity+=turbulent;
    }
    void boundField() {
        m_nu_tilda.setBoundaryFlux(m_face_flux);
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
        m_nu_tilda_gradient = math::grad(m_nu_tilda);
        m_gradient_squared.evaluate(m_nu_tilda_gradient, [](const Vec3& gradient) {
            return squaredNorm(gradient);
        });

        // Split the signed reaction coefficient without changing the positive-variable PDE.
        // a = cb1*(1-ft2)*S~ - (cw1*fw-cb1*ft2/kappa^2)*nu~/d^2.
        m_work1.evaluate(m_ft2, [](double value) { return 1.0 - value; });
        m_source.assignProduct(m_stilda, m_work1);
        m_source.assignScaled(m_cb1, m_source);
        m_work1.assignProduct(m_nu_tilda, m_inverse_distance_squared);
        m_work2.assignScaled(m_cw1, m_fw);
        m_work2.addScaled(-m_cb1 / (m_kappa * m_kappa), m_ft2);
        m_source.addProduct(-1.0, m_work2, m_work1);
        m_sink.evaluate(m_source, [this](double a) { return m_density * std::max(-a, 0.0); });
        m_source.evaluate(m_source, [](double a) { return std::max(a, 0.0); });
        m_source.assignProduct(m_nu_tilda, m_source);
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
    ScalarField& m_sink;
    ScalarField& m_work1;
    ScalarField& m_work2;
    double m_relative_change = 0.0;
    math::History<double> m_nu_tilda_old{m_nu_tilda};

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
