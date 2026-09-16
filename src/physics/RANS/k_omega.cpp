#include "api.h"
#include "babelsim/equ.h"

#include <algorithm>

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



// Wilcox 1988m 高雷诺数 k-omega。该明确命名的版本不含 2006 版交叉扩散，
// P 使用不可压缩线性涡黏性近似 mu_t * 2*S':S'。
class KOmega final : public Model {
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
    KOmega(
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
          m_beta_star(problem.physics().positive("kOmegaBetaStar", 0.09)),
          m_beta(problem.physics().positive("kOmegaBeta", 3.0 / 40.0)),
          m_gamma(problem.physics().positive("kOmegaGamma", 5.0 / 9.0)),
          m_sigma_k(problem.physics().positive("kOmegaSigmaK", 0.5)),
          m_sigma_omega(problem.physics().positive("kOmegaSigmaOmega", 0.5)),
          m_k_min(problem.physics().positive("kMin", 1e-12)),
          m_omega_min(problem.physics().positive("omegaMin", 1e-12)),
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
          m_sink_k(problem.scalarField("ransSinkK", 0.0)),
          m_sink_second(problem.scalarField("ransSinkSecond", 0.0)),
          m_work(problem.scalarField("ransWork", 0.0))
    {
        m_effective_viscosity.useCalculatedBoundary();
        m_velocity_gradient.useCalculatedBoundary();
        m_strain_measure.useCalculatedBoundary();
        m_linear_options=readLinearControl(problem, m_k);
        for (ScalarField* field : {&m_mut, &m_inverse_k, &m_inverse_omega, &m_production, &m_diffusivity_k, &m_diffusivity_omega, &m_source_k, &m_source_omega, &m_sink_k, &m_sink_second, &m_work})
            field->useCalculatedBoundary();
        boundFields();
        updateKinematics();
        updateViscosity();
        m_problem.output(m_mut);

    }

    const char* modelName() const override { return "Wilcox1988m k-omega"; }

    double tolerance() const override {return m_tolerance;}
    double relativeResidual() const override {return m_relative_residual;}
    void saveOld(double dt) override {
        math::saveOld(m_k_old, m_k, dt);
        math::saveOld(m_omega_old, m_omega, dt);
    }

    SolveResult correct() override {
        m_previous_k.assign(m_k);
        m_previous_omega.assign(m_omega);
        updateKinematics();
        updateSourcesAndDiffusivities();

        auto kEquation=assemble_k();
        equ::relax(kEquation,m_k,m_relaxation);
        const auto k_result=equ::solve(kEquation,m_k,m_linear_options);
        auto omegaEquation=assemble_omega();
        equ::relax(omegaEquation,m_omega,m_relaxation);
        const auto omega_result=equ::solve(omegaEquation,m_omega,m_linear_options);
        boundFields();
        updateSourcesAndDiffusivities();
        m_relative_residual = std::max(
            equ::relativeResidual(assemble_k(),m_k),
            equ::relativeResidual(assemble_omega(),m_omega));
        m_relative_change = std::max(
            diagnostics::relativeChange(m_k, m_previous_k),
            diagnostics::relativeChange(m_omega, m_previous_omega));
        return combine(k_result, omega_result);
    }

    double relativeChange() const override { return m_relative_change; }

private:
    equ::Matrix<double> assemble_k() {
        auto A=equ::createEquation(m_k);
        equ::ddt(A,m_density,m_k_old);
        equ::div(A,m_face_flux,m_density);
        equ::laplacian(A,m_diffusivity_k,-1.0);
        equ::reaction(A,m_sink_k);
        equ::source(A,m_source_k);
        return A;
    }
    equ::Matrix<double> assemble_omega() {
        auto A=equ::createEquation(m_omega);
        equ::ddt(A,m_density,m_omega_old);
        equ::div(A,m_face_flux,m_density);
        equ::laplacian(A,m_diffusivity_omega,-1.0);
        equ::reaction(A,m_sink_second);
        equ::source(A,m_source_omega);
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
    static SolveResult combine(const SolveResult& first, const SolveResult& second) {
    SolveResult result;
    result.status = first.status == SolveStatus::NumericalFailure ||
            second.status == SolveStatus::NumericalFailure
        ? SolveStatus::NumericalFailure
        : first.converged() && second.converged()
            ? SolveStatus::Converged : SolveStatus::MaxIterations;
    result.iterations = first.iterations + second.iterations;
    result.initial_residual = std::hypot(first.initial_residual, second.initial_residual);
    result.final_residual = std::hypot(first.final_residual, second.final_residual);
    result.relative_residual = std::max(first.relative_residual, second.relative_residual);
    return result;
}

    void boundFields() {
        m_k.setBoundaryFlux(m_face_flux);
        m_omega.setBoundaryFlux(m_face_flux);
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
        m_sink_k.assignScaled(m_beta_star * m_density, m_omega);

        m_inverse_k.evaluate(m_k, [this](double value) {
            return 1.0 / std::max(value, m_k_min);
        });
        m_work.assignProduct(m_production, m_omega);
        m_work.assignProduct(m_inverse_k, m_work);
        m_source_omega.assignScaled(m_gamma, m_work);
        m_work.assignProduct(m_omega, m_omega);
        m_sink_second.assignScaled(m_beta * m_density, m_omega);

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
    ScalarField& m_sink_k;
    ScalarField& m_sink_second;
    ScalarField& m_work;
    double m_relative_change = 0.0;
    math::History<double> m_k_old{m_k};
    math::History<double> m_omega_old{m_omega};

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
