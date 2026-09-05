#pragma once

#include "babelsim/case.h"
#include "babelsim/solver.h"

namespace babelsim::rans {

// RANS 模型的私有扩展点。模型负责自己的输运变量、源项和湍流黏度；
// SIMPLE 只读取由模型更新的有效动力黏度，不了解具体闭合方程。
class Model {
public:
    Model(
        Case& problem,
        const VectorField& velocity,
        const ScalarField& face_flux,
        ScalarField& effective_viscosity,
        double density,
        double molecular_viscosity);
    virtual ~Model() = default;

    virtual const char* modelName() const = 0;
    virtual SolveResult correct() = 0;
    virtual double relativeChange() const = 0;
    double tolerance() const { return m_tolerance; }

protected:
    SolveResult solveTransport(
        ScalarField& variable,
        const ScalarField& diffusivity,
        const ScalarField& source) const;
    void updateKinematics();
    void setEddyViscosity(const ScalarField& turbulent_viscosity);
    static SolveResult combine(const SolveResult& first, const SolveResult& second);

    Case& m_problem;
    const VectorField& m_velocity;
    const ScalarField& m_face_flux;
    ScalarField& m_effective_viscosity;
    double m_density;
    double m_molecular_viscosity;
    double m_relaxation;
    double m_tolerance;
    TensorField& m_velocity_gradient;
    ScalarField& m_strain_measure;
};

double positiveSetting(const Parameters& settings, const char* key, double fallback);
double strainMeasure(const Tensor3& gradient);
double vorticityMagnitude(const Tensor3& gradient);

Model* makeSpalartAllmaras(
    Case& problem,
    const VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity);
Model* makeKOmega(
    Case& problem,
    const VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity);
Model* makeKEpsilon(
    Case& problem,
    const VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity);

}  // babelsim::rans 命名空间
