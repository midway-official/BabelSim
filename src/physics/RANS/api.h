#pragma once
#include "babelsim/case.h"
#include "babelsim/solver.h"
#include <memory>
namespace babelsim::rans {
class Model {
public:
    virtual ~Model()=default;
    virtual void saveOld(double dt)=0;
    virtual const char* modelName() const=0;
    virtual SolveResult correct()=0;
    virtual double relativeChange() const=0;
    virtual double relativeResidual() const=0;
    virtual double tolerance() const=0;
};
Model* create(Case&,const VectorField&,const ScalarField&,ScalarField&,double,double);
void destroy(Model*) noexcept;
using Handle=std::unique_ptr<Model,void(*)(Model*)>;
inline Handle load(Case& c,const VectorField& u,const ScalarField& phi,ScalarField& mu,double rho,double molecular) {
    return Handle(create(c,u,phi,mu,rho,molecular),destroy);
}
// Coupling object owns the wiring; each selected implementation remains independent.
class Turbulence {
public:
    Turbulence(Case& problem,const VectorField& velocity,const ScalarField& phi)
        : viscosity_(&problem.scalarField("muEffective",problem.physics().positive("dynamicViscosity"))),
          model_(create(problem,velocity,phi,*viscosity_,problem.physics().positive("density"),
                        problem.physics().positive("dynamicViscosity")),destroy) {}
    explicit operator bool() const {return bool(model_);}
    const ScalarField& viscosity() const {return *viscosity_;}
    void saveOld(double dt) {if(model_) model_->saveOld(dt);}
    SolveResult correct() {
        return model_ ? model_->correct() : SolveResult{SolveStatus::Converged,0,0,0,0};
    }
    double relativeChange() const {return model_ ? model_->relativeChange() : 0.0;}
    double relativeResidual() const {return model_ ? model_->relativeResidual() : 0.0;}
    double tolerance() const {return model_ ? model_->tolerance() : 0.0;}
private:
    ScalarField* viscosity_;
    Handle model_;
};
inline Turbulence load(Case& problem,const VectorField& velocity,const ScalarField& phi) {
    return Turbulence(problem,velocity,phi);
}

void saveOld(Model&,double dt);
SolveResult correct(Model&);
double relativeChange(const Model&);
double relativeResidual(const Model&);
double tolerance(const Model&);
const char* name(const Model&);
}
