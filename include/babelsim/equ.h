#pragma once

#include "babelsim/field.h"
#include "babelsim/history.h"
#include "babelsim/methods.h"
#include "babelsim/solver_control.h"
#include <memory>

namespace babelsim::equ {
namespace detail { struct EquationAccess; }

// An integrated cell system Ax=b. Assembly calls freeze their inputs immediately.
// The field binding supplies layout, boundary conditions and the current iterate
// for deferred corrections. reset() retains allocation and removes contributions.
template<class T> class Equation {
public:
    explicit Equation(Field<T>& unknown);
    ~Equation();
    Equation(Equation&&) noexcept;
    Equation& operator=(Equation&&) noexcept;
    Equation(const Equation&) = delete;
    Equation& operator=(const Equation&) = delete;
    void reset(); // retain allocation and unknown binding, remove all terms
    Equation copy() const; // independent assembled coefficients, same unknown
    ScalarField diagonal() const; // integrated aP
    Field<T> rhs() const; // integrated b
    void reference(Index globalCell, double value); // scalar gauge at an explicit cell
    void referenceIfUnanchored(double value); // scalar constant-null-mode gauge
private:
    struct Storage;
    std::unique_ptr<Storage> storage_;
    friend struct detail::EquationAccess;
};
// Create an empty discrete equation for this unknown. No PDE terms are inferred
// from the field: subsequent assembly calls define the equation explicitly.
template<class T> Equation<T> createEquation(Field<T>& unknown) { return Equation<T>(unknown); }

// Operators add to the LHS, source adds to the RHS. For bound unknown x:
// div(eq, phi, c) adds c*div(phi*x); phi is an oriented, integrated face flux and
// becomes the boundary-flux context of the equation: a second, different phi is
// rejected instead of silently overwriting the unknown's boundary traces.
// laplacian adds multiplier*div(coefficient*grad(x)); -1 is usual LHS diffusion.
template<class T> void div(Equation<T>&, const ScalarField& flux, double scale = 1.0);
template<class T> void laplacian(Equation<T>&, double coefficient, double multiplier);
template<class T> void laplacian(Equation<T>&, const ScalarField& coefficient, double multiplier);
template<class T> void reaction(Equation<T>&, double coefficient);
template<class T> void reaction(Equation<T>&, const ScalarField& coefficient, double scale = 1.0);
template<class T> void source(Equation<T>&, T value);
template<class T> void source(Equation<T>&, const Field<T>& value, double scale = 1.0);

// Nonconservative capacity*d(x)/dt. BDF2 requires two explicit histories;
// variable step BDF2 uses previous_dt. No implicit history registration occurs.
template<class T> void ddt(Equation<T>&, double capacity, const Field<T>& previous,
    double dt, TimeMethod method = TimeMethod::Euler,
    const Field<T>* older = nullptr, double previous_dt = 0.0);
template<class T> void ddt(Equation<T>&, const ScalarField& capacity, const Field<T>& previous,
    double dt, TimeMethod method = TimeMethod::Euler,
    const Field<T>* older = nullptr, double previous_dt = 0.0);

// Configured Euler/BDF2 with explicit saved histories; first BDF2 step uses Euler.
template<class T> void ddt(Equation<T>&,double capacity,const time::History<T>&);
template<class T> void ddt(Equation<T>&,const ScalarField& capacity,const time::History<T>&);

// Algebra operates on already integrated coefficients. No volume factor here.
template<class T> void scale(Equation<T>&, double);
template<class T> void add(Equation<T>&, const Equation<T>&, double factor = 1.0);
template<class T> void addDiagonal(Equation<T>&, const ScalarField&);
template<class T> void addRhs(Equation<T>&, const Field<T>&);
template<class T> Field<T> apply(const Equation<T>&, const Field<T>&);
template<class T> Field<T> residual(const Equation<T>&, const Field<T>&); // b-Ax
// Canonical relaxation: diagonal/=alpha; b+=(new-old diagonal)*previous.
template<class T> void relax(Equation<T>&, const Field<T>& previous, double alpha);
// Sum of scalar diffusion LHS face fluxes, using the gradients and boundary
// data frozen during assembly. This is consistent with deferred corrections.
ScalarField faceFlux(const Equation<double>&, const ScalarField& solution);
template<class T> SolveResult solve(const Equation<T>&, Field<T>&);
template<class T> SolveResult solve(const Equation<T>&, Field<T>&, const LinearSolverConfig&);
// The equation already owns the binding to its unknown.  These overloads
// prevent a second target argument from accidentally disagreeing with it.
template<class T> SolveResult solve(const Equation<T>&);
template<class T> SolveResult solve(const Equation<T>&, const LinearSolverConfig&);

extern template class Equation<double>;
extern template class Equation<Vec3>;
} // namespace babelsim::equ

namespace babelsim::diagnostics {
// Observe the assembled system without reassembly or mutation: ||b-Ax|| /
// max(||Ax|| + ||b||, 1e-30). For nonlinear convergence call before relaxation.
template<class T> double relativeResidual(const equ::Equation<T>&, const Field<T>&);
}
