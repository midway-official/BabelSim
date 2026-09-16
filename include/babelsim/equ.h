#pragma once

#include "babelsim/field.h"
#include "babelsim/history.h"
#include "babelsim/methods.h"
#include "babelsim/solver_control.h"
#include <memory>

namespace babelsim::equ {
namespace detail { struct MatrixAccess; }

// An integrated cell system Ax=b. Assembly calls freeze their inputs immediately.
// The field binding supplies layout, boundary conditions and the current iterate
// for deferred corrections. clear() retains allocation and removes contributions.
template<class T> class Matrix {
public:
    explicit Matrix(Field<T>& unknown);
    ~Matrix();
    Matrix(Matrix&&) noexcept;
    Matrix& operator=(Matrix&&) noexcept;
    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;
private:
    struct Storage;
    std::unique_ptr<Storage> storage_;
    friend struct detail::MatrixAccess;
};
template<class T> Matrix<T> matrix(Field<T>& unknown) { return Matrix<T>(unknown); }
// Create an empty discrete equation for this unknown. No PDE terms are inferred
// from the field: subsequent assembly calls define the equation explicitly.
template<class T> Matrix<T> createEquation(Field<T>& unknown) { return Matrix<T>(unknown); }
template<class T> void clear(Matrix<T>&);
// Remove all assembled terms; retain the unknown binding and allocated storage.
// This does not initialize or change the unknown field's values.
template<class T> void reset(Matrix<T>& equation) { clear(equation); }
template<class T> Matrix<T> copy(const Matrix<T>&);

// Operators add to the LHS, source adds to the RHS. laplacian has its
// mathematical sign: diffusion on the LHS normally uses scale=-1.
template<class T> void div(Matrix<T>&, const ScalarField& flux, double scale = 1.0);
template<class T> void laplacian(Matrix<T>&, double coefficient, double scale);
template<class T> void laplacian(Matrix<T>&, const ScalarField& coefficient, double scale);
template<class T> void reaction(Matrix<T>&, double coefficient);
template<class T> void reaction(Matrix<T>&, const ScalarField& coefficient, double scale = 1.0);
template<class T> void source(Matrix<T>&, T value);
template<class T> void source(Matrix<T>&, const Field<T>& value, double scale = 1.0);

// Nonconservative capacity*d(x)/dt. BDF2 requires two explicit histories;
// variable step BDF2 uses previous_dt. No implicit history registration occurs.
template<class T> void ddt(Matrix<T>&, double capacity, const Field<T>& previous,
    double dt, TimeMethod method = TimeMethod::Euler,
    const Field<T>* older = nullptr, double previous_dt = 0.0);
template<class T> void ddt(Matrix<T>&, const ScalarField& capacity, const Field<T>& previous,
    double dt, TimeMethod method = TimeMethod::Euler,
    const Field<T>* older = nullptr, double previous_dt = 0.0);

// Configured Euler/BDF2 with explicit saved histories; first BDF2 step uses Euler.
template<class T> void ddt(Matrix<T>&,double capacity,const math::History<T>&);
template<class T> void ddt(Matrix<T>&,const ScalarField& capacity,const math::History<T>&);

// Algebra operates on already integrated coefficients. No volume factor here.
template<class T> void scale(Matrix<T>&, double);
template<class T> void add(Matrix<T>&, const Matrix<T>&, double factor = 1.0);
template<class T> void addDiagonal(Matrix<T>&, const ScalarField&);
template<class T> void addRhs(Matrix<T>&, const Field<T>&);
template<class T> ScalarField diagonal(const Matrix<T>&);
template<class T> ScalarField response(const Matrix<T>&); // V/aP of this matrix
template<class T> Field<T> rhs(const Matrix<T>&);
template<class T> Field<T> apply(const Matrix<T>&, const Field<T>&);
template<class T> Field<T> residual(const Matrix<T>&, const Field<T>&); // b-Ax
// Canonical relaxation: diagonal/=alpha; b+=(new-old diagonal)*previous.
template<class T> void relax(Matrix<T>&, const Field<T>& previous, double alpha);
void reference(Matrix<double>&, Index global_cell, double value);
// Supply a gauge only if the assembled operator leaves the constant mode free.
// An already anchored system is unchanged. This is not a general rank detector.
void reference(Matrix<double>&, double value);
// Sum of scalar diffusion LHS face fluxes, using the gradients and boundary
// data frozen during assembly. This is consistent with deferred corrections.
ScalarField faceFlux(const Matrix<double>&, const ScalarField& solution);
template<class T> double relativeResidual(const Matrix<T>&, const Field<T>&);
template<class T> SolveResult solve(const Matrix<T>&, Field<T>&);
template<class T> SolveResult solve(const Matrix<T>&, Field<T>&, const LinearSolverConfig&);

extern template class Matrix<double>;
extern template class Matrix<Vec3>;
} // namespace babelsim::equ
