#pragma once

#include "babelsim/field.h"
#include "babelsim/methods.h"

namespace babelsim::math {

// Known-field arithmetic. Public functions below materialize independent results.
// The operation records are the internal execution bridge, not deferred user
// expressions. Whole-field operations synchronize their inputs and outputs.
struct ScalarGradient {
    const ScalarField& field;
    OperatorOptions options{};
};

// 每个面的外法向梯度；执行时自动重构单元梯度、同步输入并修正非正交性。
struct NormalGradient {
    const ScalarField& field;
    // Optional gradient already computed by the caller.  Supplying it keeps
    // compositions such as Rhie--Chow explicit without changing the
    // mathematical definition of the face-normal derivative.
    const VectorField* gradient = nullptr;
    OperatorOptions options{};
};

struct VectorGradient {
    const VectorField& field;
    OperatorOptions options{};
};

struct FaceFlux {
    // cell 输入先插值；face 输入直接计算 Sf·value，不重复插值。
    const VectorField& velocity;
    OperatorOptions options{};
};

// 几何选择，不是 MPI 分区选择。Interior 包含跨分区的内部面，但保持物理边界值。
enum class FaceRegion { All, Interior };

// 标量扩散面通量：coefficient * Sf·grad(field)。coefficient 可位于 cell 或 face，
// FVM 执行层负责同步、插值和梯度工作区。
struct ScalarDiffusionFlux {
    const ScalarField& coefficient;
    const ScalarField& field;
    const VectorField* gradient = nullptr;
    OperatorOptions options{};
};

struct FaceDivergence {
    const ScalarField& flux;
};

struct VectorDivergence {
    const VectorField& field;
    OperatorOptions options{};
};

// 张量按 rows[i][j] 存储；返回 (div T)_i = d T_ij / dx_j。
struct TensorDivergence {
    const TensorField& field;
    OperatorOptions options{};
};

struct ScalarConvection {
    const ScalarField& flux;
    const ScalarField& field;
    OperatorOptions options{};
};

struct VectorConvection {
    const ScalarField& flux;
    const VectorField& field;
    OperatorOptions options{};
};

struct ScalarInterpolation {
    const ScalarField& field;
    OperatorOptions options{};
};

struct VectorInterpolation {
    const VectorField& field;
    OperatorOptions options{};
};

struct ScalarReconstruction {
    const ScalarField& field;
    const VectorField& gradient;
};

struct VectorReconstruction {
    const VectorField& field;
    const TensorField& gradient;
};

struct ScalarLaplacian {
    const ScalarField& field;
    double coefficient = 1.0;
    const ScalarField* coefficient_field = nullptr;
    OperatorOptions options{};
};

void evaluate(ScalarGradient operation, VectorField& result);
void evaluate(NormalGradient operation, ScalarField& result);
void evaluate(ScalarDiffusionFlux operation, ScalarField& result);
void evaluate(VectorGradient operation, TensorField& result);
void evaluate(FaceFlux operation, ScalarField& result);
void evaluate(FaceDivergence operation, ScalarField& result);
void evaluate(VectorDivergence operation, ScalarField& result);
void evaluate(TensorDivergence operation, VectorField& result);
void evaluate(ScalarConvection operation, ScalarField& result);
void evaluate(VectorConvection operation, VectorField& result);
void evaluate(ScalarInterpolation operation, ScalarField& result);
void evaluate(VectorInterpolation operation, VectorField& result);
void evaluate(ScalarReconstruction operation, ScalarField& result);
void evaluate(VectorReconstruction operation, VectorField& result);
void evaluate(ScalarLaplacian operation, ScalarField& result);
void subtract(
    const ScalarField& coefficient,
    ScalarGradient operation,
    VectorField& target);
// target 是已有通量；增量只更新所选区域，输出仍履行完整的整场同步契约。
void add(FaceFlux operation, ScalarField& target, FaceRegion region = FaceRegion::All);
void subtract(ScalarDiffusionFlux operation, ScalarField& target,
              FaceRegion region = FaceRegion::All);

// Eager mathematical operations: auto g=grad(p) is an independent snapshot.
template<class R, class Op>
Field<R> computed(const Mesh& mesh, FieldLocation location, const char* name, Op op) {
    Field<R> result(mesh,location,name);
    result.useCalculatedBoundary();
    evaluate(op,result);
    return result;
}
template<class T> void evaluate(const Field<T>& value, Field<T>& output) { output=value; }
inline VectorField grad(const ScalarField& f, const OperatorOptions& options = {}) { return computed<Vec3>(f.mesh(),FieldLocation::Cell,"grad",ScalarGradient{f,options}); }
inline TensorField grad(const VectorField& f, const OperatorOptions& options = {}) { return computed<Tensor3>(f.mesh(),FieldLocation::Cell,"grad",VectorGradient{f,options}); }
inline ScalarField normalGradient(const ScalarField& f, const OperatorOptions& options = {}) {
    return computed<double>(f.mesh(),FieldLocation::Face,"normalGradient",NormalGradient{f,nullptr,options});
}
// The second overload reuses a caller-owned cell gradient.  It still evaluates
// the same face-normal derivative, including the configured non-orthogonal
// correction and boundary treatment.
inline ScalarField normalGradient(const ScalarField& f, const VectorField& gradient, const OperatorOptions& options = {}) {
    return computed<double>(
        f.mesh(), FieldLocation::Face, "normalGradient", NormalGradient{f, &gradient,options});
}
// Oriented, area-integrated face flux: vector value dot Sf. A cell vector is
// interpolated first; an already face-centred vector is used directly.
// inletOutlet uses the input's saved flux context (owner velocity on first use).
// The resulting flux becomes the input's next boundary context.
inline ScalarField flux(const VectorField& f, const OperatorOptions& options = {}) { return computed<double>(f.mesh(),FieldLocation::Face,"flux",FaceFlux{f,options}); }
// Positive mathematical diffusive flux k*grad(f).Sf, NOT -k*grad(f).Sf.
// k may be cell- or face-centred. f must be cell-centred. The optional cell
// gradient supplies the deferred nonorthogonal correction; no Rhie-Chow here.
inline ScalarField flux(const ScalarField& k,const ScalarField& f, const OperatorOptions& options = {}) { return computed<double>(f.mesh(),FieldLocation::Face,"diffusionFlux",ScalarDiffusionFlux{k,f,nullptr,options}); }
inline ScalarField flux(const ScalarField& k,const ScalarField& f,const VectorField& gradient, const OperatorOptions& options = {}) { return computed<double>(f.mesh(),FieldLocation::Face,"diffusionFlux",ScalarDiffusionFlux{k,f,&gradient,options}); }
// Scalar input is an oriented integrated FACE flux: sum(outward flux)/cell V.
// A scalar cell field is rejected; it is not silently interpreted as a flux.
inline ScalarField div(const ScalarField& f) { return computed<double>(f.mesh(),FieldLocation::Cell,"div",FaceDivergence{f}); }
inline ScalarField div(const VectorField& f, const OperatorOptions& options = {}) { return computed<double>(f.mesh(),FieldLocation::Cell,"div",VectorDivergence{f,options}); }
inline VectorField div(const TensorField& f, const OperatorOptions& options = {}) { return computed<Vec3>(f.mesh(),FieldLocation::Cell,"div",TensorDivergence{f,options}); }
// Explicit div(phi*f): phi is face flux, f is known cell data. Unlike equ::div,
// these overloads evaluate a field and do not bind or solve for an unknown.
inline ScalarField div(const ScalarField& phi,const ScalarField& f, const OperatorOptions& options = {}) { return computed<double>(f.mesh(),FieldLocation::Cell,"div",ScalarConvection{phi,f,options}); }
inline VectorField div(const ScalarField& phi,const VectorField& f, const OperatorOptions& options = {}) { return computed<Vec3>(f.mesh(),FieldLocation::Cell,"div",VectorConvection{phi,f,options}); }
// Cell-to-face interpolation using the configured interpolation scheme only.
inline ScalarField interpolate(const ScalarField& f, const OperatorOptions& options = {}) { return computed<double>(f.mesh(),FieldLocation::Face,"interpolate",ScalarInterpolation{f,options}); }
inline VectorField interpolate(const VectorField& f, const OperatorOptions& options = {}) { return computed<Vec3>(f.mesh(),FieldLocation::Face,"interpolate",VectorInterpolation{f,options}); }
inline ScalarField reconstruct(const ScalarField& f,const VectorField& g) { return computed<double>(f.mesh(),FieldLocation::Face,"reconstruct",ScalarReconstruction{f,g}); }
inline VectorField reconstruct(const VectorField& f,const TensorField& g) { return computed<Vec3>(f.mesh(),FieldLocation::Face,"reconstruct",VectorReconstruction{f,g}); }
inline ScalarField laplacian(const ScalarField& f, const OperatorOptions& options = {}) { return computed<double>(f.mesh(),FieldLocation::Cell,"laplacian",ScalarLaplacian{f,1.0,nullptr,options}); }
inline ScalarField laplacian(double k,const ScalarField& f, const OperatorOptions& options = {}) { return computed<double>(f.mesh(),FieldLocation::Cell,"laplacian",ScalarLaplacian{f,k,nullptr,options}); }
inline ScalarField laplacian(const ScalarField& k,const ScalarField& f, const OperatorOptions& options = {}) { return computed<double>(f.mesh(),FieldLocation::Cell,"laplacian",ScalarLaplacian{f,1.0,&k,options}); }
void add(const ScalarField& increment,ScalarField& target,FaceRegion region=FaceRegion::All);
void subtract(const ScalarField& increment,ScalarField& target,FaceRegion region=FaceRegion::All);
void subtract(const ScalarField& coefficient,const VectorField& gradient,VectorField& target);

}  // babelsim::math 命名空间

namespace babelsim {
// Eager pointwise expressions. Derived cell fields carry computed boundary
// traces; assigning them into an unknown preserves that unknown's constraints.
template<class T, class Function>
auto fieldUnary(const Field<T>& field, Function function) {
    using R = decltype(function(T{}));
    Field<R> result(field.mesh(), field.location(), "math.result");
    result.useCalculatedBoundary();
    result.evaluate(field, function);
    return result;
}
template<class A, class B, class Function>
auto fieldBinary(const Field<A>& a, const Field<B>& b, Function fn) {
    using R = decltype(fn(A{}, B{}));
    Field<R> result(a.mesh(), a.location(), "math.result");
    result.useCalculatedBoundary(); result.evaluate(a,b,fn); return result;
}
template<class T> Field<T> operator+(const Field<T>& a,const Field<T>& b) {
    return fieldBinary(a,b,[](T x,T y){return x+y;});
}
template<class T> Field<T> operator-(const Field<T>& a,const Field<T>& b) {
    return fieldBinary(a,b,[](T x,T y){return x-y;});
}
template<class T> Field<T> operator*(double a,const Field<T>& b) {
    Field<T> r(b.mesh(),b.location(),"math.scale"); r.useCalculatedBoundary();
    r.evaluate(b,[a](T v){return a*v;}); return r;
}
template<class T> Field<T> operator*(const Field<T>& b,double a) {return a*b;}
template<class T> Field<T> operator-(const Field<T>& b) {return -1.0*b;}
template<class T> Field<T> operator/(const Field<T>& b,double a) {return (1.0/a)*b;}
template<class T> Field<T> operator*(const ScalarField& a,const Field<T>& b) {
    return fieldBinary(a,b,[](double x,T y){return x*y;});
}
inline ScalarField operator/(const ScalarField& a,const ScalarField& b) {
    return fieldBinary(a,b,[](double x,double y){return x/y;});
}
inline ScalarField operator+(const ScalarField& field, double value) {
    return fieldUnary(field, [value](double x) { return x + value; });
}
inline ScalarField operator+(double value, const ScalarField& field) { return field + value; }
inline ScalarField operator-(const ScalarField& field, double value) { return field + (-value); }
inline ScalarField operator-(double value, const ScalarField& field) { return value + (-field); }
inline ScalarField operator/(double value, const ScalarField& field) {
    return fieldUnary(field, [value](double x) { return value / x; });
}
namespace math {
// Apply a local mathematical function, including to computed boundary traces.
template<class T, class Function>
auto map(const Field<T>& field, Function function) { return fieldUnary(field, function); }
inline ScalarField dot(const VectorField& a,const VectorField& b) {
    return fieldBinary(a,b,[](Vec3 x,Vec3 y){return babelsim::dot(x,y);});
}
inline VectorField cross(const VectorField& a,const VectorField& b) {
    return fieldBinary(a,b,[](Vec3 x,Vec3 y){return babelsim::cross(x,y);});
}
inline ScalarField max(const ScalarField& a,double bound) {
    ScalarField r(a.mesh(),a.location(),"math.max"); r.useCalculatedBoundary();
    r.evaluate(a,[bound](double v){return std::max(v,bound);}); return r;
}
inline ScalarField sqrt(const ScalarField& a) {
    ScalarField r(a.mesh(),a.location(),"math.sqrt"); r.useCalculatedBoundary();
    r.evaluate(a,[](double v){return std::sqrt(v);}); return r;
}
inline TensorField transpose(const TensorField& a) {
    TensorField r(a.mesh(),a.location(),"transpose"); r.useCalculatedBoundary();
    r.evaluate(a,[](const Tensor3& v){return babelsim::transpose(v);}); return r;
}
inline ScalarField trace(const TensorField& a) {
    ScalarField r(a.mesh(),a.location(),"trace"); r.useCalculatedBoundary();
    r.evaluate(a,[](const Tensor3& v){return babelsim::trace(v);}); return r;
}
inline TensorField isotropic(const ScalarField& a) {
    TensorField r(a.mesh(),a.location(),"isotropic"); r.useCalculatedBoundary();
    r.evaluate(a,[](double v){Tensor3 t{}; for(int i=0;i<3;++i)t[i][i]=v; return t;}); return r;
}
double sum(const ScalarField&);
double integral(const ScalarField&); // cell volumes, or face areas
double max(const ScalarField&);
double normL2(const ScalarField&); // unweighted Euclidean norm over owned entities
double normL2(const VectorField&);
} // namespace math
} // namespace babelsim
