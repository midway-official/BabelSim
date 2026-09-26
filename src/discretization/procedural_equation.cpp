#include "babelsim/equ.h"
#include "babelsim/math.h"
#include "babelsim/operators.h"
#include "internal/fvm_execution.h"
#include "internal/compute_backend.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace babelsim::equ {
namespace impl = babelsim::detail;
template<class T> struct Equation<T>::Storage {
    explicit Storage(Field<T>& x) : unknown(&x), coefficients(x.mesh()) {}
    Field<T>* unknown;
    DiscreteEquation<T> coefficients;
    EquationControl control;
    // 本方程已绑定的唯一面通量上下文。不同通量会给出不同的边界迹，
    // 一个方程只能用一个；记录首个指针以便不一致时立即报错。
    const ScalarField* boundary_flux = nullptr;
    struct DiffusionSnapshot {
        ScalarField coefficient;
        ScalarField boundary;
        ScalarField correction;
    };
    struct DiffusionContribution {
        std::shared_ptr<const DiffusionSnapshot> snapshot;
        double factor;
    };
    std::vector<DiffusionContribution> diffusion;

};
namespace detail {
struct EquationAccess {
    template<class T> static auto& get(Equation<T>& a) { return *a.storage_; }
    template<class T> static const auto& get(const Equation<T>& a) { return *a.storage_; }
};
}
using detail::EquationAccess;
namespace {
auto& backend() { return impl::execution().backend(); }
template<class T> void cell(const Field<T>& f, const Mesh& mesh) {
    f.validateStorage();
    if (&f.mesh() != &mesh || f.location() != FieldLocation::Cell)
        throw std::invalid_argument("equ requires compatible cell fields");
}
template<class T> auto& state(const Equation<T>& a) {
    auto& s = EquationAccess::get(a);
    if (&s.unknown->mesh() != &impl::execution().mesh())
        throw std::invalid_argument("equ matrix belongs to a different run mesh");
    return s;
}
template<class T> auto& state(Equation<T>& a) {
    state(static_cast<const Equation<T>&>(a));
    return EquationAccess::get(a);
}
void finite(double v) {
    if (!std::isfinite(v)) throw std::invalid_argument("nonfinite equ coefficient");
}
template<class T> void accumulate(DiscreteEquation<T>& a, const DiscreteEquation<T>& b, double f) {
    finite(f);
    for (std::size_t i=0; i<a.diagonal.size(); ++i) {
        a.diagonal[i] += f*b.diagonal[i]; a.source[i] += f*b.source[i];
    }
    for (std::size_t i=0; i<a.upper.size(); ++i) {
        a.upper[i] += f*b.upper[i]; a.lower[i] += f*b.lower[i];
    }
}
template<class T> void finish(Field<T>& f) { backend().synchronize(f); }
template<class T> Field<T> resultField(const Equation<T>& a, const char* name) {
    return Field<T>(state(a).unknown->mesh(), FieldLocation::Cell, name);
}
// Freeze the actual explicit flux, including any limiter evaluated at assembly.
// A frozen gradient alone is insufficient when the limiter also depends on x.
ScalarField diffusionCorrection(const ScalarField& k, const ScalarField& x,
    GradientMethod gradientMethod, DiffusionMethod diffusionMethod) {
    ScalarField correction(x.mesh(), FieldLocation::Face, "equ.frozenCorrection");
    if (diffusionMethod == DiffusionMethod::Orthogonal) return correction;
    VectorField g(x.mesh(), FieldLocation::Cell, "equ.frozenGradient");
    gradient(x, g, gradientMethod); finish(g);
    diffusionFlux(k, x, g, correction, diffusionMethod);
    ScalarField orthogonal(x.mesh(), FieldLocation::Face, "equ.orthogonalFlux");
    diffusionFlux(k, x, g, orthogonal, DiffusionMethod::Orthogonal);
    correction -= orthogonal;
    return correction;
}
}
template<class T> Equation<T>::Equation(Field<T>& x, const EquationControl& control)
    : storage_(std::make_unique<Storage>(x)) {
    cell(x, impl::execution().mesh());
    if (x.calculatedBoundary()) throw std::invalid_argument("equ unknown requires physical boundaries");
    control.linear.validate();
    control.requireSpatial(control.name);
    storage_->control = control;
    for (auto& entry : storage_->control.terms) {
        entry.second.requireComplete(control.name + ".term." + entry.first);
    }
    storage_->coefficients.numerical_identity = control.name;
}
template<class T> OperatorOptions Equation<T>::options(const std::string& term) const {
    return storage_->control.options(term);
}
template<class T> const std::string& Equation<T>::name() const { return storage_->control.name; }
template<class T> Equation<T>::~Equation() = default;
template<class T> Equation<T>::Equation(Equation&&) noexcept = default;
template<class T> Equation<T>& Equation<T>::operator=(Equation&&) noexcept = default;
template<class T> void clear(Equation<T>& a) {
    state(a).coefficients.reset(); state(a).diffusion.clear(); state(a).boundary_flux = nullptr;
}
template<class T> Equation<T> copy(const Equation<T>& a) {
    Equation<T> b(*state(a).unknown, state(a).control); state(b).coefficients = state(a).coefficients;
    state(b).control=state(a).control;
    state(b).diffusion=state(a).diffusion; state(b).boundary_flux=state(a).boundary_flux; return b;
}
template<class T> void div(Equation<T>& a, const ScalarField& phi, double factor, const std::string& term, const OperatorOptions& overrides) {
    finite(factor);
    auto& s = state(a); auto& x = *s.unknown;
    if (&phi.mesh()!=&x.mesh() || phi.location()!=FieldLocation::Face)
        throw std::invalid_argument("equ::div requires a face flux");
    if (s.boundary_flux != nullptr && s.boundary_flux != &phi)
        throw std::invalid_argument("one equation requires one boundary flux context");
    s.boundary_flux = &phi;
    backend().synchronize(const_cast<ScalarField&>(phi));
    x.setBoundaryFlux(phi); finish(x);
    auto m = a.options(term); m.overlay(overrides);
    addConvection(s.coefficients, phi, x, *m.convection,
        *m.interpolation, *m.gradient, factor);
}
template<class T> void laplacian(Equation<T>& a, double k, double factor, const std::string& term, const OperatorOptions& overrides) {
    finite(k); finite(factor);
    auto& s = state(a); auto& x = *s.unknown; finish(x);
    DiscreteEquation<T> contribution(x.mesh()); auto m = a.options(term); m.overlay(overrides);
    addDiffusion(contribution, k, x, *m.gradient, *m.diffusion);
    accumulate(s.coefficients, contribution, -factor);
    if constexpr(std::is_same_v<T,double>) {
        ScalarField kf(x.mesh(),FieldLocation::Face,"equ.frozenDiffusivity",k);
        auto correction = diffusionCorrection(kf, x, *m.gradient, *m.diffusion);
        using Snapshot=typename std::remove_reference_t<decltype(s)>::DiffusionSnapshot;
        s.diffusion.push_back({std::make_shared<Snapshot>(Snapshot{std::move(kf),x,std::move(correction)}),factor});
    }

}
template<class T> void laplacian(Equation<T>& a, const ScalarField& k, double factor, const std::string& term, const OperatorOptions& overrides) {
    finite(factor); auto& s = state(a); auto& x = *s.unknown;
    if (&k.mesh()!=&x.mesh() || (k.location()!=FieldLocation::Cell && k.location()!=FieldLocation::Face))
        throw std::invalid_argument("equ diffusivity must be a compatible cell or face field");
    backend().synchronize(const_cast<ScalarField&>(k)); finish(x);
    ScalarField kf(x.mesh(), FieldLocation::Face, "equ.diffusivity");
    auto m = a.options(term); m.overlay(overrides); const ScalarField* face = &k;
    if (k.location()==FieldLocation::Cell) {
        interpolate(k, kf, m.coefficientInterpolation.value_or(*m.interpolation),
            m.coefficientGradient.value_or(*m.gradient));
        finish(kf); face=&kf;
    }
    DiscreteEquation<T> contribution(x.mesh());
    addDiffusion(contribution, *face, x, *m.gradient, *m.diffusion);
    accumulate(s.coefficients, contribution, -factor);
    if constexpr(std::is_same_v<T,double>) {
        auto correction = diffusionCorrection(*face, x, *m.gradient, *m.diffusion);
        using Snapshot=typename std::remove_reference_t<decltype(s)>::DiffusionSnapshot;
        s.diffusion.push_back({std::make_shared<Snapshot>(Snapshot{*face,x,std::move(correction)}),factor});
    }

}
template<class T> void reaction(Equation<T>& a, double c) {
    finite(c); auto& s=state(a); const auto& mesh=s.unknown->mesh();
    for (Index i: impl::meshData(mesh).owned_cells) s.coefficients.diagonal[i]+=c*mesh.cellVolume(i);
}
template<class T> void reaction(Equation<T>& a, const ScalarField& c, double factor) {
    finite(factor); auto& s=state(a); const auto& mesh=s.unknown->mesh(); cell(c,mesh);
    for (Index i: impl::meshData(mesh).owned_cells) {
        finite(impl::fieldData(c)[i]);
        s.coefficients.diagonal[i]+=factor*impl::fieldData(c)[i]*mesh.cellVolume(i);
    }
}
template<class T> void source(Equation<T>& a, T value) {
    auto& s=state(a); const auto& mesh=s.unknown->mesh();
    for (Index i: impl::meshData(mesh).owned_cells) s.coefficients.source[i]+=mesh.cellVolume(i)*value;
}
template<class T> void source(Equation<T>& a, const Field<T>& value, double factor) {
    finite(factor); auto& s=state(a); const auto& mesh=s.unknown->mesh(); cell(value,mesh);
    for (Index i: impl::meshData(mesh).owned_cells)
        s.coefficients.source[i]+=factor*mesh.cellVolume(i)*impl::fieldData(value)[i];
}
namespace {
template<class T, class Capacity> void timeTerm(Equation<T>& a, const Capacity& capacity,
    const Field<T>& previous, double dt, TimeMethod method, const Field<T>* older, double olddt) {
    auto& s=state(a); const auto& mesh=s.unknown->mesh(); cell(previous,mesh);
    if (method==TimeMethod::Steady) return;
    if (!(dt>0) || !std::isfinite(dt)) throw std::invalid_argument("equ::ddt requires positive finite dt");
    if constexpr (!std::is_same_v<Capacity,double>) cell(capacity,mesh);
    double a0=1/dt, a1=1/dt, a2=0;
    if (method==TimeMethod::BDF2) {
        if (!older) throw std::invalid_argument("BDF2 requires explicit older field; use Euler to start");
        cell(*older,mesh);
        if (olddt==0) olddt=dt;
        if (!(olddt>0) || !std::isfinite(olddt)) throw std::invalid_argument("invalid previous dt");
        const double r=dt/olddt;
        a0=(1+2*r)/((1+r)*dt); a1=(1+r)/dt; a2=-r*r/((1+r)*dt);
    }
    for (Index i: impl::meshData(mesh).owned_cells) {
        double c;
        if constexpr (std::is_same_v<Capacity,double>) c=capacity;
        else c=impl::fieldData(capacity)[i];
        if (!(c>0) || !std::isfinite(c)) throw std::invalid_argument("time capacity must be positive and finite");
        const double w=c*mesh.cellVolume(i);
        s.coefficients.diagonal[i]+=w*a0;
        s.coefficients.source[i]+=w*a1*impl::fieldData(previous)[i];
        if (a2!=0) s.coefficients.source[i]+=w*a2*impl::fieldData(*older)[i];
    }
}
}
template<class T> void ddt(Equation<T>& a,double c,const Field<T>& old,double dt,TimeMethod m,const Field<T>* older,double olddt) {
    timeTerm(a,c,old,dt,m,older,olddt);
}
template<class T> void ddt(Equation<T>& a,const ScalarField& c,const Field<T>& old,double dt,TimeMethod m,const Field<T>* older,double olddt) {
    timeTerm(a,c,old,dt,m,older,olddt);
}
template<class T,class Capacity> void historyTime(Equation<T>& a,const Capacity& c,const time::History<T>& h) {
    if(numericalMethods().time==TimeMethod::Steady) return;
    if(h.levels()==0) throw std::logic_error("saveOld must precede a transient assembly");
    const auto configured=numericalMethods().time;
    const auto method=configured==TimeMethod::BDF2 && h.levels()<2 ? TimeMethod::Euler : configured;
    ddt(a,c,h.previous(),h.dt(),method,&h.older(),h.previousDt());
}
template<class T> void ddt(Equation<T>& a,double c,const time::History<T>& h) {historyTime(a,c,h);}
template<class T> void ddt(Equation<T>& a,const ScalarField& c,const time::History<T>& h) {historyTime(a,c,h);}
template<class T> void scale(Equation<T>& a,double f) {
    finite(f); auto& q=state(a).coefficients;
    for (auto& v:q.diagonal) v*=f;
    for (auto& v:q.upper) v*=f;
    for (auto& v:q.lower) v*=f;
    for (auto& v:q.source) v=f*v;
    for (auto& contribution:state(a).diffusion) contribution.factor*=f;
}
template<class T> void add(Equation<T>& a,const Equation<T>& b,double f) {
    auto& sa=state(a); const auto& sb=state(b);
    if (sa.unknown != sb.unknown) throw std::invalid_argument("matrix addition requires same unknown binding");
    if (sb.boundary_flux != nullptr) {
        if (sa.boundary_flux != nullptr && sa.boundary_flux != sb.boundary_flux)
            throw std::invalid_argument("one equation requires one boundary flux context");
        sa.boundary_flux = sb.boundary_flux;
    }
    // Copy first so self-addition does not invalidate iterators.
    auto contributions=sb.diffusion;
    accumulate(sa.coefficients,sb.coefficients,f);
    for(auto& contribution:contributions) {
        contribution.factor*=f; sa.diffusion.push_back(std::move(contribution));
    }
}
template<class T> void addDiagonal(Equation<T>& a,const ScalarField& d) {
    auto& s=state(a); cell(d,s.unknown->mesh());
    for (Index i:impl::meshData(d.mesh()).owned_cells) s.coefficients.diagonal[i]+=impl::fieldData(d)[i];
}
template<class T> void addRhs(Equation<T>& a,const Field<T>& b) {
    auto& s=state(a); cell(b,s.unknown->mesh());
    for (Index i:impl::meshData(b.mesh()).owned_cells) s.coefficients.source[i]+=impl::fieldData(b)[i];
}
template<class T> ScalarField diagonal(const Equation<T>& a) {
    const auto& s=state(a); ScalarField d(s.unknown->mesh(),FieldLocation::Cell,"equ.diagonal");
    for (Index i:impl::meshData(d.mesh()).owned_cells) impl::fieldData(d)[i]=s.coefficients.diagonal[i];
    finish(d); return d;
}
template<class T> Field<T> rhs(const Equation<T>& a) {
    auto b=resultField(a,"equ.rhs"); const auto& q=state(a).coefficients;
    for (Index i:impl::meshData(b.mesh()).owned_cells) impl::fieldData(b)[i]=q.source[i];
    finish(b); return b;
}
template<class T> Field<T> apply(const Equation<T>& a,const Field<T>& x) {
    const auto& s=state(a); const auto& mesh=s.unknown->mesh(); cell(x,mesh);
    backend().synchronize(const_cast<Field<T>&>(x)); auto y=resultField(a,"equ.apply");
    auto* values=impl::fieldData(y); const auto* xv=impl::fieldData(x); const auto& q=s.coefficients;
    for (Index i:impl::meshData(mesh).owned_cells) values[i]=q.diagonal[i]*xv[i];
    for (Index f=0;f<mesh.faceCount();++f) {
        const Index o=mesh.owner(f), n=mesh.neighbour(f); if(n==invalid_index) continue;
        if(impl::isOwned(mesh,o)) values[o]+=q.upper[f]*xv[n];
        if(impl::isOwned(mesh,n)) values[n]+=q.lower[f]*xv[o];
    }
    finish(y); return y;
}
template<class T> Field<T> residual(const Equation<T>& a,const Field<T>& x) {
    auto r=apply(a,x); const auto& q=state(a).coefficients;
    for(Index i:impl::meshData(x.mesh()).owned_cells) impl::fieldData(r)[i]=q.source[i]-impl::fieldData(r)[i];
    finish(r); return r;
}
template<class T> void relax(Equation<T>& a,const Field<T>& previous,double alpha) {
    if (!(alpha>0 && alpha<=1)) throw std::invalid_argument("relaxation must be in (0,1]");
    auto& s=state(a); cell(previous,s.unknown->mesh()); auto& q=s.coefficients;
    for(Index i:impl::meshData(previous.mesh()).owned_cells) {
        const double extra=q.diagonal[i]*(1/alpha-1);
        q.diagonal[i]+=extra; q.source[i]+=extra*impl::fieldData(previous)[i];
    }
}
void reference(Equation<double>& a,Index global_cell,double value) {
    finite(value); auto& s=state(a); auto& q=s.coefficients; const auto& mesh=s.unknown->mesh();
    bool found=false;
    for(Index i:impl::meshData(mesh).owned_cells) {
        if(impl::globalCellId(mesh,i)!=global_cell) continue;
        if (!(q.diagonal[i]>0)) throw std::invalid_argument("reference requires a positive diagonal");
        q.source[i]+=q.diagonal[i]*value; q.diagonal[i]*=2; found=true;
    }
    if(backend().all(!found)) throw std::invalid_argument("reference global cell does not exist");
}
void reference(Equation<double>& a,double value) {
    finite(value);
    const auto& s=state(a); const auto& mesh=s.unknown->mesh(); const auto& q=s.coefficients;
    auto rowSum=q.diagonal;
    auto magnitude=q.diagonal;
    for(auto& v:magnitude) v=std::abs(v);
    for(Index f=0;f<mesh.faceCount();++f) {
        const Index o=mesh.owner(f),n=mesh.neighbour(f);
        if(n==invalid_index) continue;
        if(impl::isOwned(mesh,o)) {rowSum[o]+=q.upper[f];magnitude[o]+=std::abs(q.upper[f]);}
        if(impl::isOwned(mesh,n)) {rowSum[n]+=q.lower[f];magnitude[n]+=std::abs(q.lower[f]);}
    }
    bool constantMode=true;
    for(Index i:impl::meshData(mesh).owned_cells)
        constantMode=constantMode && std::isfinite(rowSum[i]) &&
            std::abs(rowSum[i])<=64*std::numeric_limits<double>::epsilon()*magnitude[i];
    if(backend().all(constantMode)) reference(a,0,value);
}

ScalarField faceFlux(const Equation<double>& a,const ScalarField& solution) {
    const auto& s=state(a); cell(solution,s.unknown->mesh());
    ScalarField result(solution.mesh(),FieldLocation::Face,"equ.diffusionFlux");
    // Orthogonal flux ignores this gradient; the explicit part is already frozen.
    VectorField zeroGradient(solution.mesh(),FieldLocation::Cell);
    for(const auto& contribution:s.diffusion) {
        const auto& snapshot=*contribution.snapshot;
        auto value=snapshot.boundary;
        value=solution; // keep the assembled boundary constraints
        finish(value);
        ScalarField increment(solution.mesh(),FieldLocation::Face,"equ.diffusionFluxTerm");
        diffusionFlux(snapshot.coefficient,value,zeroGradient,increment,
            DiffusionMethod::Orthogonal);
        increment += snapshot.correction;
        result.addScaled(contribution.factor,increment);
    }
    finish(result); return result;
}
template<class T> SolveResult solveConfigured(const Equation<T>& a,Field<T>& x,const LinearSolverConfig* config) {
    const auto& s=state(a);
    if(s.unknown!=&x) throw std::invalid_argument("solve target differs from matrix unknown");
    if constexpr(std::is_same_v<T,double>) return backend().solve(s.coefficients,x,config ? *config : s.control.linear);
    else {
        const auto components=backend().solve(s.coefficients,x,config ? *config : s.control.linear);
        SolveResult result; result.status=SolveStatus::Converged;
        for(const auto& c:components) {
            if(!c.healthy()) result.status=SolveStatus::NumericalFailure;
            else if(!c.converged() && result.status!=SolveStatus::NumericalFailure) result.status=SolveStatus::MaxIterations;
            result.iterations+=c.iterations;
            result.initial_residual=std::hypot(result.initial_residual,c.initial_residual);
            result.final_residual=std::hypot(result.final_residual,c.final_residual);
            result.relative_residual=std::max(result.relative_residual,c.relative_residual);
            result.performance+=c.performance;
        }
        return result;
    }
}
template<class T> SolveResult solve(const Equation<T>& a,Field<T>& x) { return solveConfigured(a,x,nullptr); }
template<class T> SolveResult solve(const Equation<T>& a,Field<T>& x,const LinearSolverConfig& config) {
    config.validate(); return solveConfigured(a,x,&config);
}
template<class T> SolveResult solve(const Equation<T>& a) {
    return solveConfigured(a, *EquationAccess::get(a).unknown, nullptr);
}
template<class T> SolveResult solve(const Equation<T>& a, const LinearSolverConfig& config) {
    config.validate();
    return solveConfigured(a, *EquationAccess::get(a).unknown, &config);
}
template<class T> void Equation<T>::reset() { clear(*this); }
template<class T> Equation<T> Equation<T>::copy() const { return equ::copy(*this); }
template<class T> ScalarField Equation<T>::diagonal() const { return equ::diagonal(*this); }
template<class T> Field<T> Equation<T>::rhs() const { return equ::rhs(*this); }
template<class T> void Equation<T>::reference(Index cellIndex, double value) {
    if constexpr (std::is_same_v<T, double>) equ::reference(*this, cellIndex, value);
    else throw std::invalid_argument("reference requires a scalar equation");
}
template<class T> void Equation<T>::referenceIfUnanchored(double value) {
    if constexpr (std::is_same_v<T, double>) equ::reference(*this, value);
    else throw std::invalid_argument("reference requires a scalar equation");
}
#define INSTANTIATE(T) \
 template class Equation<T>; \
 template void clear(Equation<T>&); \
 template Equation<T> copy(const Equation<T>&); \
 template void div(Equation<T>&,const ScalarField&,double,const std::string&,const OperatorOptions&); \
 template void laplacian(Equation<T>&,double,double,const std::string&,const OperatorOptions&); \
 template void laplacian(Equation<T>&,const ScalarField&,double,const std::string&,const OperatorOptions&); \
 template void reaction(Equation<T>&,double); \
 template void reaction(Equation<T>&,const ScalarField&,double); \
 template void source(Equation<T>&,T); \
 template void source(Equation<T>&,const Field<T>&,double); \
 template void ddt(Equation<T>&,double,const Field<T>&,double,TimeMethod,const Field<T>*,double); \
 template void ddt(Equation<T>&,const ScalarField&,const Field<T>&,double,TimeMethod,const Field<T>*,double); \
 template void ddt(Equation<T>&,double,const time::History<T>&); \
 template void ddt(Equation<T>&,const ScalarField&,const time::History<T>&); \
 template void scale(Equation<T>&,double); \
 template void add(Equation<T>&,const Equation<T>&,double); \
 template void addDiagonal(Equation<T>&,const ScalarField&); \
 template void addRhs(Equation<T>&,const Field<T>&); \
 template ScalarField diagonal(const Equation<T>&); \
 template Field<T> rhs(const Equation<T>&); \
 template Field<T> apply(const Equation<T>&,const Field<T>&); \
 template Field<T> residual(const Equation<T>&,const Field<T>&); \
 template void relax(Equation<T>&,const Field<T>&,double); \
 template SolveResult solve(const Equation<T>&,Field<T>&); \
 template SolveResult solve(const Equation<T>&,Field<T>&,const LinearSolverConfig&); \
 template SolveResult solve(const Equation<T>&); \
 template SolveResult solve(const Equation<T>&,const LinearSolverConfig&);
INSTANTIATE(double)
INSTANTIATE(Vec3)
#undef INSTANTIATE
} // namespace babelsim::equ

namespace babelsim::diagnostics {
template<class T> double relativeResidual(const equ::Equation<T>& equation, const Field<T>& x) {
    auto ax = equ::apply(equation, x); auto b = equation.rhs(); auto r = b - ax;
    return math::normL2(r) / std::max(math::normL2(ax) + math::normL2(b), 1e-30);
}
template double relativeResidual(const equ::Equation<double>&, const ScalarField&);
template double relativeResidual(const equ::Equation<Vec3>&, const VectorField&);
}
