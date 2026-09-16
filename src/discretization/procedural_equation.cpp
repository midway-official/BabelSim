#include "babelsim/equ.h"
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
template<class T> struct Matrix<T>::Storage {
    explicit Storage(Field<T>& x) : unknown(&x), coefficients(x.mesh()) {}
    Field<T>* unknown;
    DiscreteEquation<T> coefficients;
    struct DiffusionSnapshot {
        ScalarField coefficient;
        ScalarField boundary;
        VectorField gradient;
    };
    struct DiffusionContribution {
        std::shared_ptr<const DiffusionSnapshot> snapshot;
        double factor;
        DiffusionMethod method;
    };
    std::vector<DiffusionContribution> diffusion;

};
namespace detail {
struct MatrixAccess {
    template<class T> static auto& get(Matrix<T>& a) { return *a.storage_; }
    template<class T> static const auto& get(const Matrix<T>& a) { return *a.storage_; }
};
}
using detail::MatrixAccess;
namespace {
auto& backend() { return impl::execution().backend(); }
template<class T> void cell(const Field<T>& f, const Mesh& mesh) {
    f.validateStorage();
    if (&f.mesh() != &mesh || f.location() != FieldLocation::Cell)
        throw std::invalid_argument("equ requires compatible cell fields");
}
template<class T> auto& state(const Matrix<T>& a) {
    auto& s = MatrixAccess::get(a);
    if (&s.unknown->mesh() != &impl::execution().mesh())
        throw std::invalid_argument("equ matrix belongs to a different run mesh");
    return s;
}
template<class T> auto& state(Matrix<T>& a) {
    state(static_cast<const Matrix<T>&>(a));
    return MatrixAccess::get(a);
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
template<class T> Field<T> resultField(const Matrix<T>& a, const char* name) {
    return Field<T>(state(a).unknown->mesh(), FieldLocation::Cell, name);
}
}
template<class T> Matrix<T>::Matrix(Field<T>& x) : storage_(std::make_unique<Storage>(x)) {
    cell(x, impl::execution().mesh());
    if (x.calculatedBoundary()) throw std::invalid_argument("equ unknown requires physical boundaries");
}
template<class T> Matrix<T>::~Matrix() = default;
template<class T> Matrix<T>::Matrix(Matrix&&) noexcept = default;
template<class T> Matrix<T>& Matrix<T>::operator=(Matrix&&) noexcept = default;
template<class T> void clear(Matrix<T>& a) { state(a).coefficients.reset(); state(a).diffusion.clear(); }
template<class T> Matrix<T> copy(const Matrix<T>& a) {
    Matrix<T> b(*state(a).unknown); state(b).coefficients = state(a).coefficients; state(b).diffusion=state(a).diffusion; return b;
}
template<class T> void div(Matrix<T>& a, const ScalarField& phi, double factor) {
    finite(factor);
    auto& s = state(a); auto& x = *s.unknown;
    if (&phi.mesh()!=&x.mesh() || phi.location()!=FieldLocation::Face)
        throw std::invalid_argument("equ::div requires a face flux");
    backend().synchronize(const_cast<ScalarField&>(phi));
    x.setBoundaryFlux(phi); finish(x);
    const auto& m = numericalMethods();
    addConvection(s.coefficients, phi, x, m.convectionFor(x.name()),
        m.interpolationFor(x.name()), m.gradientFor(x.name()), factor);
}
template<class T> void laplacian(Matrix<T>& a, double k, double factor) {
    finite(k); finite(factor);
    auto& s = state(a); auto& x = *s.unknown; finish(x);
    DiscreteEquation<T> contribution(x.mesh()); const auto& m = numericalMethods();
    addDiffusion(contribution, k, x, m.gradientFor(x.name()), m.diffusionFor(x.name()));
    accumulate(s.coefficients, contribution, -factor);
    if constexpr(std::is_same_v<T,double>) {
        ScalarField kf(x.mesh(),FieldLocation::Face,"equ.frozenDiffusivity",k);
        VectorField g(x.mesh(),FieldLocation::Cell,"equ.frozenGradient");
        gradient(x,g,m.gradientFor(x.name())); finish(g);
        using Snapshot=typename std::remove_reference_t<decltype(s)>::DiffusionSnapshot;
        s.diffusion.push_back({std::make_shared<Snapshot>(Snapshot{std::move(kf),x,std::move(g)}),factor,m.diffusionFor(x.name())});
    }

}
template<class T> void laplacian(Matrix<T>& a, const ScalarField& k, double factor) {
    finite(factor); auto& s = state(a); auto& x = *s.unknown;
    if (&k.mesh()!=&x.mesh() || (k.location()!=FieldLocation::Cell && k.location()!=FieldLocation::Face))
        throw std::invalid_argument("equ diffusivity must be a compatible cell or face field");
    backend().synchronize(const_cast<ScalarField&>(k)); finish(x);
    ScalarField kf(x.mesh(), FieldLocation::Face, "equ.diffusivity");
    const auto& m = numericalMethods(); const ScalarField* face = &k;
    if (k.location()==FieldLocation::Cell) {
        interpolate(k, kf, m.interpolationFor(x.name()), m.gradientFor(x.name()));
        finish(kf); face=&kf;
    }
    DiscreteEquation<T> contribution(x.mesh());
    addDiffusion(contribution, *face, x, m.gradientFor(x.name()), m.diffusionFor(x.name()));
    accumulate(s.coefficients, contribution, -factor);
    if constexpr(std::is_same_v<T,double>) {
        VectorField g(x.mesh(),FieldLocation::Cell,"equ.frozenGradient");
        gradient(x,g,m.gradientFor(x.name())); finish(g);
        using Snapshot=typename std::remove_reference_t<decltype(s)>::DiffusionSnapshot;
        s.diffusion.push_back({std::make_shared<Snapshot>(Snapshot{*face,x,std::move(g)}),factor,m.diffusionFor(x.name())});
    }

}
template<class T> void reaction(Matrix<T>& a, double c) {
    finite(c); auto& s=state(a); const auto& mesh=s.unknown->mesh();
    for (Index i: impl::meshData(mesh).owned_cells) s.coefficients.diagonal[i]+=c*mesh.cellVolume(i);
}
template<class T> void reaction(Matrix<T>& a, const ScalarField& c, double factor) {
    finite(factor); auto& s=state(a); const auto& mesh=s.unknown->mesh(); cell(c,mesh);
    for (Index i: impl::meshData(mesh).owned_cells) {
        finite(impl::fieldData(c)[i]);
        s.coefficients.diagonal[i]+=factor*impl::fieldData(c)[i]*mesh.cellVolume(i);
    }
}
template<class T> void source(Matrix<T>& a, T value) {
    auto& s=state(a); const auto& mesh=s.unknown->mesh();
    for (Index i: impl::meshData(mesh).owned_cells) s.coefficients.source[i]+=mesh.cellVolume(i)*value;
}
template<class T> void source(Matrix<T>& a, const Field<T>& value, double factor) {
    finite(factor); auto& s=state(a); const auto& mesh=s.unknown->mesh(); cell(value,mesh);
    for (Index i: impl::meshData(mesh).owned_cells)
        s.coefficients.source[i]+=factor*mesh.cellVolume(i)*impl::fieldData(value)[i];
}
namespace {
template<class T, class Capacity> void timeTerm(Matrix<T>& a, const Capacity& capacity,
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
template<class T> void ddt(Matrix<T>& a,double c,const Field<T>& old,double dt,TimeMethod m,const Field<T>* older,double olddt) {
    timeTerm(a,c,old,dt,m,older,olddt);
}
template<class T> void ddt(Matrix<T>& a,const ScalarField& c,const Field<T>& old,double dt,TimeMethod m,const Field<T>* older,double olddt) {
    timeTerm(a,c,old,dt,m,older,olddt);
}
template<class T,class Capacity> void historyTime(Matrix<T>& a,const Capacity& c,const math::History<T>& h) {
    if(numericalMethods().time==TimeMethod::Steady) return;
    if(h.levels()==0) throw std::logic_error("saveOld must precede a transient assembly");
    const auto configured=numericalMethods().time;
    const auto method=configured==TimeMethod::BDF2 && h.levels()<2 ? TimeMethod::Euler : configured;
    ddt(a,c,h.previous(),h.dt(),method,&h.older(),h.previousDt());
}
template<class T> void ddt(Matrix<T>& a,double c,const math::History<T>& h) {historyTime(a,c,h);}
template<class T> void ddt(Matrix<T>& a,const ScalarField& c,const math::History<T>& h) {historyTime(a,c,h);}
template<class T> void scale(Matrix<T>& a,double f) {
    finite(f); auto& q=state(a).coefficients;
    for (auto& v:q.diagonal) v*=f;
    for (auto& v:q.upper) v*=f;
    for (auto& v:q.lower) v*=f;
    for (auto& v:q.source) v=f*v;
    for (auto& contribution:state(a).diffusion) contribution.factor*=f;
}
template<class T> void add(Matrix<T>& a,const Matrix<T>& b,double f) {
    auto& sa=state(a); const auto& sb=state(b);
    if (sa.unknown != sb.unknown) throw std::invalid_argument("matrix addition requires same unknown binding");
    // Copy first so self-addition does not invalidate iterators.
    auto contributions=sb.diffusion;
    accumulate(sa.coefficients,sb.coefficients,f);
    for(auto& contribution:contributions) {
        contribution.factor*=f; sa.diffusion.push_back(std::move(contribution));
    }
}
template<class T> void addDiagonal(Matrix<T>& a,const ScalarField& d) {
    auto& s=state(a); cell(d,s.unknown->mesh());
    for (Index i:impl::meshData(d.mesh()).owned_cells) s.coefficients.diagonal[i]+=impl::fieldData(d)[i];
}
template<class T> void addRhs(Matrix<T>& a,const Field<T>& b) {
    auto& s=state(a); cell(b,s.unknown->mesh());
    for (Index i:impl::meshData(b.mesh()).owned_cells) s.coefficients.source[i]+=impl::fieldData(b)[i];
}
template<class T> ScalarField diagonal(const Matrix<T>& a) {
    const auto& s=state(a); ScalarField d(s.unknown->mesh(),FieldLocation::Cell,"equ.diagonal");
    for (Index i:impl::meshData(d.mesh()).owned_cells) impl::fieldData(d)[i]=s.coefficients.diagonal[i];
    finish(d); return d;
}
template<class T> ScalarField response(const Matrix<T>& a) {
    auto d=diagonal(a);
    for (Index i:impl::meshData(d.mesh()).owned_cells) {
        const double ap=impl::fieldData(d)[i];
        if (!(ap>0) || !std::isfinite(ap)) throw std::runtime_error("response requires positive finite diagonal");
        impl::fieldData(d)[i]=d.mesh().cellVolume(i)/ap;
    }
    finish(d); return d;
}
template<class T> Field<T> rhs(const Matrix<T>& a) {
    auto b=resultField(a,"equ.rhs"); const auto& q=state(a).coefficients;
    for (Index i:impl::meshData(b.mesh()).owned_cells) impl::fieldData(b)[i]=q.source[i];
    finish(b); return b;
}
template<class T> Field<T> apply(const Matrix<T>& a,const Field<T>& x) {
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
template<class T> Field<T> residual(const Matrix<T>& a,const Field<T>& x) {
    auto r=apply(a,x); const auto& q=state(a).coefficients;
    for(Index i:impl::meshData(x.mesh()).owned_cells) impl::fieldData(r)[i]=q.source[i]-impl::fieldData(r)[i];
    finish(r); return r;
}
template<class T> void relax(Matrix<T>& a,const Field<T>& previous,double alpha) {
    if (!(alpha>0 && alpha<=1)) throw std::invalid_argument("relaxation must be in (0,1]");
    auto& s=state(a); cell(previous,s.unknown->mesh()); auto& q=s.coefficients;
    for(Index i:impl::meshData(previous.mesh()).owned_cells) {
        const double extra=q.diagonal[i]*(1/alpha-1);
        q.diagonal[i]+=extra; q.source[i]+=extra*impl::fieldData(previous)[i];
    }
}
void reference(Matrix<double>& a,Index global_cell,double value) {
    finite(value); auto& s=state(a); auto& q=s.coefficients; const auto& mesh=s.unknown->mesh();
    bool found=false;
    for(Index i:impl::meshData(mesh).owned_cells) {
        if(impl::globalCellId(mesh,i)!=global_cell) continue;
        if (!(q.diagonal[i]>0)) throw std::invalid_argument("reference requires a positive diagonal");
        q.source[i]+=q.diagonal[i]*value; q.diagonal[i]*=2; found=true;
    }
    if(backend().all(!found)) throw std::invalid_argument("reference global cell does not exist");
}
void reference(Matrix<double>& a,double value) {
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

ScalarField faceFlux(const Matrix<double>& a,const ScalarField& solution) {
    const auto& s=state(a); cell(solution,s.unknown->mesh());
    ScalarField result(solution.mesh(),FieldLocation::Face,"equ.diffusionFlux");
    for(const auto& contribution:s.diffusion) {
        const auto& snapshot=*contribution.snapshot;
        auto value=snapshot.boundary;
        value=solution; // keep the assembled boundary constraints
        finish(value);
        ScalarField increment(solution.mesh(),FieldLocation::Face,"equ.diffusionFluxTerm");
        diffusionFlux(snapshot.coefficient,value,snapshot.gradient,increment,
            contribution.method);
        result.addScaled(contribution.factor,increment);
    }
    finish(result); return result;
}
template<class T> double relativeResidual(const Matrix<T>& a,const Field<T>& x) {
    auto ax=apply(a,x); auto b=rhs(a); auto r=b-ax;
    return math::normL2(r)/std::max(math::normL2(ax)+math::normL2(b),1e-30);
}
template<class T> SolveResult solveConfigured(const Matrix<T>& a,Field<T>& x,const LinearSolverConfig* config) {
    const auto& s=state(a);
    if(s.unknown!=&x) throw std::invalid_argument("solve target differs from matrix unknown");
    if constexpr(std::is_same_v<T,double>) return config ? backend().solve(s.coefficients,x,*config) : backend().solve(s.coefficients,x);
    else {
        const auto components=config ? backend().solve(s.coefficients,x,*config) : backend().solve(s.coefficients,x);
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
template<class T> SolveResult solve(const Matrix<T>& a,Field<T>& x) { return solveConfigured(a,x,nullptr); }
template<class T> SolveResult solve(const Matrix<T>& a,Field<T>& x,const LinearSolverConfig& config) {
    config.validate(); return solveConfigured(a,x,&config);
}
#define INSTANTIATE(T) \
 template class Matrix<T>; \
 template void clear(Matrix<T>&); \
 template Matrix<T> copy(const Matrix<T>&); \
 template void div(Matrix<T>&,const ScalarField&,double); \
 template void laplacian(Matrix<T>&,double,double); \
 template void laplacian(Matrix<T>&,const ScalarField&,double); \
 template void reaction(Matrix<T>&,double); \
 template void reaction(Matrix<T>&,const ScalarField&,double); \
 template void source(Matrix<T>&,T); \
 template void source(Matrix<T>&,const Field<T>&,double); \
 template void ddt(Matrix<T>&,double,const Field<T>&,double,TimeMethod,const Field<T>*,double); \
 template void ddt(Matrix<T>&,const ScalarField&,const Field<T>&,double,TimeMethod,const Field<T>*,double); \
 template void ddt(Matrix<T>&,double,const math::History<T>&); \
 template void ddt(Matrix<T>&,const ScalarField&,const math::History<T>&); \
 template void scale(Matrix<T>&,double); \
 template void add(Matrix<T>&,const Matrix<T>&,double); \
 template void addDiagonal(Matrix<T>&,const ScalarField&); \
 template void addRhs(Matrix<T>&,const Field<T>&); \
 template ScalarField diagonal(const Matrix<T>&); \
 template ScalarField response(const Matrix<T>&); \
 template Field<T> rhs(const Matrix<T>&); \
 template Field<T> apply(const Matrix<T>&,const Field<T>&); \
 template Field<T> residual(const Matrix<T>&,const Field<T>&); \
 template void relax(Matrix<T>&,const Field<T>&,double); \
 template double relativeResidual(const Matrix<T>&,const Field<T>&); \
 template SolveResult solve(const Matrix<T>&,Field<T>&); \
 template SolveResult solve(const Matrix<T>&,Field<T>&,const LinearSolverConfig&);
INSTANTIATE(double)
INSTANTIATE(Vec3)
#undef INSTANTIATE
} // namespace babelsim::equ
