#include "babelsim/math.h"
#include "internal/fvm_execution.h"
#include "internal/compute_backend.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"

namespace babelsim::math {
namespace {
template<class T, class Function>
double reduceField(const Field<T>& f, Function fn, bool maximum, bool weighted) {
    f.validateStorage();
    if (&f.mesh() != &detail::execution().mesh() || f.location() == FieldLocation::Vertex)
        throw std::invalid_argument("math reduction requires a cell/face field on the run mesh");
    const auto& ids = f.location() == FieldLocation::Cell
        ? detail::meshData(f.mesh()).owned_cells : detail::meshData(f.mesh()).owned_faces;
    double local = maximum ? -std::numeric_limits<double>::infinity() : 0.0;
    for (Index i : ids) {
        const double weight = !weighted ? 1.0 : f.location() == FieldLocation::Cell
            ? f.mesh().cellVolume(i) : f.mesh().faceArea(i);
        const double value = weight * fn(detail::fieldData(f)[i]);
        if (maximum) local = std::max(local, value); else local += value;
    }
    double global;
    auto& backend = detail::execution().backend();
    if (maximum) backend.maximum(&local, &global, 1); else backend.sum(&local, &global, 1);
    return global;
}
}
double sum(const ScalarField& f) { return reduceField(f, [](double v){return v;}, false, false); }
double integral(const ScalarField& f) { return reduceField(f, [](double v){return v;}, false, true); }
double max(const ScalarField& f) { return reduceField(f, [](double v){return v;}, true, false); }
double normL2(const ScalarField& f) {
    return std::sqrt(reduceField(f, [](double v){return v*v;}, false, false));
}
double normL2(const VectorField& f) {
    return std::sqrt(reduceField(f, [](Vec3 v){return squaredNorm(v);}, false, false));
}

void add(const ScalarField& increment, ScalarField& target, FaceRegion region) {
    if (&increment.mesh()!=&target.mesh() || target.location()!=FieldLocation::Face ||
        increment.location()!=FieldLocation::Face || &target.mesh()!=&detail::execution().mesh())
        throw std::invalid_argument("face increment requires matching run mesh face fields");
    if (region!=FaceRegion::All && region!=FaceRegion::Interior)
        throw std::invalid_argument("invalid face region");
    auto& backend=detail::execution().backend();
    backend.synchronize(const_cast<ScalarField&>(increment));
    for (Index face:detail::meshData(target.mesh()).owned_faces) {
        if(region==FaceRegion::Interior && target.mesh().boundaryFace(face)) continue;
        detail::fieldData(target)[face]+=detail::fieldData(increment)[face];
    }
    backend.synchronize(target);
}
void subtract(const ScalarField& increment,ScalarField& target,FaceRegion region) {
    const auto negative=-increment; add(negative,target,region);
}
void subtract(const ScalarField& coefficient,const VectorField& gradient,VectorField& target) {
    target.addProduct(-1.0,coefficient,gradient);
    detail::execution().backend().synchronize(target);
}


} // namespace babelsim::math
