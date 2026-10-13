#include "babelsim/math.h"
#include "internal/fvm_execution.h"
#include "internal/compute_backend.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"

namespace babelsim::math {
namespace {
enum class Reduction { Sum, Maximum };
enum class Weight { Unit, Measure };

// Global statistics count each entity once. Ghost cells and replicated faces
// must never enter the local reduction, even when their values are synchronized.
template <class T, class Function>
double reduceField(const Field<T>& field, Function valueOf, Reduction reduction, Weight weighting) {
    field.validateStorage();
    if (&field.mesh() != &detail::execution().mesh() || field.location() == FieldLocation::Vertex)
        throw std::invalid_argument("math reduction requires a cell/face field on the run mesh");
    const auto& owned_entities = field.location() == FieldLocation::Cell
                                     ? detail::meshData(field.mesh()).owned_cells
                                     : detail::meshData(field.mesh()).owned_faces;
    double local = reduction == Reduction::Maximum ? -std::numeric_limits<double>::infinity() : 0.0;
    for (Index i : owned_entities) {
        const double weight = weighting == Weight::Unit                 ? 1.0
                              : field.location() == FieldLocation::Cell ? field.mesh().cellVolume(i)
                                                                        : field.mesh().faceArea(i);
        const double value = weight * valueOf(detail::fieldData(field)[i]);
        if (reduction == Reduction::Maximum)
            local = std::max(local, value);
        else
            local += value;
    }
    double global;
    auto& backend = detail::execution().backend();
    if (reduction == Reduction::Maximum)
        backend.maximum(&local, &global, 1);
    else
        backend.sum(&local, &global, 1);
    return global;
}

template <class T>
std::vector<double>
reduceNorms(std::initializer_list<std::reference_wrapper<const Field<T>>> fields) {
    std::vector<double> local(fields.size(), 0.0), global(fields.size(), 0.0);
    std::size_t index = 0;
    for (const auto& entry : fields) {
        const auto& field = entry.get();
        field.validateStorage();
        if (&field.mesh() != &detail::execution().mesh() ||
            field.location() == FieldLocation::Vertex)
            throw std::invalid_argument("math norm batch requires run mesh cell/face fields");
        const auto& owned_entities = field.location() == FieldLocation::Cell
                                         ? detail::meshData(field.mesh()).owned_cells
                                         : detail::meshData(field.mesh()).owned_faces;
        const auto* values = detail::fieldData(field);
        for (Index i : owned_entities) {
            if constexpr (std::is_same_v<T, double>)
                local[index] += values[i] * values[i];
            else
                local[index] += squaredNorm(values[i]);
        }
        ++index;
    }
    if (!fields.size())
        return global;
    detail::execution().backend().sum(local.data(), global.data(), static_cast<int>(fields.size()));
    for (double& value : global)
        value = std::sqrt(value);
    return global;
}
}
double sum(const ScalarField& field) {
    return reduceField(
        field, [](double v) { return v; }, Reduction::Sum, Weight::Unit);
}
double integral(const ScalarField& field) {
    return reduceField(
        field, [](double v) { return v; }, Reduction::Sum, Weight::Measure);
}
double max(const ScalarField& field) {
    return reduceField(
        field, [](double v) { return v; }, Reduction::Maximum, Weight::Unit);
}
double normL2(const ScalarField& field) {
    return std::sqrt(reduceField(
        field, [](double v) { return v * v; }, Reduction::Sum, Weight::Unit));
}
double normL2(const VectorField& field) {
    return std::sqrt(reduceField(
        field, [](Vec3 v) { return squaredNorm(v); }, Reduction::Sum, Weight::Unit));
}
std::vector<double>
normsL2(std::initializer_list<std::reference_wrapper<const ScalarField>> fields) {
    return reduceNorms(fields);
}
std::vector<double>
normsL2(std::initializer_list<std::reference_wrapper<const VectorField>> fields) {
    return reduceNorms(fields);
}

void add(const ScalarField& increment, ScalarField& target, FaceRegion region) {
    if (&increment.mesh() != &target.mesh() || target.location() != FieldLocation::Face ||
        increment.location() != FieldLocation::Face ||
        &target.mesh() != &detail::execution().mesh())
        throw std::invalid_argument("face increment requires matching run mesh face fields");
    if (region != FaceRegion::All && region != FaceRegion::Interior)
        throw std::invalid_argument("invalid face region");
    auto& backend = detail::execution().backend();
    backend.synchronize(const_cast<ScalarField&>(increment));
    for (Index face : detail::meshData(target.mesh()).owned_faces) {
        if (region == FaceRegion::Interior && target.mesh().boundaryFace(face))
            continue;
        detail::fieldData(target)[face] += detail::fieldData(increment)[face];
    }
    backend.synchronize(target);
}
void subtract(const ScalarField& increment, ScalarField& target, FaceRegion region) {
    const auto negative = -increment;
    add(negative, target, region);
}
void subtract(const ScalarField& coefficient, const VectorField& gradient, VectorField& target) {
    target.addProduct(-1.0, coefficient, gradient);
    detail::execution().backend().synchronize(target);
}

} // namespace babelsim::math
