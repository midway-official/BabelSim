#pragma once

#include "babelsim/operators.h"
#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "internal/boundary_evaluation.h"

#include <limits>
#include <stdexcept>

namespace babelsim::detail::operators {

// The caller exchanges the reconstructed field before using remote cells.
void gradientForExchange(const ScalarField&, VectorField&, GradientMethod);
void gradientForExchange(const VectorField&, TensorField&, GradientMethod);

// Private face sampling rules shared by explicit operators and matrix assembly.
// These functions only read synchronized values; communication stays in the backend.
template <typename T>
void requireField(const Field<T>& field, const Mesh& mesh, FieldLocation location,
                  const char* label) {
    if (&field.mesh() != &mesh || field.location() != location) {
        throw std::invalid_argument(std::string(label) + " field location or mesh is invalid");
    }
    field.validateStorage();
}

inline double ownerWeight(const Mesh& mesh, Index face) {
    const auto f = static_cast<std::size_t>(face);
    return detail::meshData(mesh).face_owner_weights[f];
}

template <typename T>
T internalFaceValue(const Field<T>& field, Index face) {
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    const double weight = ownerWeight(mesh, face);
    return weight * detail::fieldData(field)[owner] +
           (1.0 - weight) * detail::fieldData(field)[neighbour];
}

inline double correctedInternalFaceValue(const ScalarField& field, const Vec3* field_gradient,
                                         Index face) {
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    const double weight = ownerWeight(mesh, face);
    const Vec3 face_gradient =
        weight * field_gradient[owner] + (1.0 - weight) * field_gradient[neighbour];
    return internalFaceValue(field, face) +
           dot(detail::meshData(mesh).face_skewness[f], face_gradient);
}

inline Vec3 correctedInternalFaceValue(const VectorField& field, const Tensor3* field_gradient,
                                       Index face) {
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    const double weight = ownerWeight(mesh, face);
    Vec3 value = internalFaceValue(field, face);
    for (std::size_t component = 0; component < 3; ++component) {
        const Vec3 face_gradient = weight * field_gradient[owner][component] +
                                   (1.0 - weight) * field_gradient[neighbour][component];
        value[component] += dot(detail::meshData(mesh).face_skewness[f], face_gradient);
    }
    return value;
}

inline double linearUpwindCorrection(const ScalarField& field, const VectorField& field_gradient,
                                     Index cell, Index face) {
    return dot(detail::fieldData(field_gradient)[cell],
               detail::meshData(field.mesh()).face_centres[static_cast<std::size_t>(face)] -
                   detail::meshData(field.mesh()).cell_centres[static_cast<std::size_t>(cell)]);
}

inline Vec3 linearUpwindCorrection(const VectorField& field, const TensorField& field_gradient,
                                   Index cell, Index face) {
    const Vec3 offset =
        detail::meshData(field.mesh()).face_centres[static_cast<std::size_t>(face)] -
        detail::meshData(field.mesh()).cell_centres[static_cast<std::size_t>(cell)];
    Vec3 correction{};
    for (std::size_t component = 0; component < 3; ++component) {
        correction[component] = dot(detail::fieldData(field_gradient)[cell][component], offset);
    }
    return correction;
}

inline double
reconstructedBoundaryValue(const ScalarField& field, const Vec3* field_gradient, Index face,
                           double outward_flux = std::numeric_limits<double>::quiet_NaN()) {
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const auto& condition = detail::FieldAccess::condition(field, face, outward_flux);
    if (condition.type == BoundaryType::FixedValue ||
        (condition.type == BoundaryType::InletOutlet && outward_flux < 0.0)) {
        return condition.value;
    }
    const Vec3 normal = mesh.faceNormal(face);
    const double prescribed_normal_gradient =
        condition.type == BoundaryType::FixedGradient ? condition.value : 0.0;
    const Vec3 constrained_gradient =
        field_gradient[owner] +
        (prescribed_normal_gradient - dot(field_gradient[owner], normal)) * normal;
    return detail::fieldData(field)[owner] +
           dot(constrained_gradient,
               detail::meshData(mesh).face_centres[f] -
                   detail::meshData(mesh).cell_centres[static_cast<std::size_t>(owner)]);
}

inline Vec3
reconstructedBoundaryValue(const VectorField& field, const Tensor3* field_gradient, Index face,
                           double outward_flux = std::numeric_limits<double>::quiet_NaN()) {
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const auto& condition = detail::FieldAccess::condition(field, face, outward_flux);
    if (condition.type == BoundaryType::FixedValue ||
        (condition.type == BoundaryType::InletOutlet && outward_flux < 0.0)) {
        return condition.value;
    }
    const Vec3 normal = mesh.faceNormal(face);
    const Vec3 offset = detail::meshData(mesh).face_centres[f] -
                        detail::meshData(mesh).cell_centres[static_cast<std::size_t>(owner)];
    Vec3 value = detail::fieldData(field)[owner];
    for (std::size_t component = 0; component < 3; ++component) {
        const double prescribed_normal_gradient =
            condition.type == BoundaryType::FixedGradient ? condition.value[component] : 0.0;
        const Vec3 constrained_gradient =
            field_gradient[owner][component] +
            (prescribed_normal_gradient - dot(field_gradient[owner][component], normal)) * normal;
        value[component] += dot(constrained_gradient, offset);
    }
    if (condition.type == BoundaryType::Symmetry) {
        value -= dot(value, normal) * normal;
    }
    return value;
}

inline double correctedFaceValue(const ScalarField& field, const Vec3* field_gradient, Index face,
                                 double outward_flux = std::numeric_limits<double>::quiet_NaN()) {
    return field.mesh().boundaryFace(face)
               ? reconstructedBoundaryValue(field, field_gradient, face, outward_flux)
               : correctedInternalFaceValue(field, field_gradient, face);
}

inline Vec3 correctedFaceValue(const VectorField& field, const Tensor3* field_gradient, Index face,
                               double outward_flux = std::numeric_limits<double>::quiet_NaN()) {
    return field.mesh().boundaryFace(face)
               ? reconstructedBoundaryValue(field, field_gradient, face, outward_flux)
               : correctedInternalFaceValue(field, field_gradient, face);
}

inline double correctedInternalFaceValue(const ScalarField& field,
                                         const VectorField& field_gradient, Index face) {
    return correctedInternalFaceValue(field, detail::fieldData(field_gradient), face);
}

inline Vec3 correctedInternalFaceValue(const VectorField& field, const TensorField& field_gradient,
                                       Index face) {
    return correctedInternalFaceValue(field, detail::fieldData(field_gradient), face);
}

inline double
reconstructedBoundaryValue(const ScalarField& field, const VectorField& field_gradient, Index face,
                           double outward_flux = std::numeric_limits<double>::quiet_NaN()) {
    return reconstructedBoundaryValue(field, detail::fieldData(field_gradient), face, outward_flux);
}

inline Vec3
reconstructedBoundaryValue(const VectorField& field, const TensorField& field_gradient, Index face,
                           double outward_flux = std::numeric_limits<double>::quiet_NaN()) {
    return reconstructedBoundaryValue(field, detail::fieldData(field_gradient), face, outward_flux);
}

inline double correctedFaceValue(const ScalarField& field, const VectorField& field_gradient,
                                 Index face,
                                 double outward_flux = std::numeric_limits<double>::quiet_NaN()) {
    return correctedFaceValue(field, detail::fieldData(field_gradient), face, outward_flux);
}

inline Vec3 correctedFaceValue(const VectorField& field, const TensorField& field_gradient,
                               Index face,
                               double outward_flux = std::numeric_limits<double>::quiet_NaN()) {
    return correctedFaceValue(field, detail::fieldData(field_gradient), face, outward_flux);
}

template <typename T>
T interpolatedFaceValue(const Field<T>& field, Index face,
                        double flux_value = std::numeric_limits<double>::quiet_NaN()) {
    const Mesh& mesh = field.mesh();
    return mesh.boundaryFace(face) ? boundaryFaceValue(field, face, flux_value)
                                   : internalFaceValue(field, face);
}

template <typename T>
struct GradientFieldType;

template <>
struct GradientFieldType<double> {
    using Type = VectorField;
};

template <>
struct GradientFieldType<Vec3> {
    using Type = TensorField;
};

inline bool correctedDiffusion(DiffusionMethod method) {
    return method == DiffusionMethod::Corrected || method == DiffusionMethod::LimitedCorrected;
}

} // namespace babelsim::detail::operators
