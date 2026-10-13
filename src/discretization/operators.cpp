#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/operators.h"
#include "internal/boundary_evaluation.h"
#include "operator_support.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace babelsim {
using namespace detail::operators;
namespace {

template <typename T>
void interpolateImpl(const Field<T>& cell, Field<T>& face, InterpolationMethod method) {
    const Mesh& mesh = cell.mesh();
    requireField(cell, mesh, FieldLocation::Cell, "source");
    requireField(face, mesh, FieldLocation::Face, "destination");
    if (method != InterpolationMethod::Linear) {
        throw std::invalid_argument("unsupported interpolation method");
    }
    for (Index f = 0; f < mesh.faceCount(); ++f) {
        detail::fieldData(face)[f] = interpolatedFaceValue(cell, f);
    }
}

template <typename T>
void addFaceDivergence(const Mesh& mesh, Index face, const T& integrated_flux, Field<T>& result) {
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    detail::fieldData(result)[owner] +=
        integrated_flux *
        detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(owner)];
    if (neighbour != invalid_index) {
        detail::fieldData(result)[neighbour] -=
            integrated_flux *
            detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
    }
}

double limitNonOrthogonalCorrection(double orthogonal_flux, double correction,
                                    DiffusionMethod method) {
    if (method == DiffusionMethod::Corrected) {
        return correction;
    }
    if (method != DiffusionMethod::LimitedCorrected) {
        return 0.0;
    }
    // 与 limited corrected 0.5 的稳定性目标一致：显式交叉扩散不超过隐式正交通量。
    const double bound = std::abs(orthogonal_flux);
    return std::clamp(correction, -bound, bound);
}

} // 匿名命名空间

void interpolate(const ScalarField& cell, ScalarField& face, InterpolationMethod method,
                 GradientMethod gradient_method) {
    if (method == InterpolationMethod::Linear) {
        interpolateImpl(cell, face, method);
        return;
    }
    if (method != InterpolationMethod::Corrected) {
        throw std::invalid_argument("unsupported interpolation method");
    }
    const Mesh& mesh = cell.mesh();
    requireField(face, mesh, FieldLocation::Face, "destination");
    VectorField cell_gradient(mesh, FieldLocation::Cell, "grad(" + cell.name() + ')');
    gradient(cell, cell_gradient, gradient_method);
    for (Index index = 0; index < mesh.faceCount(); ++index) {
        detail::fieldData(face)[index] = correctedFaceValue(cell, cell_gradient, index);
    }
}

void interpolate(const VectorField& cell, VectorField& face, InterpolationMethod method,
                 GradientMethod gradient_method) {
    if (method == InterpolationMethod::Linear) {
        interpolateImpl(cell, face, method);
        return;
    }
    if (method != InterpolationMethod::Corrected) {
        throw std::invalid_argument("unsupported interpolation method");
    }
    const Mesh& mesh = cell.mesh();
    requireField(face, mesh, FieldLocation::Face, "destination");
    TensorField cell_gradient(mesh, FieldLocation::Cell, "grad(" + cell.name() + ')');
    gradient(cell, cell_gradient, gradient_method);
    for (Index index = 0; index < mesh.faceCount(); ++index) {
        detail::fieldData(face)[index] = correctedFaceValue(cell, cell_gradient, index);
    }
}

void reconstruct(const ScalarField& cell, const VectorField& cell_gradient, ScalarField& face) {
    const Mesh& mesh = cell.mesh();
    requireField(cell, mesh, FieldLocation::Cell, "source");
    requireField(cell_gradient, mesh, FieldLocation::Cell, "gradient");
    requireField(face, mesh, FieldLocation::Face, "reconstruction");
    for (Index index = 0; index < mesh.faceCount(); ++index) {
        detail::fieldData(face)[index] = correctedFaceValue(cell, cell_gradient, index);
    }
}

void reconstruct(const VectorField& cell, const TensorField& cell_gradient, VectorField& face) {
    const Mesh& mesh = cell.mesh();
    requireField(cell, mesh, FieldLocation::Cell, "source");
    requireField(cell_gradient, mesh, FieldLocation::Cell, "gradient");
    requireField(face, mesh, FieldLocation::Face, "reconstruction");
    for (Index index = 0; index < mesh.faceCount(); ++index) {
        detail::fieldData(face)[index] = correctedFaceValue(cell, cell_gradient, index);
    }
}

double integratedNormalGradient(const ScalarField& scalar, const VectorField& scalar_gradient,
                                Index face, DiffusionMethod method) {
    const Mesh& mesh = scalar.mesh();
    requireField(scalar, mesh, FieldLocation::Cell, "scalar");
    requireField(scalar_gradient, mesh, FieldLocation::Cell, "scalar gradient");
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner.at(f);
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    if (neighbour != invalid_index) {
        const double orthogonal =
            detail::meshData(mesh).face_orthogonal_coefficients[f] *
            (detail::fieldData(scalar)[neighbour] - detail::fieldData(scalar)[owner]);
        if (!correctedDiffusion(method))
            return orthogonal;
        const double weight = ownerWeight(mesh, face);
        const Vec3 face_gradient = weight * detail::fieldData(scalar_gradient)[owner] +
                                   (1.0 - weight) * detail::fieldData(scalar_gradient)[neighbour];
        return orthogonal + limitNonOrthogonalCorrection(
                                orthogonal,
                                dot(detail::meshData(mesh).face_non_orthogonal[f], face_gradient),
                                method);
    }

    const auto& condition = detail::FieldAccess::condition(scalar, face);
    if (condition.type == BoundaryType::FixedGradient) {
        return condition.value * detail::meshData(mesh).face_areas[f];
    }
    if (condition.type != BoundaryType::FixedValue) {
        return 0.0;
    }
    const double orthogonal = detail::meshData(mesh).face_orthogonal_coefficients[f] *
                              (condition.value - detail::fieldData(scalar)[owner]);
    return correctedDiffusion(method)
               ? orthogonal +
                     limitNonOrthogonalCorrection(orthogonal,
                                                  dot(detail::meshData(mesh).face_non_orthogonal[f],
                                                      detail::fieldData(scalar_gradient)[owner]),
                                                  method)
               : orthogonal;
}

Vec3 integratedNormalGradient(const VectorField& vector, const TensorField& vector_gradient,
                              Index face, DiffusionMethod method) {
    const Mesh& mesh = vector.mesh();
    requireField(vector, mesh, FieldLocation::Cell, "vector");
    requireField(vector_gradient, mesh, FieldLocation::Cell, "vector gradient");
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner.at(f);
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    Vec3 result{};
    if (neighbour != invalid_index) {
        const double weight = ownerWeight(mesh, face);
        for (std::size_t component = 0; component < 3; ++component) {
            const double orthogonal = detail::meshData(mesh).face_orthogonal_coefficients[f] *
                                      (detail::fieldData(vector)[neighbour][component] -
                                       detail::fieldData(vector)[owner][component]);
            result[component] = orthogonal;
            if (correctedDiffusion(method)) {
                const Vec3 face_gradient =
                    weight * detail::fieldData(vector_gradient)[owner][component] +
                    (1.0 - weight) * detail::fieldData(vector_gradient)[neighbour][component];
                result[component] += limitNonOrthogonalCorrection(
                    orthogonal, dot(detail::meshData(mesh).face_non_orthogonal[f], face_gradient),
                    method);
            }
        }
        return result;
    }

    const auto& condition = detail::FieldAccess::condition(vector, face);
    if (condition.type == BoundaryType::FixedGradient) {
        return detail::meshData(mesh).face_areas[f] * condition.value;
    }
    if (condition.type != BoundaryType::FixedValue && condition.type != BoundaryType::Symmetry) {
        return {};
    }
    const Vec3 boundary_value = boundaryFaceValue(vector, face);
    const Vec3 normal = mesh.faceNormal(face);
    for (std::size_t component = 0; component < 3; ++component) {
        const double orthogonal =
            detail::meshData(mesh).face_orthogonal_coefficients[f] *
            (boundary_value[component] - detail::fieldData(vector)[owner][component]);
        result[component] = orthogonal;
        if (correctedDiffusion(method)) {
            Vec3 boundary_gradient = detail::fieldData(vector_gradient)[owner][component];
            if (condition.type == BoundaryType::Symmetry) {
                boundary_gradient -= dot(boundary_gradient, normal) * normal;
            }
            result[component] += limitNonOrthogonalCorrection(
                orthogonal, dot(detail::meshData(mesh).face_non_orthogonal[f], boundary_gradient),
                method);
        }
    }
    return result;
}

void flux(const VectorField& velocity, ScalarField& face_flux, InterpolationMethod method,
          GradientMethod gradient_method) {
    const Mesh& mesh = velocity.mesh();
    requireField(face_flux, mesh, FieldLocation::Face, "flux");
    if (velocity.location() == FieldLocation::Face) {
        requireField(velocity, mesh, FieldLocation::Face, "face vector");
        for (Index face : detail::meshData(mesh).owned_faces)
            detail::fieldData(face_flux)[face] =
                dot(detail::fieldData(velocity)[face], mesh.faceAreaVector(face));
        return;
    }
    requireField(velocity, mesh, FieldLocation::Cell, "velocity");
    if (method != InterpolationMethod::Linear && method != InterpolationMethod::Corrected) {
        throw std::invalid_argument("unsupported flux interpolation method");
    }
    // Read the input's existing direction context, never the output buffer.
    // On first use, seed directions from owner velocities before reconstructing
    // gradients, which also need a defined inletOutlet branch.
    if (!detail::FieldAccess::hasBoundaryFlux(velocity)) {
        ScalarField initialFlux(mesh, FieldLocation::Face, "flux.initialDirection");
        for (Index face = 0; face < mesh.faceCount(); ++face)
            if (mesh.boundaryFace(face))
                detail::fieldData(initialFlux)[face] =
                    dot(detail::fieldData(velocity)[mesh.owner(face)], mesh.faceAreaVector(face));
        const_cast<VectorField&>(velocity).setBoundaryFlux(initialFlux);
    }
    std::optional<TensorField> velocity_gradient;
    if (method == InterpolationMethod::Corrected) {
        velocity_gradient.emplace(mesh, FieldLocation::Cell, "grad(" + velocity.name() + ')');
        gradient(velocity, *velocity_gradient, gradient_method);
    }
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Vec3 face_velocity = velocity_gradient
                                       ? correctedFaceValue(velocity, *velocity_gradient, face)
                                       : interpolatedFaceValue(velocity, face);
        detail::fieldData(face_flux)[face] =
            dot(face_velocity, detail::meshData(mesh).face_area_vectors[f]);
    }
}

void diffusionFlux(const ScalarField& face_diffusivity, const ScalarField& scalar,
                   const VectorField& scalar_gradient, ScalarField& face_flux,
                   DiffusionMethod diffusion_method) {
    const Mesh& mesh = scalar.mesh();
    requireField(face_diffusivity, mesh, FieldLocation::Face, "face diffusivity");
    requireField(scalar, mesh, FieldLocation::Cell, "diffusion field");
    requireField(scalar_gradient, mesh, FieldLocation::Cell, "diffusion gradient");
    requireField(face_flux, mesh, FieldLocation::Face, "diffusion flux");
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const double diffusivity = detail::fieldData(face_diffusivity)[face];
        if (!(diffusivity >= 0.0) || !std::isfinite(diffusivity)) {
            throw std::invalid_argument("face diffusivity must be non-negative and finite");
        }
        detail::fieldData(face_flux)[face] =
            diffusivity * integratedNormalGradient(scalar, scalar_gradient, face, diffusion_method);
    }
}

void divergence(const VectorField& vector, ScalarField& result, InterpolationMethod method,
                GradientMethod gradient_method) {
    const Mesh& mesh = vector.mesh();
    requireField(vector, mesh, FieldLocation::Cell, "vector");
    requireField(result, mesh, FieldLocation::Cell, "divergence");
    if (method != InterpolationMethod::Linear && method != InterpolationMethod::Corrected) {
        throw std::invalid_argument("unsupported divergence interpolation method");
    }
    std::optional<TensorField> vector_gradient;
    if (method == InterpolationMethod::Corrected) {
        vector_gradient.emplace(mesh, FieldLocation::Cell, "grad(" + vector.name() + ')');
        gradient(vector, *vector_gradient, gradient_method);
    }
    result.fill(0.0);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Vec3 face_value = vector_gradient ? correctedFaceValue(vector, *vector_gradient, face)
                                                : interpolatedFaceValue(vector, face);
        const double integrated_flux = dot(face_value, detail::meshData(mesh).face_area_vectors[f]);
        addFaceDivergence(mesh, face, integrated_flux, result);
    }
}

void divergence(const ScalarField& face_flux, ScalarField& result) {
    const Mesh& mesh = face_flux.mesh();
    requireField(face_flux, mesh, FieldLocation::Face, "face flux");
    requireField(result, mesh, FieldLocation::Cell, "divergence");
    result.fill(0.0);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        addFaceDivergence(mesh, face, detail::fieldData(face_flux)[face], result);
    }
}

template <typename T>
void convectionImpl(const ScalarField& face_flux, const Field<T>& transported, Field<T>& result,
                    ConvectionMethod method, InterpolationMethod interpolation_method,
                    GradientMethod gradient_method) {
    const Mesh& mesh = transported.mesh();
    requireField(face_flux, mesh, FieldLocation::Face, "face flux");
    requireField(transported, mesh, FieldLocation::Cell, "transported");
    requireField(result, mesh, FieldLocation::Cell, "convection");
    if (method != ConvectionMethod::Upwind && method != ConvectionMethod::LinearUpwind &&
        method != ConvectionMethod::Central) {
        throw std::invalid_argument("unsupported convection method");
    }
    if (interpolation_method != InterpolationMethod::Linear &&
        interpolation_method != InterpolationMethod::Corrected) {
        throw std::invalid_argument("unsupported convection interpolation method");
    }
    using GradientField = typename GradientFieldType<T>::Type;
    std::optional<GradientField> transported_gradient;
    if (method == ConvectionMethod::LinearUpwind ||
        (method == ConvectionMethod::Central &&
         interpolation_method == InterpolationMethod::Corrected)) {
        transported_gradient.emplace(mesh, FieldLocation::Cell, "grad(" + transported.name() + ')');
        gradient(transported, *transported_gradient, gradient_method);
    }
    result.fill(T{});
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const std::size_t index = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[index];
        const Index neighbour = detail::meshData(mesh).face_neighbour[index];
        const double flux_value = detail::fieldData(face_flux)[face];
        T face_value{};
        if (method == ConvectionMethod::Upwind || method == ConvectionMethod::LinearUpwind) {
            if (neighbour != invalid_index) {
                const Index upwind = flux_value >= 0.0 ? owner : neighbour;
                face_value = detail::fieldData(transported)[upwind];
                if (method == ConvectionMethod::LinearUpwind) {
                    face_value +=
                        linearUpwindCorrection(transported, *transported_gradient, upwind, face);
                }
            } else {
                if (flux_value >= 0.0) {
                    face_value = detail::fieldData(transported)[owner];
                    if (method == ConvectionMethod::LinearUpwind) {
                        face_value +=
                            linearUpwindCorrection(transported, *transported_gradient, owner, face);
                    }
                } else {
                    face_value = boundaryFaceValue(transported, face, flux_value);
                }
            }
        } else if (transported_gradient) {
            face_value = correctedFaceValue(transported, *transported_gradient, face, flux_value);
        } else {
            face_value = interpolatedFaceValue(transported, face, flux_value);
        }
        addFaceDivergence(mesh, face, flux_value * face_value, result);
    }
}

void convection(const ScalarField& face_flux, const ScalarField& transported, ScalarField& result,
                ConvectionMethod method, InterpolationMethod interpolation_method,
                GradientMethod gradient_method) {
    convectionImpl(face_flux, transported, result, method, interpolation_method, gradient_method);
}

void convection(const ScalarField& face_flux, const VectorField& transported, VectorField& result,
                ConvectionMethod method, InterpolationMethod interpolation_method,
                GradientMethod gradient_method) {
    convectionImpl(face_flux, transported, result, method, interpolation_method, gradient_method);
}

void laplacianImpl(const ScalarField& scalar, ScalarField& result, double constant_diffusivity,
                   const ScalarField* face_diffusivity, GradientMethod gradient_method,
                   DiffusionMethod diffusion_method) {
    const Mesh& mesh = scalar.mesh();
    requireField(scalar, mesh, FieldLocation::Cell, "scalar");
    requireField(result, mesh, FieldLocation::Cell, "laplacian");
    if (!(constant_diffusivity >= 0.0) || !std::isfinite(constant_diffusivity)) {
        throw std::invalid_argument("diffusivity must be non-negative and finite");
    }
    if (face_diffusivity != nullptr) {
        requireField(*face_diffusivity, mesh, FieldLocation::Face, "face diffusivity");
    }
    std::optional<VectorField> cell_gradient;
    if (correctedDiffusion(diffusion_method)) {
        cell_gradient.emplace(mesh, FieldLocation::Cell, "grad(" + scalar.name() + ')');
        gradient(scalar, *cell_gradient, gradient_method);
    } else if (diffusion_method != DiffusionMethod::Orthogonal) {
        throw std::invalid_argument("unsupported diffusion method");
    }

    result.fill(0.0);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        double integrated_flux = 0.0;
        if (neighbour != invalid_index) {
            integrated_flux =
                detail::meshData(mesh).face_orthogonal_coefficients[f] *
                (detail::fieldData(scalar)[neighbour] - detail::fieldData(scalar)[owner]);
            if (correctedDiffusion(diffusion_method)) {
                integrated_flux =
                    integratedNormalGradient(scalar, *cell_gradient, face, diffusion_method);
            }
        } else {
            const auto& condition = detail::FieldAccess::condition(scalar, face);
            if (condition.type == BoundaryType::FixedValue) {
                integrated_flux =
                    detail::meshData(mesh).face_orthogonal_coefficients[f] *
                    (boundaryFaceValue(scalar, face, -1.0) - detail::fieldData(scalar)[owner]);
                if (correctedDiffusion(diffusion_method)) {
                    integrated_flux =
                        integratedNormalGradient(scalar, *cell_gradient, face, diffusion_method);
                }
            } else if (condition.type == BoundaryType::FixedGradient) {
                integrated_flux = condition.value * detail::meshData(mesh).face_areas[f];
            }
        }
        const double diffusivity = face_diffusivity == nullptr
                                       ? constant_diffusivity
                                       : detail::fieldData((*face_diffusivity))[face];
        if (!(diffusivity >= 0.0) || !std::isfinite(diffusivity)) {
            throw std::invalid_argument("face diffusivity must be non-negative and finite");
        }
        addFaceDivergence(mesh, face, diffusivity * integrated_flux, result);
    }
}

void laplacian(const ScalarField& scalar, ScalarField& result, GradientMethod gradient_method,
               DiffusionMethod diffusion_method) {
    laplacianImpl(scalar, result, 1.0, nullptr, gradient_method, diffusion_method);
}

void laplacian(double diffusivity, const ScalarField& scalar, ScalarField& result,
               GradientMethod gradient_method, DiffusionMethod diffusion_method) {
    laplacianImpl(scalar, result, diffusivity, nullptr, gradient_method, diffusion_method);
}

void laplacian(const ScalarField& face_diffusivity, const ScalarField& scalar, ScalarField& result,
               GradientMethod gradient_method, DiffusionMethod diffusion_method) {
    laplacianImpl(scalar, result, 1.0, &face_diffusivity, gradient_method, diffusion_method);
}

void divergence(const TensorField& tensor, VectorField& result,
                InterpolationMethod interpolation_method, GradientMethod gradient_method) {
    const Mesh& mesh = tensor.mesh();
    requireField(tensor, mesh, FieldLocation::Cell, "tensor");
    requireField(result, mesh, FieldLocation::Cell, "tensor divergence");
    VectorField row(mesh, FieldLocation::Cell, tensor.name());
    VectorField face_values(mesh, FieldLocation::Face);
    row.useCalculatedBoundary();
    result.fill({});
    // Reuse the vector reconstruction for each tensor row, including all boundary traces.
    for (int component = 0; component < 3; ++component) {
        row.evaluate(tensor, [component](const Tensor3& value) { return value[component]; });
        interpolate(row, face_values, interpolation_method, gradient_method);
        for (Index face = 0; face < mesh.faceCount(); ++face) {
            const Index owner = mesh.owner(face), neighbour = mesh.neighbour(face);
            const double contribution =
                dot(detail::fieldData(static_cast<const VectorField&>(face_values))[face],
                    mesh.faceAreaVector(face));
            detail::fieldData(result)[owner][component] += contribution / mesh.cellVolume(owner);
            if (neighbour != invalid_index)
                detail::fieldData(result)[neighbour][component] -=
                    contribution / mesh.cellVolume(neighbour);
        }
    }
}

} // babelsim 命名空间
