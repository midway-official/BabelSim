#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/operators.h"
#include "internal/boundary_evaluation.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace babelsim {
namespace {

template <typename T>
void requireField(
    const Field<T>& field,
    const Mesh& mesh,
    FieldLocation location,
    const char* label)
{
    if (&field.mesh() != &mesh || field.location() != location) {
        throw std::invalid_argument(std::string(label) + " field location or mesh is invalid");
    }
    field.validateStorage();
}

double ownerWeight(const Mesh& mesh, Index face) {
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
    return weight * detail::fieldData(field)[owner] + (1.0 - weight) * detail::fieldData(field)[neighbour];
}

double correctedInternalFaceValue(
    const ScalarField& field,
    const VectorField& field_gradient,
    Index face)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    const double weight = ownerWeight(mesh, face);
    const Vec3 face_gradient =
        weight * detail::fieldData(field_gradient)[owner] +
        (1.0 - weight) * detail::fieldData(field_gradient)[neighbour];
    return internalFaceValue(field, face) +
        dot(detail::meshData(mesh).face_skewness[f], face_gradient);
}

Vec3 correctedInternalFaceValue(
    const VectorField& field,
    const TensorField& field_gradient,
    Index face)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    const double weight = ownerWeight(mesh, face);
    Vec3 value = internalFaceValue(field, face);
    for (std::size_t component = 0; component < 3; ++component) {
        const Vec3 face_gradient =
            weight * detail::fieldData(field_gradient)[owner][component] +
            (1.0 - weight) * detail::fieldData(field_gradient)[neighbour][component];
        value[component] += dot(detail::meshData(mesh).face_skewness[f], face_gradient);
    }
    return value;
}

double linearUpwindCorrection(
    const ScalarField& field,
    const VectorField& field_gradient,
    Index cell,
    Index face)
{
    return dot(
        detail::fieldData(field_gradient)[cell],
        detail::meshData(field.mesh()).face_centres[static_cast<std::size_t>(face)] -
            detail::meshData(field.mesh()).cell_centres[static_cast<std::size_t>(cell)]);
}

Vec3 linearUpwindCorrection(
    const VectorField& field,
    const TensorField& field_gradient,
    Index cell,
    Index face)
{
    const Vec3 offset =
        detail::meshData(field.mesh()).face_centres[static_cast<std::size_t>(face)] -
        detail::meshData(field.mesh()).cell_centres[static_cast<std::size_t>(cell)];
    Vec3 correction{};
    for (std::size_t component = 0; component < 3; ++component) {
        correction[component] = dot(
            detail::fieldData(field_gradient)[cell][component], offset);
    }
    return correction;
}

double reconstructedBoundaryValue(
    const ScalarField& field,
    const VectorField& field_gradient,
    Index face,
    double outward_flux = 0.0)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const auto& condition = field.boundary(detail::meshData(mesh).face_patch[f]);
    if (condition.type == BoundaryType::FixedValue ||
        (condition.type == BoundaryType::InletOutlet && outward_flux < 0.0)) {
        return condition.value;
    }
    const Vec3 normal = mesh.faceNormal(face);
    const double prescribed_normal_gradient =
        condition.type == BoundaryType::FixedGradient ? condition.value : 0.0;
    const Vec3 constrained_gradient = detail::fieldData(field_gradient)[owner] +
        (prescribed_normal_gradient - dot(detail::fieldData(field_gradient)[owner], normal)) * normal;
    return detail::fieldData(field)[owner] + dot(
        constrained_gradient,
        detail::meshData(mesh).face_centres[f] - detail::meshData(mesh).cell_centres[static_cast<std::size_t>(owner)]);
}

Vec3 reconstructedBoundaryValue(
    const VectorField& field,
    const TensorField& field_gradient,
    Index face,
    double outward_flux = 0.0)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const auto& condition = field.boundary(detail::meshData(mesh).face_patch[f]);
    if (condition.type == BoundaryType::FixedValue ||
        (condition.type == BoundaryType::InletOutlet && outward_flux < 0.0)) {
        return condition.value;
    }
    const Vec3 normal = mesh.faceNormal(face);
    const Vec3 offset =
        detail::meshData(mesh).face_centres[f] - detail::meshData(mesh).cell_centres[static_cast<std::size_t>(owner)];
    Vec3 value = detail::fieldData(field)[owner];
    for (std::size_t component = 0; component < 3; ++component) {
        const double prescribed_normal_gradient =
            condition.type == BoundaryType::FixedGradient
            ? condition.value[component] : 0.0;
        const Vec3 constrained_gradient = detail::fieldData(field_gradient)[owner][component] +
            (prescribed_normal_gradient -
             dot(detail::fieldData(field_gradient)[owner][component], normal)) * normal;
        value[component] += dot(constrained_gradient, offset);
    }
    if (condition.type == BoundaryType::Symmetry) {
        value -= dot(value, normal) * normal;
    }
    return value;
}

double correctedFaceValue(
    const ScalarField& field,
    const VectorField& field_gradient,
    Index face,
    double outward_flux = 0.0)
{
    return field.mesh().boundaryFace(face)
        ? reconstructedBoundaryValue(field, field_gradient, face, outward_flux)
        : correctedInternalFaceValue(field, field_gradient, face);
}

Vec3 correctedFaceValue(
    const VectorField& field,
    const TensorField& field_gradient,
    Index face,
    double outward_flux = 0.0)
{
    return field.mesh().boundaryFace(face)
        ? reconstructedBoundaryValue(field, field_gradient, face, outward_flux)
        : correctedInternalFaceValue(field, field_gradient, face);
}

template <typename T>
T interpolatedFaceValue(const Field<T>& field, Index face, double flux_value = 0.0) {
    const Mesh& mesh = field.mesh();
    return mesh.boundaryFace(face)
        ? boundaryFaceValue(field, face, flux_value)
        : internalFaceValue(field, face);
}

template <typename T>
void interpolateImpl(
    const Field<T>& cell,
    Field<T>& face,
    InterpolationMethod method)
{
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

Vec3 toVec3(const Eigen::Vector3d& value) {
    return {value.x(), value.y(), value.z()};
}

Eigen::Vector3d toEigen(const Vec3& value) {
    return {value.x, value.y, value.z};
}

void greenGaussGradient(const ScalarField& scalar, VectorField& result) {
    const Mesh& mesh = scalar.mesh();
    result.fill({});
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const double value = interpolatedFaceValue(scalar, face);
        const Vec3 contribution = value * detail::meshData(mesh).face_area_vectors[f];
        detail::fieldData(result)[owner] += contribution *
            detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(owner)];
        if (neighbour != invalid_index) {
            detail::fieldData(result)[neighbour] -= contribution *
                detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
        }
    }

    // 用初始高斯梯度把交点面值重构到真实面中心。一次显式修正已能恢复光滑网格上的
    // 二阶面值，同时只把计算模板扩展到两层邻居，符合当前 MPI halo 宽度。
    if (mesh.orthogonalGeometry()) {
        return;
    }
    const VectorField initial_gradient = result;
    result.fill({});
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const double value = correctedFaceValue(scalar, initial_gradient, face);
        const Vec3 contribution = value * detail::meshData(mesh).face_area_vectors[f];
        detail::fieldData(result)[owner] += contribution *
            detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(owner)];
        if (neighbour != invalid_index) {
            detail::fieldData(result)[neighbour] -= contribution *
                detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
        }
    }
}

void leastSquaresGradient(const ScalarField& scalar, VectorField& result) {
    const Mesh& mesh = scalar.mesh();
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
        Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
        const auto c = static_cast<std::size_t>(cell);
        for (Index face : detail::meshData(mesh).cell_faces[c]) {
            const auto f = static_cast<std::size_t>(face);
            const Index owner = detail::meshData(mesh).face_owner[f];
            const Index neighbour = detail::meshData(mesh).face_neighbour[f];
            Vec3 delta{};
            double difference = 0.0;
            if (neighbour != invalid_index) {
                const Index other = owner == cell ? neighbour : owner;
                delta = detail::meshData(mesh).cell_centres[static_cast<std::size_t>(other)] -
                    detail::meshData(mesh).cell_centres[c];
                difference = detail::fieldData(scalar)[other] - detail::fieldData(scalar)[cell];
            } else {
                const auto type = scalar.boundary(detail::meshData(mesh).face_patch[f]).type;
                const Vec3 offset = detail::meshData(mesh).face_centres[f] - detail::meshData(mesh).cell_centres[c];
                delta = type == BoundaryType::FixedValue
                    ? 2.0 * offset
                    : 2.0 * boundaryNormalDistance(mesh, face) *
                        mesh.faceNormal(face);
                difference = 2.0 * (boundaryFaceValue(scalar, face) - detail::fieldData(scalar)[cell]);
            }
            const double distance_squared = squaredNorm(delta);
            if (!(distance_squared > 0.0)) {
                throw std::runtime_error("least-squares neighbour distance is invalid");
            }
            const Eigen::Vector3d d = toEigen(delta);
            const double weight = 1.0 / distance_squared;
            normal.noalias() += weight * d * d.transpose();
            rhs.noalias() += weight * d * difference;
        }
        Eigen::LDLT<Eigen::Matrix3d> decomposition(normal);
        if (decomposition.info() != Eigen::Success) {
            throw std::runtime_error("least-squares gradient matrix is singular");
        }
        const Eigen::Vector3d value = decomposition.solve(rhs);
        if (decomposition.info() != Eigen::Success || !value.allFinite()) {
            throw std::runtime_error("least-squares gradient solve failed");
        }
        detail::fieldData(result)[cell] = toVec3(value);
    }
}

void greenGaussGradient(const VectorField& vector, TensorField& result) {
    const Mesh& mesh = vector.mesh();
    result.fill({});
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const Vec3 value = interpolatedFaceValue(vector, face);
        for (std::size_t component = 0; component < 3; ++component) {
            const Vec3 contribution = value[component] * detail::meshData(mesh).face_area_vectors[f];
            detail::fieldData(result)[owner][component] += contribution *
                detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(owner)];
            if (neighbour != invalid_index) {
                detail::fieldData(result)[neighbour][component] -= contribution *
                    detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
            }
        }
    }
    if (mesh.orthogonalGeometry()) {
        return;
    }
    const TensorField initial_gradient = result;
    result.fill({});
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const Vec3 value = correctedFaceValue(vector, initial_gradient, face);
        for (std::size_t component = 0; component < 3; ++component) {
            const Vec3 contribution = value[component] * detail::meshData(mesh).face_area_vectors[f];
            detail::fieldData(result)[owner][component] += contribution *
                detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(owner)];
            if (neighbour != invalid_index) {
                detail::fieldData(result)[neighbour][component] -= contribution *
                    detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
            }
        }
    }
}

void leastSquaresGradient(const VectorField& vector, TensorField& result) {
    const Mesh& mesh = vector.mesh();
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d rhs = Eigen::Matrix3d::Zero();
        const auto c = static_cast<std::size_t>(cell);
        for (Index face : detail::meshData(mesh).cell_faces[c]) {
            const auto f = static_cast<std::size_t>(face);
            const Index owner = detail::meshData(mesh).face_owner[f];
            const Index neighbour = detail::meshData(mesh).face_neighbour[f];
            Vec3 delta{};
            Vec3 difference{};
            if (neighbour != invalid_index) {
                const Index other = owner == cell ? neighbour : owner;
                delta = detail::meshData(mesh).cell_centres[static_cast<std::size_t>(other)] -
                    detail::meshData(mesh).cell_centres[c];
                difference = detail::fieldData(vector)[other] - detail::fieldData(vector)[cell];
            } else {
                const auto type = vector.boundary(detail::meshData(mesh).face_patch[f]).type;
                const Vec3 offset = detail::meshData(mesh).face_centres[f] - detail::meshData(mesh).cell_centres[c];
                delta = type == BoundaryType::FixedValue
                    ? 2.0 * offset
                    : 2.0 * boundaryNormalDistance(mesh, face) *
                        mesh.faceNormal(face);
                difference = 2.0 * (boundaryFaceValue(vector, face) - detail::fieldData(vector)[cell]);
            }
            const double distance_squared = squaredNorm(delta);
            if (!(distance_squared > 0.0)) {
                throw std::runtime_error("least-squares neighbour distance is invalid");
            }
            const Eigen::Vector3d d = toEigen(delta);
            const double weight = 1.0 / distance_squared;
            normal.noalias() += weight * d * d.transpose();
            for (std::size_t component = 0; component < 3; ++component) {
                rhs.col(static_cast<Eigen::Index>(component)).noalias() +=
                    weight * d * difference[component];
            }
        }
        Eigen::LDLT<Eigen::Matrix3d> decomposition(normal);
        if (decomposition.info() != Eigen::Success) {
            throw std::runtime_error("least-squares vector-gradient matrix is singular");
        }
        const Eigen::Matrix3d value = decomposition.solve(rhs);
        if (decomposition.info() != Eigen::Success || !value.allFinite()) {
            throw std::runtime_error("least-squares vector-gradient solve failed");
        }
        for (std::size_t component = 0; component < 3; ++component) {
            detail::fieldData(result)[cell][component] = toVec3(
                value.col(static_cast<Eigen::Index>(component)));
        }
    }
}

void addFaceDivergence(
    const Mesh& mesh,
    Index face,
    double integrated_flux,
    ScalarField& result)
{
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    detail::fieldData(result)[owner] += integrated_flux *
        detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(owner)];
    if (neighbour != invalid_index) {
        detail::fieldData(result)[neighbour] -= integrated_flux *
            detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
    }
}

void addFaceDivergence(
    const Mesh& mesh,
    Index face,
    const Vec3& integrated_flux,
    VectorField& result)
{
    const std::size_t index = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[index];
    const Index neighbour = detail::meshData(mesh).face_neighbour[index];
    detail::fieldData(result)[owner] += integrated_flux * detail::meshData(mesh).cell_inverse_volumes[
        static_cast<std::size_t>(owner)];
    if (neighbour != invalid_index) {
        detail::fieldData(result)[neighbour] -= integrated_flux * detail::meshData(mesh).cell_inverse_volumes[
            static_cast<std::size_t>(neighbour)];
    }
}

bool correctedDiffusion(DiffusionMethod method) {
    return method == DiffusionMethod::Corrected ||
        method == DiffusionMethod::LimitedCorrected;
}

double limitNonOrthogonalCorrection(
    double orthogonal_flux,
    double correction,
    DiffusionMethod method)
{
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

template <typename T>
void requireEquation(const DiscreteEquation<T>& equation, const Mesh& mesh) {
    if (equation.mesh != &mesh) {
        throw std::invalid_argument("equation and field use different meshes");
    }
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

template <typename T>
void addConvectionImpl(
    DiscreteEquation<T>& equation,
    const ScalarField& face_flux,
    const Field<T>& transported,
    ConvectionMethod method,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method,
    double flux_scale)
{
    const Mesh& mesh = transported.mesh();
    requireEquation(equation, mesh);
    requireField(face_flux, mesh, FieldLocation::Face, "face flux");
    requireField(transported, mesh, FieldLocation::Cell, "transported");

    if (interpolation_method != InterpolationMethod::Linear &&
        interpolation_method != InterpolationMethod::Corrected) {
        throw std::invalid_argument("unsupported convection interpolation method");
    }
    if (!std::isfinite(flux_scale)) {
        throw std::invalid_argument("convection flux scale must be finite");
    }
    using GradientField = typename GradientFieldType<T>::Type;
    std::optional<GradientField> transported_gradient;
    if (method == ConvectionMethod::LinearUpwind ||
        (method == ConvectionMethod::Central &&
         interpolation_method == InterpolationMethod::Corrected)) {
        transported_gradient.emplace(
            mesh, FieldLocation::Cell, "grad(" + transported.name() + ')');
        gradient(transported, *transported_gradient, gradient_method);
    }

    for (Index face : detail::meshData(mesh).owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const double F = flux_scale * detail::fieldData(face_flux)[face];
        if (neighbour != invalid_index) {
            if (method == ConvectionMethod::Upwind ||
                method == ConvectionMethod::LinearUpwind) {
                equation.diagonal[static_cast<std::size_t>(owner)] += std::max(F, 0.0);
                equation.upper[f] += std::min(F, 0.0);
                equation.diagonal[static_cast<std::size_t>(neighbour)] +=
                    std::max(-F, 0.0);
                equation.lower[f] += std::min(-F, 0.0);
                if (method == ConvectionMethod::LinearUpwind) {
                    const Index upwind = F >= 0.0 ? owner : neighbour;
                    const T correction = linearUpwindCorrection(
                        transported, *transported_gradient, upwind, face);
                    equation.source[static_cast<std::size_t>(owner)] -= F * correction;
                    equation.source[static_cast<std::size_t>(neighbour)] += F * correction;
                }
            } else if (method == ConvectionMethod::Central) {
                const double weight = ownerWeight(mesh, face);
                equation.diagonal[static_cast<std::size_t>(owner)] += F * weight;
                equation.upper[f] += F * (1.0 - weight);
                equation.lower[f] -= F * weight;
                equation.diagonal[static_cast<std::size_t>(neighbour)] -=
                    F * (1.0 - weight);
                if (transported_gradient) {
                    const T correction =
                        correctedInternalFaceValue(
                            transported, *transported_gradient, face) -
                        internalFaceValue(transported, face);
                    equation.source[static_cast<std::size_t>(owner)] -= F * correction;
                    equation.source[static_cast<std::size_t>(neighbour)] += F * correction;
                }
            } else {
                throw std::invalid_argument("unsupported convection method");
            }
            continue;
        }

        const auto& condition = transported.boundary(detail::meshData(mesh).face_patch[f]);
        if (method == ConvectionMethod::Upwind ||
            method == ConvectionMethod::LinearUpwind) {
            if (F >= 0.0) {
                equation.diagonal[static_cast<std::size_t>(owner)] += F;
                if (method == ConvectionMethod::LinearUpwind) {
                    equation.source[static_cast<std::size_t>(owner)] -= F *
                        linearUpwindCorrection(
                            transported, *transported_gradient, owner, face);
                }
            } else {
                equation.source[static_cast<std::size_t>(owner)] -=
                    F * boundaryFaceValue(transported, face, F);
            }
            continue;
        }
        if (method != ConvectionMethod::Central) {
            throw std::invalid_argument("unsupported convection method");
        }
        switch (condition.type) {
            case BoundaryType::FixedValue:
                equation.source[static_cast<std::size_t>(owner)] -=
                    F * condition.value;
                break;
            case BoundaryType::FixedGradient: {
                const double distance = boundaryNormalDistance(mesh, face);
                equation.diagonal[static_cast<std::size_t>(owner)] += F;
                equation.source[static_cast<std::size_t>(owner)] -=
                    F * distance * condition.value;
                break;
            }
            case BoundaryType::ZeroGradient:
            case BoundaryType::Symmetry:
                equation.diagonal[static_cast<std::size_t>(owner)] += F;
                break;
            case BoundaryType::InletOutlet:
                if (F >= 0.0) {
                    equation.diagonal[static_cast<std::size_t>(owner)] += F;
                } else {
                    equation.source[static_cast<std::size_t>(owner)] -=
                        F * condition.value;
                }
                break;
        }
        if (transported_gradient) {
            const T correction =
                reconstructedBoundaryValue(
                    transported, *transported_gradient, face, F) -
                boundaryFaceValue(transported, face, F);
            equation.source[static_cast<std::size_t>(owner)] -= F * correction;
        }
    }
}

template <typename T>
void addTimeDerivativeImpl(
    DiscreteEquation<T>& equation,
    const Field<T>& previous,
    double dt,
    double density,
    TimeMethod method,
    const Field<T>* older)
{
    const Mesh& mesh = previous.mesh();
    requireEquation(equation, mesh);
    requireField(previous, mesh, FieldLocation::Cell, "previous");
    if (method == TimeMethod::Steady) {
        return;
    }
    if (!(dt > 0.0) || !(density > 0.0) || !std::isfinite(dt) ||
        !std::isfinite(density)) {
        throw std::invalid_argument("time step and density must be positive and finite");
    }
    if (method == TimeMethod::BDF2) {
        if (older == nullptr) {
            throw std::invalid_argument("BDF2 requires two previous fields");
        }
        requireField(*older, mesh, FieldLocation::Cell, "older");
    } else if (method != TimeMethod::Euler) {
        throw std::invalid_argument("unsupported time method");
    }

    for (Index cell : detail::meshData(mesh).owned_cells) {
        const auto c = static_cast<std::size_t>(cell);
        const double coefficient = density * detail::meshData(mesh).cell_volumes[c] / dt;
        if (method == TimeMethod::Euler) {
            equation.diagonal[c] += coefficient;
            equation.source[c] += coefficient * detail::fieldData(previous)[cell];
        } else {
            equation.diagonal[c] += 1.5 * coefficient;
            equation.source[c] += coefficient *
                (2.0 * detail::fieldData(previous)[cell] - 0.5 * detail::fieldData((*older))[cell]);
        }
    }
}

template <typename T>
void addFieldTimeDerivativeImpl(
    DiscreteEquation<T>& equation,
    const Field<T>& previous,
    double dt,
    const ScalarField& capacity,
    TimeMethod method,
    const Field<T>* older)
{
    const Mesh& mesh = previous.mesh();
    requireEquation(equation, mesh);
    requireField(previous, mesh, FieldLocation::Cell, "previous");
    requireField(capacity, mesh, FieldLocation::Cell, "time coefficient");
    if (method == TimeMethod::Steady) return;
    if (!(dt > 0.0) || !std::isfinite(dt)) {
        throw std::invalid_argument("time step must be positive and finite");
    }
    if (method == TimeMethod::BDF2) {
        if (older == nullptr) throw std::invalid_argument("BDF2 requires two previous fields");
        requireField(*older, mesh, FieldLocation::Cell, "older");
    } else if (method != TimeMethod::Euler) {
        throw std::invalid_argument("unsupported time method");
    }
    for (Index cell : detail::meshData(mesh).owned_cells) {
        const std::size_t index = static_cast<std::size_t>(cell);
        if (!(detail::fieldData(capacity)[cell] > 0.0) || !std::isfinite(detail::fieldData(capacity)[cell])) {
            throw std::invalid_argument("time coefficient must be positive and finite");
        }
        const double coefficient = detail::fieldData(capacity)[cell] * detail::meshData(mesh).cell_volumes[index] / dt;
        if (method == TimeMethod::Euler) {
            equation.diagonal[index] += coefficient;
            equation.source[index] += coefficient * detail::fieldData(previous)[cell];
        } else {
            equation.diagonal[index] += 1.5 * coefficient;
            equation.source[index] += coefficient *
                (2.0 * detail::fieldData(previous)[cell] - 0.5 * detail::fieldData((*older))[cell]);
        }
    }
}

}  // 匿名命名空间

void interpolate(
    const ScalarField& cell,
    ScalarField& face,
    InterpolationMethod method,
    GradientMethod gradient_method)
{
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

void interpolate(
    const VectorField& cell,
    VectorField& face,
    InterpolationMethod method,
    GradientMethod gradient_method)
{
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

void reconstruct(
    const ScalarField& cell,
    const VectorField& cell_gradient,
    ScalarField& face)
{
    const Mesh& mesh = cell.mesh();
    requireField(cell, mesh, FieldLocation::Cell, "source");
    requireField(cell_gradient, mesh, FieldLocation::Cell, "gradient");
    requireField(face, mesh, FieldLocation::Face, "reconstruction");
    for (Index index = 0; index < mesh.faceCount(); ++index) {
        detail::fieldData(face)[index] = correctedFaceValue(cell, cell_gradient, index);
    }
}

void reconstruct(
    const VectorField& cell,
    const TensorField& cell_gradient,
    VectorField& face)
{
    const Mesh& mesh = cell.mesh();
    requireField(cell, mesh, FieldLocation::Cell, "source");
    requireField(cell_gradient, mesh, FieldLocation::Cell, "gradient");
    requireField(face, mesh, FieldLocation::Face, "reconstruction");
    for (Index index = 0; index < mesh.faceCount(); ++index) {
        detail::fieldData(face)[index] = correctedFaceValue(cell, cell_gradient, index);
    }
}

void gradient(
    const ScalarField& scalar,
    VectorField& result,
    GradientMethod method)
{
    const Mesh& mesh = scalar.mesh();
    requireField(scalar, mesh, FieldLocation::Cell, "scalar");
    requireField(result, mesh, FieldLocation::Cell, "gradient");
    switch (method) {
        case GradientMethod::GreenGauss:
            greenGaussGradient(scalar, result);
            return;
        case GradientMethod::LeastSquares:
            leastSquaresGradient(scalar, result);
            return;
    }
    throw std::invalid_argument("unsupported gradient method");
}

void gradient(
    const VectorField& vector,
    TensorField& result,
    GradientMethod method)
{
    const Mesh& mesh = vector.mesh();
    requireField(vector, mesh, FieldLocation::Cell, "vector");
    requireField(result, mesh, FieldLocation::Cell, "vector gradient");
    switch (method) {
        case GradientMethod::GreenGauss:
            greenGaussGradient(vector, result);
            return;
        case GradientMethod::LeastSquares:
            leastSquaresGradient(vector, result);
            return;
    }
    throw std::invalid_argument("unsupported vector-gradient method");
}

double integratedNormalGradient(
    const ScalarField& scalar,
    const VectorField& scalar_gradient,
    Index face,
    DiffusionMethod method)
{
    const Mesh& mesh = scalar.mesh();
    requireField(scalar, mesh, FieldLocation::Cell, "scalar");
    requireField(scalar_gradient, mesh, FieldLocation::Cell, "scalar gradient");
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner.at(f);
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    if (neighbour != invalid_index) {
        const double orthogonal = detail::meshData(mesh).face_orthogonal_coefficients[f] *
            (detail::fieldData(scalar)[neighbour] - detail::fieldData(scalar)[owner]);
        if (!correctedDiffusion(method)) return orthogonal;
        const double weight = ownerWeight(mesh, face);
        const Vec3 face_gradient =
            weight * detail::fieldData(scalar_gradient)[owner] +
            (1.0 - weight) * detail::fieldData(scalar_gradient)[neighbour];
        return orthogonal + limitNonOrthogonalCorrection(
            orthogonal,
            dot(detail::meshData(mesh).face_non_orthogonal[f], face_gradient), method);
    }

    const auto& condition = scalar.boundary(detail::meshData(mesh).face_patch[f]);
    if (condition.type == BoundaryType::FixedGradient) {
        return condition.value * detail::meshData(mesh).face_areas[f];
    }
    if (condition.type != BoundaryType::FixedValue) {
        return 0.0;
    }
    const double orthogonal = detail::meshData(mesh).face_orthogonal_coefficients[f] *
        (condition.value - detail::fieldData(scalar)[owner]);
    return correctedDiffusion(method)
        ? orthogonal + limitNonOrthogonalCorrection(
              orthogonal,
              dot(detail::meshData(mesh).face_non_orthogonal[f], detail::fieldData(scalar_gradient)[owner]), method)
        : orthogonal;
}

Vec3 integratedNormalGradient(
    const VectorField& vector,
    const TensorField& vector_gradient,
    Index face,
    DiffusionMethod method)
{
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
                (detail::fieldData(vector)[neighbour][component] - detail::fieldData(vector)[owner][component]);
            result[component] = orthogonal;
            if (correctedDiffusion(method)) {
                const Vec3 face_gradient =
                    weight * detail::fieldData(vector_gradient)[owner][component] +
                    (1.0 - weight) * detail::fieldData(vector_gradient)[neighbour][component];
                result[component] += limitNonOrthogonalCorrection(
                    orthogonal,
                    dot(detail::meshData(mesh).face_non_orthogonal[f], face_gradient), method);
            }
        }
        return result;
    }

    const auto& condition = vector.boundary(detail::meshData(mesh).face_patch[f]);
    if (condition.type == BoundaryType::FixedGradient) {
        return detail::meshData(mesh).face_areas[f] * condition.value;
    }
    if (condition.type != BoundaryType::FixedValue &&
        condition.type != BoundaryType::Symmetry) {
        return {};
    }
    const Vec3 boundary_value = boundaryFaceValue(vector, face);
    const Vec3 normal = mesh.faceNormal(face);
    for (std::size_t component = 0; component < 3; ++component) {
        const double orthogonal = detail::meshData(mesh).face_orthogonal_coefficients[f] *
            (boundary_value[component] - detail::fieldData(vector)[owner][component]);
        result[component] = orthogonal;
        if (correctedDiffusion(method)) {
            Vec3 boundary_gradient = detail::fieldData(vector_gradient)[owner][component];
            if (condition.type == BoundaryType::Symmetry) {
                boundary_gradient -= dot(boundary_gradient, normal) * normal;
            }
            result[component] += limitNonOrthogonalCorrection(
                orthogonal,
                dot(detail::meshData(mesh).face_non_orthogonal[f], boundary_gradient), method);
        }
    }
    return result;
}


void flux(
    const VectorField& velocity,
    ScalarField& face_flux,
    InterpolationMethod method,
    GradientMethod gradient_method)
{
    const Mesh& mesh = velocity.mesh();
    requireField(face_flux, mesh, FieldLocation::Face, "flux");
    if (velocity.location() == FieldLocation::Face) {
        requireField(velocity, mesh, FieldLocation::Face, "face vector");
        for (Index face : detail::meshData(mesh).owned_faces)
            detail::fieldData(face_flux)[face] = dot(
                detail::fieldData(velocity)[face], mesh.faceAreaVector(face));
        return;
    }
    requireField(velocity, mesh, FieldLocation::Cell, "velocity");
    if (method != InterpolationMethod::Linear &&
        method != InterpolationMethod::Corrected) {
        throw std::invalid_argument("unsupported flux interpolation method");
    }
    std::optional<TensorField> velocity_gradient;
    if (method == InterpolationMethod::Corrected) {
        velocity_gradient.emplace(mesh, FieldLocation::Cell, "grad(" + velocity.name() + ')');
        gradient(velocity, *velocity_gradient, gradient_method);
    }
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        // 入口出口边界需要上一轮面通量判定流向；普通边界会忽略该参数。
        const double previous_flux = detail::fieldData(face_flux)[face];
        const Vec3 face_velocity = velocity_gradient
            ? correctedFaceValue(
                  velocity, *velocity_gradient, face, previous_flux)
            : interpolatedFaceValue(velocity, face, previous_flux);
        detail::fieldData(face_flux)[face] = dot(
            face_velocity, detail::meshData(mesh).face_area_vectors[f]);
    }
}

void diffusionFlux(
    const ScalarField& face_diffusivity,
    const ScalarField& scalar,
    const VectorField& scalar_gradient,
    ScalarField& face_flux,
    DiffusionMethod diffusion_method)
{
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
        detail::fieldData(face_flux)[face] = diffusivity * integratedNormalGradient(
            scalar, scalar_gradient, face, diffusion_method);
    }
}

void divergence(
    const VectorField& vector,
    ScalarField& result,
    InterpolationMethod method,
    GradientMethod gradient_method)
{
    const Mesh& mesh = vector.mesh();
    requireField(vector, mesh, FieldLocation::Cell, "vector");
    requireField(result, mesh, FieldLocation::Cell, "divergence");
    if (method != InterpolationMethod::Linear &&
        method != InterpolationMethod::Corrected) {
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
        const Vec3 face_value = vector_gradient
            ? correctedFaceValue(vector, *vector_gradient, face)
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
void convectionImpl(
    const ScalarField& face_flux,
    const Field<T>& transported,
    Field<T>& result,
    ConvectionMethod method,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method)
{
    const Mesh& mesh = transported.mesh();
    requireField(face_flux, mesh, FieldLocation::Face, "face flux");
    requireField(transported, mesh, FieldLocation::Cell, "transported");
    requireField(result, mesh, FieldLocation::Cell, "convection");
    if (method != ConvectionMethod::Upwind &&
        method != ConvectionMethod::LinearUpwind &&
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
        transported_gradient.emplace(
            mesh, FieldLocation::Cell, "grad(" + transported.name() + ')');
        gradient(transported, *transported_gradient, gradient_method);
    }
    result.fill(T{});
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const std::size_t index = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[index];
        const Index neighbour = detail::meshData(mesh).face_neighbour[index];
        const double flux_value = detail::fieldData(face_flux)[face];
        T face_value{};
        if (method == ConvectionMethod::Upwind ||
            method == ConvectionMethod::LinearUpwind) {
            if (neighbour != invalid_index) {
                const Index upwind = flux_value >= 0.0 ? owner : neighbour;
                face_value = detail::fieldData(transported)[upwind];
                if (method == ConvectionMethod::LinearUpwind) {
                    face_value += linearUpwindCorrection(
                        transported, *transported_gradient, upwind, face);
                }
            } else {
                if (flux_value >= 0.0) {
                    face_value = detail::fieldData(transported)[owner];
                    if (method == ConvectionMethod::LinearUpwind) {
                        face_value += linearUpwindCorrection(
                            transported, *transported_gradient, owner, face);
                    }
                } else {
                    face_value = boundaryFaceValue(transported, face, flux_value);
                }
            }
        } else if (transported_gradient) {
            face_value = correctedFaceValue(
                transported, *transported_gradient, face, flux_value);
        } else {
            face_value = interpolatedFaceValue(transported, face, flux_value);
        }
        addFaceDivergence(mesh, face, flux_value * face_value, result);
    }
}

void convection(
    const ScalarField& face_flux,
    const ScalarField& transported,
    ScalarField& result,
    ConvectionMethod method,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method)
{
    convectionImpl(
        face_flux, transported, result, method, interpolation_method, gradient_method);
}

void convection(
    const ScalarField& face_flux,
    const VectorField& transported,
    VectorField& result,
    ConvectionMethod method,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method)
{
    convectionImpl(
        face_flux, transported, result, method, interpolation_method, gradient_method);
}

void laplacianImpl(
    const ScalarField& scalar,
    ScalarField& result,
    double constant_diffusivity,
    const ScalarField* face_diffusivity,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
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
        cell_gradient.emplace(
            mesh, FieldLocation::Cell, "grad(" + scalar.name() + ')');
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
            integrated_flux = detail::meshData(mesh).face_orthogonal_coefficients[f] *
                (detail::fieldData(scalar)[neighbour] - detail::fieldData(scalar)[owner]);
            if (correctedDiffusion(diffusion_method)) {
                integrated_flux = integratedNormalGradient(
                    scalar, *cell_gradient, face, diffusion_method);
            }
        } else {
            const auto& condition = scalar.boundary(detail::meshData(mesh).face_patch[f]);
            if (condition.type == BoundaryType::FixedValue) {
                integrated_flux = detail::meshData(mesh).face_orthogonal_coefficients[f] *
                    (boundaryFaceValue(scalar, face, -1.0) - detail::fieldData(scalar)[owner]);
                if (correctedDiffusion(diffusion_method)) {
                    integrated_flux = integratedNormalGradient(
                        scalar, *cell_gradient, face, diffusion_method);
                }
            } else if (condition.type == BoundaryType::FixedGradient) {
                integrated_flux = condition.value * detail::meshData(mesh).face_areas[f];
            }
        }
        const double diffusivity = face_diffusivity == nullptr
            ? constant_diffusivity : detail::fieldData((*face_diffusivity))[face];
        if (!(diffusivity >= 0.0) || !std::isfinite(diffusivity)) {
            throw std::invalid_argument("face diffusivity must be non-negative and finite");
        }
        addFaceDivergence(mesh, face, diffusivity * integrated_flux, result);
    }
}

void laplacian(
    const ScalarField& scalar,
    ScalarField& result,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    laplacianImpl(scalar, result, 1.0, nullptr, gradient_method, diffusion_method);
}

void laplacian(
    double diffusivity,
    const ScalarField& scalar,
    ScalarField& result,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    laplacianImpl(scalar, result, diffusivity, nullptr, gradient_method, diffusion_method);
}

void laplacian(
    const ScalarField& face_diffusivity,
    const ScalarField& scalar,
    ScalarField& result,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    laplacianImpl(
        scalar, result, 1.0, &face_diffusivity, gradient_method, diffusion_method);
}

void addConvection(
    ScalarDiscreteEquation& equation,
    const ScalarField& face_flux,
    const ScalarField& transported,
    ConvectionMethod method,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method,
    double flux_scale)
{
    addConvectionImpl(
        equation, face_flux, transported, method,
        interpolation_method, gradient_method, flux_scale);
}

void addConvection(
    VectorDiscreteEquation& equation,
    const ScalarField& face_flux,
    const VectorField& transported,
    ConvectionMethod method,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method,
    double flux_scale)
{
    addConvectionImpl(
        equation, face_flux, transported, method,
        interpolation_method, gradient_method, flux_scale);
}

namespace {

struct FaceDiffusivity {
    double constant = 0.0;
    const ScalarField* field = nullptr;

    double operator[](Index face) const {
        return field == nullptr ? constant : detail::fieldData((*field))[face];
    }
};

void addScalarDiffusion(
    ScalarDiscreteEquation& equation,
    FaceDiffusivity diffusivity,
    const ScalarField& scalar,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    const Mesh& mesh = scalar.mesh();
    requireEquation(equation, mesh);
    requireField(scalar, mesh, FieldLocation::Cell, "scalar");

    std::optional<VectorField> cell_gradient;
    if (correctedDiffusion(diffusion_method)) {
        cell_gradient.emplace(
            mesh, FieldLocation::Cell, "grad(" + scalar.name() + ')');
        gradient(scalar, *cell_gradient, gradient_method);
    } else if (diffusion_method != DiffusionMethod::Orthogonal) {
        throw std::invalid_argument("unsupported diffusion method");
    }

    for (Index face : detail::meshData(mesh).owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const double gamma = diffusivity[face];
        if (gamma < 0.0 || !std::isfinite(gamma)) {
            throw std::invalid_argument("face diffusivity must be non-negative and finite");
        }
        if (gamma == 0.0) {
            continue;
        }
        const double coefficient =
            gamma * detail::meshData(mesh).face_orthogonal_coefficients[f];
        if (neighbour != invalid_index) {
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.diagonal[static_cast<std::size_t>(neighbour)] += coefficient;
            equation.upper[f] -= coefficient;
            equation.lower[f] -= coefficient;
            if (correctedDiffusion(diffusion_method)) {
                const double orthogonal =
                    detail::meshData(mesh).face_orthogonal_coefficients[f] *
                    (detail::fieldData(scalar)[neighbour] - detail::fieldData(scalar)[owner]);
                const double correction = gamma *
                    (integratedNormalGradient(
                         scalar, *cell_gradient, face, diffusion_method) -
                     orthogonal);
                equation.source[static_cast<std::size_t>(owner)] += correction;
                equation.source[static_cast<std::size_t>(neighbour)] -= correction;
            }
            continue;
        }

        const auto& condition = scalar.boundary(detail::meshData(mesh).face_patch[f]);
        if (condition.type == BoundaryType::FixedValue) {
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.source[static_cast<std::size_t>(owner)] +=
                coefficient * condition.value;
            if (correctedDiffusion(diffusion_method)) {
                const double orthogonal =
                    detail::meshData(mesh).face_orthogonal_coefficients[f] *
                    (condition.value - detail::fieldData(scalar)[owner]);
                equation.source[static_cast<std::size_t>(owner)] += gamma *
                    (integratedNormalGradient(
                         scalar, *cell_gradient, face, diffusion_method) -
                     orthogonal);
            }
        } else if (condition.type == BoundaryType::FixedGradient) {
            equation.source[static_cast<std::size_t>(owner)] +=
                gamma * condition.value * detail::meshData(mesh).face_areas[f];
        }
    }
}

void addVectorDiffusion(
    VectorDiscreteEquation& equation,
    FaceDiffusivity diffusivity,
    const VectorField& vector,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    const Mesh& mesh = vector.mesh();
    requireEquation(equation, mesh);
    requireField(vector, mesh, FieldLocation::Cell, "vector");
    std::optional<TensorField> cell_gradient;
    if (correctedDiffusion(diffusion_method)) {
        cell_gradient.emplace(
            mesh, FieldLocation::Cell, "grad(" + vector.name() + ')');
        gradient(vector, *cell_gradient, gradient_method);
    } else if (diffusion_method != DiffusionMethod::Orthogonal) {
        throw std::invalid_argument("unsupported vector diffusion method");
    }

    for (Index face : detail::meshData(mesh).owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const double gamma = diffusivity[face];
        if (gamma < 0.0 || !std::isfinite(gamma)) {
            throw std::invalid_argument("face diffusivity must be non-negative and finite");
        }
        if (gamma == 0.0) {
            continue;
        }
        const double coefficient =
            gamma * detail::meshData(mesh).face_orthogonal_coefficients[f];
        if (neighbour != invalid_index) {
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.diagonal[static_cast<std::size_t>(neighbour)] += coefficient;
            equation.upper[f] -= coefficient;
            equation.lower[f] -= coefficient;
            if (correctedDiffusion(diffusion_method)) {
                const Vec3 orthogonal =
                    detail::meshData(mesh).face_orthogonal_coefficients[f] *
                    (detail::fieldData(vector)[neighbour] - detail::fieldData(vector)[owner]);
                const Vec3 correction = gamma *
                    (integratedNormalGradient(
                         vector, *cell_gradient, face, diffusion_method) -
                     orthogonal);
                equation.source[static_cast<std::size_t>(owner)] += correction;
                equation.source[static_cast<std::size_t>(neighbour)] -= correction;
            }
            continue;
        }

        const auto& condition = vector.boundary(detail::meshData(mesh).face_patch[f]);
        if (condition.type == BoundaryType::FixedValue ||
            condition.type == BoundaryType::Symmetry) {
            const Vec3 boundary_value = boundaryFaceValue(vector, face, -1.0);
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.source[static_cast<std::size_t>(owner)] +=
                coefficient * boundary_value;
            if (correctedDiffusion(diffusion_method)) {
                const Vec3 orthogonal =
                    detail::meshData(mesh).face_orthogonal_coefficients[f] *
                    (boundary_value - detail::fieldData(vector)[owner]);
                const Vec3 correction = gamma *
                    (integratedNormalGradient(
                         vector, *cell_gradient, face, diffusion_method) -
                     orthogonal);
                equation.source[static_cast<std::size_t>(owner)] += correction;
            }
        } else if (condition.type == BoundaryType::FixedGradient) {
            equation.source[static_cast<std::size_t>(owner)] +=
                gamma * detail::meshData(mesh).face_areas[f] * condition.value;
        }
    }
}

}  // 匿名命名空间

void addDiffusion(
    ScalarDiscreteEquation& equation,
    double diffusivity,
    const ScalarField& scalar,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    if (diffusivity < 0.0 || !std::isfinite(diffusivity)) {
        throw std::invalid_argument("diffusivity must be non-negative and finite");
    }
    if (diffusivity != 0.0) {
        addScalarDiffusion(
            equation, {diffusivity, nullptr}, scalar,
            gradient_method, diffusion_method);
    }
}

void addDiffusion(
    ScalarDiscreteEquation& equation,
    const ScalarField& face_diffusivity,
    const ScalarField& scalar,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    requireField(
        face_diffusivity, scalar.mesh(), FieldLocation::Face,
        "face diffusivity");
    addScalarDiffusion(
        equation, {0.0, &face_diffusivity}, scalar,
        gradient_method, diffusion_method);
}

void addDiffusion(
    VectorDiscreteEquation& equation,
    double diffusivity,
    const VectorField& vector,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    if (diffusivity < 0.0 || !std::isfinite(diffusivity)) {
        throw std::invalid_argument("diffusivity must be non-negative and finite");
    }
    if (diffusivity != 0.0) {
        addVectorDiffusion(
            equation, {diffusivity, nullptr}, vector,
            gradient_method, diffusion_method);
    }
}

void addDiffusion(
    VectorDiscreteEquation& equation,
    const ScalarField& face_diffusivity,
    const VectorField& vector,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    requireField(
        face_diffusivity, vector.mesh(), FieldLocation::Face,
        "face diffusivity");
    addVectorDiffusion(
        equation, {0.0, &face_diffusivity}, vector,
        gradient_method, diffusion_method);
}

void addTimeDerivative(
    ScalarDiscreteEquation& equation,
    const ScalarField& previous,
    double dt,
    double density,
    TimeMethod method,
    const ScalarField* older)
{
    addTimeDerivativeImpl(equation, previous, dt, density, method, older);
}

void addTimeDerivative(
    ScalarDiscreteEquation& equation,
    const ScalarField& previous,
    double dt,
    const ScalarField& volumetric_capacity,
    TimeMethod method,
    const ScalarField* older)
{
    addFieldTimeDerivativeImpl(
        equation, previous, dt, volumetric_capacity, method, older);
}

void addTimeDerivative(
    VectorDiscreteEquation& equation,
    const VectorField& previous,
    double dt,
    const ScalarField& density,
    TimeMethod method,
    const VectorField* older)
{
    addFieldTimeDerivativeImpl(equation, previous, dt, density, method, older);
}

void addTimeDerivative(
    VectorDiscreteEquation& equation,
    const VectorField& previous,
    double dt,
    double density,
    TimeMethod method,
    const VectorField* older)
{
    addTimeDerivativeImpl(equation, previous, dt, density, method, older);
}

}  // babelsim 命名空间
