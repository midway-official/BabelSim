#include "babelsim/operators.h"

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
    return mesh.face_owner_weights[f];
}

template <typename T>
T internalFaceValue(const Field<T>& field, Index face) {
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = mesh.face_owner[f];
    const Index neighbour = mesh.face_neighbour[f];
    const double weight = ownerWeight(mesh, face);
    return weight * field[owner] + (1.0 - weight) * field[neighbour];
}

double correctedInternalFaceValue(
    const ScalarField& field,
    const VectorField& field_gradient,
    Index face)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = mesh.face_owner[f];
    const Index neighbour = mesh.face_neighbour[f];
    const double weight = ownerWeight(mesh, face);
    const Vec3 face_gradient =
        weight * field_gradient[owner] +
        (1.0 - weight) * field_gradient[neighbour];
    return internalFaceValue(field, face) +
        dot(mesh.face_skewness[f], face_gradient);
}

Vec3 correctedInternalFaceValue(
    const VectorField& field,
    const TensorField& field_gradient,
    Index face)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = mesh.face_owner[f];
    const Index neighbour = mesh.face_neighbour[f];
    const double weight = ownerWeight(mesh, face);
    Vec3 value = internalFaceValue(field, face);
    for (std::size_t component = 0; component < 3; ++component) {
        const Vec3 face_gradient =
            weight * field_gradient[owner][component] +
            (1.0 - weight) * field_gradient[neighbour][component];
        value[component] += dot(mesh.face_skewness[f], face_gradient);
    }
    return value;
}

double reconstructedBoundaryValue(
    const ScalarField& field,
    const VectorField& field_gradient,
    Index face,
    double outward_flux = 0.0)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = mesh.face_owner[f];
    const auto& condition = field.boundary(mesh.face_patch[f]);
    if (condition.type == BoundaryType::FixedValue ||
        (condition.type == BoundaryType::InletOutlet && outward_flux < 0.0)) {
        return condition.value;
    }
    const Vec3 normal = mesh.faceNormal(face);
    const double prescribed_normal_gradient =
        condition.type == BoundaryType::FixedGradient ? condition.value : 0.0;
    const Vec3 constrained_gradient = field_gradient[owner] +
        (prescribed_normal_gradient - dot(field_gradient[owner], normal)) * normal;
    return field[owner] + dot(
        constrained_gradient,
        mesh.face_centres[f] - mesh.cell_centres[static_cast<std::size_t>(owner)]);
}

Vec3 reconstructedBoundaryValue(
    const VectorField& field,
    const TensorField& field_gradient,
    Index face,
    double outward_flux = 0.0)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = mesh.face_owner[f];
    const auto& condition = field.boundary(mesh.face_patch[f]);
    if (condition.type == BoundaryType::FixedValue ||
        (condition.type == BoundaryType::InletOutlet && outward_flux < 0.0)) {
        return condition.value;
    }
    const Vec3 normal = mesh.faceNormal(face);
    const Vec3 offset =
        mesh.face_centres[f] - mesh.cell_centres[static_cast<std::size_t>(owner)];
    Vec3 value = field[owner];
    for (std::size_t component = 0; component < 3; ++component) {
        const double prescribed_normal_gradient =
            condition.type == BoundaryType::FixedGradient
            ? condition.value[component] : 0.0;
        const Vec3 constrained_gradient = field_gradient[owner][component] +
            (prescribed_normal_gradient -
             dot(field_gradient[owner][component], normal)) * normal;
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
        face[f] = interpolatedFaceValue(cell, f);
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
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        const double value = interpolatedFaceValue(scalar, face);
        const Vec3 contribution = value * mesh.face_area_vectors[f];
        result[owner] += contribution *
            mesh.cell_inverse_volumes[static_cast<std::size_t>(owner)];
        if (neighbour != invalid_index) {
            result[neighbour] -= contribution *
                mesh.cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
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
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        const double value = correctedFaceValue(scalar, initial_gradient, face);
        const Vec3 contribution = value * mesh.face_area_vectors[f];
        result[owner] += contribution *
            mesh.cell_inverse_volumes[static_cast<std::size_t>(owner)];
        if (neighbour != invalid_index) {
            result[neighbour] -= contribution *
                mesh.cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
        }
    }
}

void leastSquaresGradient(const ScalarField& scalar, VectorField& result) {
    const Mesh& mesh = scalar.mesh();
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
        Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
        const auto c = static_cast<std::size_t>(cell);
        for (Index face : mesh.cell_faces[c]) {
            const auto f = static_cast<std::size_t>(face);
            const Index owner = mesh.face_owner[f];
            const Index neighbour = mesh.face_neighbour[f];
            Vec3 delta{};
            double difference = 0.0;
            if (neighbour != invalid_index) {
                const Index other = owner == cell ? neighbour : owner;
                delta = mesh.cell_centres[static_cast<std::size_t>(other)] -
                    mesh.cell_centres[c];
                difference = scalar[other] - scalar[cell];
            } else {
                const auto type = scalar.boundary(mesh.face_patch[f]).type;
                const Vec3 offset = mesh.face_centres[f] - mesh.cell_centres[c];
                delta = type == BoundaryType::FixedValue
                    ? 2.0 * offset
                    : 2.0 * boundaryNormalDistance(mesh, face) *
                        mesh.faceNormal(face);
                difference = 2.0 * (boundaryFaceValue(scalar, face) - scalar[cell]);
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
        result[cell] = toVec3(value);
    }
}

void greenGaussGradient(const VectorField& vector, TensorField& result) {
    const Mesh& mesh = vector.mesh();
    result.fill({});
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        const Vec3 value = interpolatedFaceValue(vector, face);
        for (std::size_t component = 0; component < 3; ++component) {
            const Vec3 contribution = value[component] * mesh.face_area_vectors[f];
            result[owner][component] += contribution *
                mesh.cell_inverse_volumes[static_cast<std::size_t>(owner)];
            if (neighbour != invalid_index) {
                result[neighbour][component] -= contribution *
                    mesh.cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
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
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        const Vec3 value = correctedFaceValue(vector, initial_gradient, face);
        for (std::size_t component = 0; component < 3; ++component) {
            const Vec3 contribution = value[component] * mesh.face_area_vectors[f];
            result[owner][component] += contribution *
                mesh.cell_inverse_volumes[static_cast<std::size_t>(owner)];
            if (neighbour != invalid_index) {
                result[neighbour][component] -= contribution *
                    mesh.cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
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
        for (Index face : mesh.cell_faces[c]) {
            const auto f = static_cast<std::size_t>(face);
            const Index owner = mesh.face_owner[f];
            const Index neighbour = mesh.face_neighbour[f];
            Vec3 delta{};
            Vec3 difference{};
            if (neighbour != invalid_index) {
                const Index other = owner == cell ? neighbour : owner;
                delta = mesh.cell_centres[static_cast<std::size_t>(other)] -
                    mesh.cell_centres[c];
                difference = vector[other] - vector[cell];
            } else {
                const auto type = vector.boundary(mesh.face_patch[f]).type;
                const Vec3 offset = mesh.face_centres[f] - mesh.cell_centres[c];
                delta = type == BoundaryType::FixedValue
                    ? 2.0 * offset
                    : 2.0 * boundaryNormalDistance(mesh, face) *
                        mesh.faceNormal(face);
                difference = 2.0 * (boundaryFaceValue(vector, face) - vector[cell]);
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
            result[cell][component] = toVec3(
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
    const Index owner = mesh.face_owner[f];
    const Index neighbour = mesh.face_neighbour[f];
    result[owner] += integrated_flux *
        mesh.cell_inverse_volumes[static_cast<std::size_t>(owner)];
    if (neighbour != invalid_index) {
        result[neighbour] -= integrated_flux *
            mesh.cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
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
void requireEquation(const Equation<T>& equation, const Mesh& mesh) {
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
    Equation<T>& equation,
    const ScalarField& face_flux,
    const Field<T>& transported,
    ConvectionMethod method,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method)
{
    const Mesh& mesh = transported.mesh();
    requireEquation(equation, mesh);
    requireField(face_flux, mesh, FieldLocation::Face, "face flux");
    requireField(transported, mesh, FieldLocation::Cell, "transported");

    if (interpolation_method != InterpolationMethod::Linear &&
        interpolation_method != InterpolationMethod::Corrected) {
        throw std::invalid_argument("unsupported convection interpolation method");
    }
    using GradientField = typename GradientFieldType<T>::Type;
    std::optional<GradientField> transported_gradient;
    if (method == ConvectionMethod::Central &&
        interpolation_method == InterpolationMethod::Corrected) {
        transported_gradient.emplace(
            mesh, FieldLocation::Cell, "grad(" + transported.name() + ')');
        gradient(transported, *transported_gradient, gradient_method);
    }

    for (Index face : mesh.owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        const double F = face_flux[face];
        if (neighbour != invalid_index) {
            if (method == ConvectionMethod::Upwind) {
                equation.diagonal[static_cast<std::size_t>(owner)] += std::max(F, 0.0);
                equation.upper[f] += std::min(F, 0.0);
                equation.diagonal[static_cast<std::size_t>(neighbour)] +=
                    std::max(-F, 0.0);
                equation.lower[f] += std::min(-F, 0.0);
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

        const auto& condition = transported.boundary(mesh.face_patch[f]);
        if (method == ConvectionMethod::Upwind) {
            if (F >= 0.0) {
                equation.diagonal[static_cast<std::size_t>(owner)] += F;
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
    Equation<T>& equation,
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

    for (Index cell : mesh.owned_cells) {
        const auto c = static_cast<std::size_t>(cell);
        const double coefficient = density * mesh.cell_volumes[c] / dt;
        if (method == TimeMethod::Euler) {
            equation.diagonal[c] += coefficient;
            equation.source[c] += coefficient * previous[cell];
        } else {
            equation.diagonal[c] += 1.5 * coefficient;
            equation.source[c] += coefficient *
                (2.0 * previous[cell] - 0.5 * (*older)[cell]);
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
        face[index] = correctedFaceValue(cell, cell_gradient, index);
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
        face[index] = correctedFaceValue(cell, cell_gradient, index);
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
        face[index] = correctedFaceValue(cell, cell_gradient, index);
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
        face[index] = correctedFaceValue(cell, cell_gradient, index);
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
    const Index owner = mesh.face_owner.at(f);
    const Index neighbour = mesh.face_neighbour[f];
    if (neighbour != invalid_index) {
        const double orthogonal = mesh.face_orthogonal_coefficients[f] *
            (scalar[neighbour] - scalar[owner]);
        if (!correctedDiffusion(method)) return orthogonal;
        const double weight = ownerWeight(mesh, face);
        const Vec3 face_gradient =
            weight * scalar_gradient[owner] +
            (1.0 - weight) * scalar_gradient[neighbour];
        return orthogonal + limitNonOrthogonalCorrection(
            orthogonal,
            dot(mesh.face_non_orthogonal[f], face_gradient), method);
    }

    const auto& condition = scalar.boundary(mesh.face_patch[f]);
    if (condition.type == BoundaryType::FixedGradient) {
        return condition.value * mesh.face_areas[f];
    }
    if (condition.type != BoundaryType::FixedValue) {
        return 0.0;
    }
    const double orthogonal = mesh.face_orthogonal_coefficients[f] *
        (condition.value - scalar[owner]);
    return correctedDiffusion(method)
        ? orthogonal + limitNonOrthogonalCorrection(
              orthogonal,
              dot(mesh.face_non_orthogonal[f], scalar_gradient[owner]), method)
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
    const Index owner = mesh.face_owner.at(f);
    const Index neighbour = mesh.face_neighbour[f];
    Vec3 result{};
    if (neighbour != invalid_index) {
        const double weight = ownerWeight(mesh, face);
        for (std::size_t component = 0; component < 3; ++component) {
            const double orthogonal = mesh.face_orthogonal_coefficients[f] *
                (vector[neighbour][component] - vector[owner][component]);
            result[component] = orthogonal;
            if (correctedDiffusion(method)) {
                const Vec3 face_gradient =
                    weight * vector_gradient[owner][component] +
                    (1.0 - weight) * vector_gradient[neighbour][component];
                result[component] += limitNonOrthogonalCorrection(
                    orthogonal,
                    dot(mesh.face_non_orthogonal[f], face_gradient), method);
            }
        }
        return result;
    }

    const auto& condition = vector.boundary(mesh.face_patch[f]);
    if (condition.type == BoundaryType::FixedGradient) {
        return mesh.face_areas[f] * condition.value;
    }
    if (condition.type != BoundaryType::FixedValue &&
        condition.type != BoundaryType::Symmetry) {
        return {};
    }
    const Vec3 boundary_value = boundaryFaceValue(vector, face);
    const Vec3 normal = mesh.faceNormal(face);
    for (std::size_t component = 0; component < 3; ++component) {
        const double orthogonal = mesh.face_orthogonal_coefficients[f] *
            (boundary_value[component] - vector[owner][component]);
        result[component] = orthogonal;
        if (correctedDiffusion(method)) {
            Vec3 boundary_gradient = vector_gradient[owner][component];
            if (condition.type == BoundaryType::Symmetry) {
                boundary_gradient -= dot(boundary_gradient, normal) * normal;
            }
            result[component] += limitNonOrthogonalCorrection(
                orthogonal,
                dot(mesh.face_non_orthogonal[f], boundary_gradient), method);
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
    requireField(velocity, mesh, FieldLocation::Cell, "velocity");
    requireField(face_flux, mesh, FieldLocation::Face, "flux");
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
        const double previous_flux = face_flux[face];
        const Vec3 face_velocity = velocity_gradient
            ? correctedFaceValue(
                  velocity, *velocity_gradient, face, previous_flux)
            : interpolatedFaceValue(velocity, face, previous_flux);
        face_flux[face] = dot(
            face_velocity, mesh.face_area_vectors[f]);
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
        const double integrated_flux = dot(face_value, mesh.face_area_vectors[f]);
        addFaceDivergence(mesh, face, integrated_flux, result);
    }
}

void divergence(const ScalarField& face_flux, ScalarField& result) {
    const Mesh& mesh = face_flux.mesh();
    requireField(face_flux, mesh, FieldLocation::Face, "face flux");
    requireField(result, mesh, FieldLocation::Cell, "divergence");
    result.fill(0.0);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        addFaceDivergence(mesh, face, face_flux[face], result);
    }
}

void laplacian(
    const ScalarField& scalar,
    ScalarField& result,
    GradientMethod gradient_method,
    DiffusionMethod diffusion_method)
{
    const Mesh& mesh = scalar.mesh();
    requireField(scalar, mesh, FieldLocation::Cell, "scalar");
    requireField(result, mesh, FieldLocation::Cell, "laplacian");
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
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        double integrated_flux = 0.0;
        if (neighbour != invalid_index) {
            integrated_flux = mesh.face_orthogonal_coefficients[f] *
                (scalar[neighbour] - scalar[owner]);
            if (correctedDiffusion(diffusion_method)) {
                integrated_flux = integratedNormalGradient(
                    scalar, *cell_gradient, face, diffusion_method);
            }
        } else {
            const auto& condition = scalar.boundary(mesh.face_patch[f]);
            if (condition.type == BoundaryType::FixedValue) {
                integrated_flux = mesh.face_orthogonal_coefficients[f] *
                    (boundaryFaceValue(scalar, face, -1.0) - scalar[owner]);
                if (correctedDiffusion(diffusion_method)) {
                    integrated_flux = integratedNormalGradient(
                        scalar, *cell_gradient, face, diffusion_method);
                }
            } else if (condition.type == BoundaryType::FixedGradient) {
                integrated_flux = condition.value * mesh.face_areas[f];
            }
        }
        addFaceDivergence(mesh, face, integrated_flux, result);
    }
}

void addConvection(
    ScalarEquation& equation,
    const ScalarField& face_flux,
    const ScalarField& transported,
    ConvectionMethod method,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method)
{
    addConvectionImpl(
        equation, face_flux, transported, method,
        interpolation_method, gradient_method);
}

void addConvection(
    VectorEquation& equation,
    const ScalarField& face_flux,
    const VectorField& transported,
    ConvectionMethod method,
    InterpolationMethod interpolation_method,
    GradientMethod gradient_method)
{
    addConvectionImpl(
        equation, face_flux, transported, method,
        interpolation_method, gradient_method);
}

namespace {

struct FaceDiffusivity {
    double constant = 0.0;
    const ScalarField* field = nullptr;

    double operator[](Index face) const {
        return field == nullptr ? constant : (*field)[face];
    }
};

void addScalarDiffusion(
    ScalarEquation& equation,
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

    for (Index face : mesh.owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        const double gamma = diffusivity[face];
        if (gamma < 0.0 || !std::isfinite(gamma)) {
            throw std::invalid_argument("face diffusivity must be non-negative and finite");
        }
        if (gamma == 0.0) {
            continue;
        }
        const double coefficient =
            gamma * mesh.face_orthogonal_coefficients[f];
        if (neighbour != invalid_index) {
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.diagonal[static_cast<std::size_t>(neighbour)] += coefficient;
            equation.upper[f] -= coefficient;
            equation.lower[f] -= coefficient;
            if (correctedDiffusion(diffusion_method)) {
                const double orthogonal =
                    mesh.face_orthogonal_coefficients[f] *
                    (scalar[neighbour] - scalar[owner]);
                const double correction = gamma *
                    (integratedNormalGradient(
                         scalar, *cell_gradient, face, diffusion_method) -
                     orthogonal);
                equation.source[static_cast<std::size_t>(owner)] += correction;
                equation.source[static_cast<std::size_t>(neighbour)] -= correction;
            }
            continue;
        }

        const auto& condition = scalar.boundary(mesh.face_patch[f]);
        if (condition.type == BoundaryType::FixedValue) {
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.source[static_cast<std::size_t>(owner)] +=
                coefficient * condition.value;
            if (correctedDiffusion(diffusion_method)) {
                const double orthogonal =
                    mesh.face_orthogonal_coefficients[f] *
                    (condition.value - scalar[owner]);
                equation.source[static_cast<std::size_t>(owner)] += gamma *
                    (integratedNormalGradient(
                         scalar, *cell_gradient, face, diffusion_method) -
                     orthogonal);
            }
        } else if (condition.type == BoundaryType::FixedGradient) {
            equation.source[static_cast<std::size_t>(owner)] +=
                gamma * condition.value * mesh.face_areas[f];
        }
    }
}

void addVectorDiffusion(
    VectorEquation& equation,
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

    for (Index face : mesh.owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = mesh.face_owner[f];
        const Index neighbour = mesh.face_neighbour[f];
        const double gamma = diffusivity[face];
        if (gamma < 0.0 || !std::isfinite(gamma)) {
            throw std::invalid_argument("face diffusivity must be non-negative and finite");
        }
        if (gamma == 0.0) {
            continue;
        }
        const double coefficient =
            gamma * mesh.face_orthogonal_coefficients[f];
        if (neighbour != invalid_index) {
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.diagonal[static_cast<std::size_t>(neighbour)] += coefficient;
            equation.upper[f] -= coefficient;
            equation.lower[f] -= coefficient;
            if (correctedDiffusion(diffusion_method)) {
                const Vec3 orthogonal =
                    mesh.face_orthogonal_coefficients[f] *
                    (vector[neighbour] - vector[owner]);
                const Vec3 correction = gamma *
                    (integratedNormalGradient(
                         vector, *cell_gradient, face, diffusion_method) -
                     orthogonal);
                equation.source[static_cast<std::size_t>(owner)] += correction;
                equation.source[static_cast<std::size_t>(neighbour)] -= correction;
            }
            continue;
        }

        const auto& condition = vector.boundary(mesh.face_patch[f]);
        if (condition.type == BoundaryType::FixedValue ||
            condition.type == BoundaryType::Symmetry) {
            const Vec3 boundary_value = boundaryFaceValue(vector, face, -1.0);
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.source[static_cast<std::size_t>(owner)] +=
                coefficient * boundary_value;
            if (correctedDiffusion(diffusion_method)) {
                const Vec3 orthogonal =
                    mesh.face_orthogonal_coefficients[f] *
                    (boundary_value - vector[owner]);
                const Vec3 correction = gamma *
                    (integratedNormalGradient(
                         vector, *cell_gradient, face, diffusion_method) -
                     orthogonal);
                equation.source[static_cast<std::size_t>(owner)] += correction;
            }
        } else if (condition.type == BoundaryType::FixedGradient) {
            equation.source[static_cast<std::size_t>(owner)] +=
                gamma * mesh.face_areas[f] * condition.value;
        }
    }
}

}  // 匿名命名空间

void addDiffusion(
    ScalarEquation& equation,
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
    ScalarEquation& equation,
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
    VectorEquation& equation,
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
    VectorEquation& equation,
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
    ScalarEquation& equation,
    const ScalarField& previous,
    double dt,
    double density,
    TimeMethod method,
    const ScalarField* older)
{
    addTimeDerivativeImpl(equation, previous, dt, density, method, older);
}

void addTimeDerivative(
    VectorEquation& equation,
    const VectorField& previous,
    double dt,
    double density,
    TimeMethod method,
    const VectorField* older)
{
    addTimeDerivativeImpl(equation, previous, dt, density, method, older);
}

}  // babelsim 命名空间
