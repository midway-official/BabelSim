#include "babelsim/operators.h"
#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "operator_support.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <cmath>
#include <stdexcept>

namespace babelsim {
using namespace detail::operators;
namespace {

Vec3 toVec3(const Eigen::Vector3d& value) {
    return {value.x(), value.y(), value.z()};
}

Eigen::Vector3d toEigen(const Vec3& value) {
    return {value.x, value.y, value.z};
}

void greenGaussGradient(const ScalarField& scalar, VectorField& result) {
    const Mesh& mesh = scalar.mesh();
    result.fill({});
    Vec3* result_values = detail::fieldData(result);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const double value = interpolatedFaceValue(scalar, face);
        const Vec3 contribution = value * detail::meshData(mesh).face_area_vectors[f];
        result_values[owner] +=
            contribution *
            detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(owner)];
        if (neighbour != invalid_index) {
            result_values[neighbour] -=
                contribution *
                detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
        }
    }

    // 用初始高斯梯度把交点面值重构到真实面中心。一次显式修正已能恢复光滑网格上的
    // 二阶面值。单独梯度需要两层邻居，消费 ghost 梯度的复合算子需要三层。
    if (mesh.orthogonalGeometry()) {
        return;
    }
    // Reconstruction only reads cell gradients; avoid copying boundary
    // traces and flux metadata of a complete calculated Field.
    const std::vector<Vec3> initial_gradient = detail::fieldValues(result);
    result.fill({});
    result_values = detail::fieldData(result);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const double value = correctedFaceValue(scalar, initial_gradient.data(), face);
        const Vec3 contribution = value * detail::meshData(mesh).face_area_vectors[f];
        result_values[owner] +=
            contribution *
            detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(owner)];
        if (neighbour != invalid_index) {
            result_values[neighbour] -=
                contribution *
                detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
        }
    }
}

template <typename T>
void leastSquaresSample(const Field<T>& field, Index cell, Index face, Vec3& delta, T& difference) {
    const Mesh& mesh = field.mesh();
    const auto c = static_cast<std::size_t>(cell);
    const auto f = static_cast<std::size_t>(face);
    const Index owner = detail::meshData(mesh).face_owner[f];
    const Index neighbour = detail::meshData(mesh).face_neighbour[f];
    if (neighbour != invalid_index) {
        const Index other = owner == cell ? neighbour : owner;
        delta = detail::meshData(mesh).cell_centres[static_cast<std::size_t>(other)] -
                detail::meshData(mesh).cell_centres[c];
        difference = detail::fieldData(field)[other] - detail::fieldData(field)[cell];
    } else {
        // Mirror the boundary sample. Fixed values use the centre-to-face
        // displacement; normal-gradient conditions use its normal projection.
        const auto type = detail::FieldAccess::condition(field, face).type;
        const Vec3 offset =
            detail::meshData(mesh).face_centres[f] - detail::meshData(mesh).cell_centres[c];
        delta = type == BoundaryType::FixedValue
                    ? 2.0 * offset
                    : 2.0 * boundaryNormalDistance(mesh, face) * mesh.faceNormal(face);
        difference = 2.0 * (boundaryFaceValue(field, face) - detail::fieldData(field)[cell]);
    }
}

void leastSquaresGradient(const ScalarField& scalar, VectorField& result) {
    const Mesh& mesh = scalar.mesh();
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
        Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
        for (Index face : mesh.cellFaces(cell)) {
            Vec3 delta{};
            double difference = 0.0;
            leastSquaresSample(scalar, cell, face, delta, difference);
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
    Tensor3* result_values = detail::fieldData(result);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const Vec3 value = interpolatedFaceValue(vector, face);
        for (std::size_t component = 0; component < 3; ++component) {
            const Vec3 contribution =
                value[component] * detail::meshData(mesh).face_area_vectors[f];
            result_values[owner][component] +=
                contribution *
                detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(owner)];
            if (neighbour != invalid_index) {
                result_values[neighbour][component] -=
                    contribution * detail::meshData(mesh)
                                       .cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
            }
        }
    }
    if (mesh.orthogonalGeometry()) {
        return;
    }
    const std::vector<Tensor3> initial_gradient = detail::fieldValues(result);
    result.fill({});
    result_values = detail::fieldData(result);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const Vec3 value = correctedFaceValue(vector, initial_gradient.data(), face);
        for (std::size_t component = 0; component < 3; ++component) {
            const Vec3 contribution =
                value[component] * detail::meshData(mesh).face_area_vectors[f];
            result_values[owner][component] +=
                contribution *
                detail::meshData(mesh).cell_inverse_volumes[static_cast<std::size_t>(owner)];
            if (neighbour != invalid_index) {
                result_values[neighbour][component] -=
                    contribution * detail::meshData(mesh)
                                       .cell_inverse_volumes[static_cast<std::size_t>(neighbour)];
            }
        }
    }
}

void leastSquaresGradient(const VectorField& vector, TensorField& result) {
    const Mesh& mesh = vector.mesh();
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d rhs = Eigen::Matrix3d::Zero();
        for (Index face : mesh.cellFaces(cell)) {
            Vec3 delta{};
            Vec3 difference{};
            leastSquaresSample(vector, cell, face, delta, difference);
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
            detail::fieldData(result)[cell][component] =
                toVec3(value.col(static_cast<Eigen::Index>(component)));
        }
    }
}

} // namespace

void gradient(const ScalarField& scalar, VectorField& result, GradientMethod method) {
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

void gradient(const VectorField& vector, TensorField& result, GradientMethod method) {
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

} // namespace babelsim
