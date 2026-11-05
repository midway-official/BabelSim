#include "babelsim/operators.h"
#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "operator_support.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace babelsim::detail {
struct LeastSquaresGeometry {
    std::vector<Vec3> weighted_displacements;
    std::vector<bool> skew_boundaries;
    std::vector<Eigen::LDLT<Eigen::Matrix3d>> factors;
};
} // namespace babelsim::detail

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

Vec3 boundaryDisplacement(const Mesh& mesh, Index cell, Index face, bool fixedValue) {
    return fixedValue ? 2.0 * (mesh.faceCentre(face) - mesh.cellCentre(cell))
                      : 2.0 * boundaryNormalDistance(mesh, face) * mesh.faceNormal(face);
}

double sampleWeight(const Vec3& displacement) {
    const double distanceSquared = squaredNorm(displacement);
    if (!(distanceSquared > 0.0))
        throw std::runtime_error("least-squares neighbour distance is invalid");
    return 1.0 / distanceSquared;
}

const detail::LeastSquaresGeometry& leastSquaresGeometry(const Mesh& mesh) {
    auto& cache = *detail::meshData(mesh).least_squares_cache;
    std::call_once(cache.initialized, [&] {
        auto geometry = std::make_shared<detail::LeastSquaresGeometry>();
        geometry->weighted_displacements.resize(mesh.faceCount());
        geometry->skew_boundaries.resize(mesh.faceCount());
        geometry->factors.resize(mesh.cellCount());

        std::vector<Vec3> displacements(mesh.faceCount());
        std::vector<double> weights(mesh.faceCount());
        for (Index face = 0; face < mesh.faceCount(); ++face) {
            const Index owner = mesh.owner(face), neighbour = mesh.neighbour(face);
            auto& displacement = displacements[face];
            if (neighbour == invalid_index) {
                displacement = boundaryDisplacement(mesh, owner, face, false);
                const auto fixed = boundaryDisplacement(mesh, owner, face, true);
                geometry->skew_boundaries[face] = fixed.x != displacement.x ||
                                                  fixed.y != displacement.y ||
                                                  fixed.z != displacement.z;
            } else {
                displacement = mesh.cellCentre(neighbour) - mesh.cellCentre(owner);
            }
            weights[face] = sampleWeight(displacement);
            geometry->weighted_displacements[face] = weights[face] * displacement;
        }
        for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
            Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
            for (Index face : mesh.cellFaces(cell)) {
                const auto displacement = toEigen(displacements[face]);
                normal.noalias() += weights[face] * displacement * displacement.transpose();
            }
            geometry->factors[cell].compute(normal);
        }
        cache.geometry = std::move(geometry);
    });
    return *cache.geometry;
}

void addLeastSquaresRhs(Eigen::Vector3d& rhs, const Eigen::Vector3d& weighted, double difference) {
    rhs.noalias() += weighted * difference;
}

void addLeastSquaresRhs(Eigen::Matrix3d& rhs, const Eigen::Vector3d& weighted,
                        const Vec3& difference) {
    for (std::size_t component = 0; component < 3; ++component)
        rhs.col(static_cast<Eigen::Index>(component)).noalias() += weighted * difference[component];
}

template <class Rhs>
Rhs solveLeastSquares(const Eigen::LDLT<Eigen::Matrix3d>& factor, const Rhs& rhs) {
    if (factor.info() != Eigen::Success)
        throw std::runtime_error("least-squares gradient matrix is singular");
    const Rhs value = factor.solve(rhs);
    if (!value.allFinite())
        throw std::runtime_error("least-squares gradient solve failed");
    return value;
}

template <class T, class Rhs>
Rhs leastSquaresWithBoundarySamples(const Field<T>& field, Index cell) {
    const Mesh& mesh = field.mesh();
    const auto* values = detail::fieldData(field);
    Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
    Rhs rhs = Rhs::Zero();
    for (Index face : mesh.cellFaces(cell)) {
        Vec3 displacement;
        T difference;
        if (mesh.boundaryFace(face)) {
            const bool fixed =
                detail::FieldAccess::condition(field, face).type == BoundaryType::FixedValue;
            displacement = boundaryDisplacement(mesh, cell, face, fixed);
            difference = 2.0 * (boundaryFaceValue(field, face) - values[cell]);
        } else {
            const Index other = mesh.owner(face) == cell ? mesh.neighbour(face) : mesh.owner(face);
            displacement = mesh.cellCentre(other) - mesh.cellCentre(cell);
            difference = values[other] - values[cell];
        }
        const auto direction = toEigen(displacement);
        const double weight = sampleWeight(displacement);
        normal.noalias() += weight * direction * direction.transpose();
        addLeastSquaresRhs(rhs, weight * direction, difference);
    }
    return solveLeastSquares(Eigen::LDLT<Eigen::Matrix3d>(normal), rhs);
}

template <class T, class Gradient>
void leastSquaresGradient(const Field<T>& field, Field<Gradient>& result, bool ownedOnly = false) {
    using Rhs = std::conditional_t<std::is_same_v<T, double>, Eigen::Vector3d, Eigen::Matrix3d>;
    const Mesh& mesh = field.mesh();
    if (mesh.cellCount() <= 0)
        return;
    const auto& geometry = leastSquaresGeometry(mesh);
    const auto& meshData = detail::meshData(mesh);
    const auto* values = detail::fieldData(field);
    auto* gradients = detail::fieldData(result);
    const auto reconstructCell = [&](Index cell) {
        Rhs rhs = Rhs::Zero();
        bool needsBoundarySamples = false;
        for (Index face : mesh.cellFaces(cell)) {
            const Index owner = meshData.face_owner[face],
                        neighbour = meshData.face_neighbour[face];
            const auto weighted = toEigen(geometry.weighted_displacements[face]);
            if (neighbour == invalid_index) {
                // A skew fixed-value face uses a different mirror point. Mixed
                // boundaries can switch at runtime, so use their current type.
                if (geometry.skew_boundaries[face] &&
                    detail::FieldAccess::condition(field, face).type == BoundaryType::FixedValue) {
                    needsBoundarySamples = true;
                    break;
                }
                addLeastSquaresRhs(rhs, weighted,
                                   2.0 * (boundaryFaceValue(field, face) - values[cell]));
            } else if (owner == cell) {
                addLeastSquaresRhs(rhs, weighted, values[neighbour] - values[cell]);
            } else {
                addLeastSquaresRhs(rhs, -weighted, values[owner] - values[cell]);
            }
        }
        const auto value = needsBoundarySamples
                               ? leastSquaresWithBoundarySamples<T, Rhs>(field, cell)
                               : solveLeastSquares(geometry.factors[cell], rhs);
        if constexpr (std::is_same_v<T, double>) {
            gradients[cell] = toVec3(value);
        } else {
            for (std::size_t component = 0; component < 3; ++component)
                gradients[cell][component] =
                    toVec3(value.col(static_cast<Eigen::Index>(component)));
        }
    };
    if (ownedOnly) {
        for (Index cell : meshData.owned_cells)
            reconstructCell(cell);
    } else {
        for (Index cell = 0; cell < mesh.cellCount(); ++cell)
            reconstructCell(cell);
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

} // namespace

namespace detail::operators {
template <class T, class Gradient>
void gradientForExchangeImpl(const Field<T>& field, Field<Gradient>& result,
                             GradientMethod method) {
    if (method != GradientMethod::LeastSquares) {
        gradient(field, result, method);
        return;
    }
    requireField(field, field.mesh(), FieldLocation::Cell, "gradient input");
    requireField(result, field.mesh(), FieldLocation::Cell, "gradient result");
    leastSquaresGradient(field, result, true);
}

void gradientForExchange(const ScalarField& field, VectorField& result, GradientMethod method) {
    gradientForExchangeImpl(field, result, method);
}

void gradientForExchange(const VectorField& field, TensorField& result, GradientMethod method) {
    gradientForExchangeImpl(field, result, method);
}
} // namespace detail::operators

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
