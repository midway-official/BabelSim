#pragma once

#include "babelsim/mesh.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace babelsim {

enum class FieldLocation {
    Cell,
    Face,
    Vertex,
};

enum class BoundaryType {
    FixedValue,
    Dirichlet = FixedValue,
    FixedGradient,
    Neumann = FixedGradient,
    ZeroGradient,
    InletOutlet,
    Symmetry,
    Mirror = Symmetry,
};

template <typename T>
struct BoundaryCondition {
    BoundaryType type = BoundaryType::ZeroGradient;
    // 固定值、固定外法向梯度，或入口值。
    T value{};

    static BoundaryCondition fixedValue(T value) {
        return {BoundaryType::FixedValue, std::move(value)};
    }
    static BoundaryCondition fixedGradient(T gradient) {
        return {BoundaryType::FixedGradient, std::move(gradient)};
    }
    static BoundaryCondition zeroGradient() {
        return {BoundaryType::ZeroGradient, T{}};
    }
    static BoundaryCondition symmetry() {
        return {BoundaryType::Symmetry, T{}};
    }
    static BoundaryCondition inletOutlet(T inlet_value) {
        return {BoundaryType::InletOutlet, std::move(inlet_value)};
    }
};

template <typename T>
class Field {
public:
    Field(
        const Mesh& mesh,
        FieldLocation location,
        std::string name = {},
        T initial = T{})
        : mesh_(&mesh),
          location_(location),
          name_(std::move(name)),
          values_(entityCount(mesh, location), std::move(initial))
    {
        if (location == FieldLocation::Cell) {
            boundaries_.resize(mesh.patches.size());
        }
    }

    const Mesh& mesh() const { return *mesh_; }
    FieldLocation location() const { return location_; }
    const std::string& name() const { return name_; }
    std::size_t size() const { return values_.size(); }

    T* data() { return values_.data(); }
    const T* data() const { return values_.data(); }
    T& operator[](Index index) { return values_[static_cast<std::size_t>(index)]; }
    const T& operator[](Index index) const {
        return values_[static_cast<std::size_t>(index)];
    }
    std::vector<T>& values() { return values_; }
    const std::vector<T>& values() const { return values_; }
    void fill(const T& value) { std::fill(values_.begin(), values_.end(), value); }

    void setBoundary(Index patch, BoundaryCondition<T> condition) {
        requireCellBoundary(patch);
        boundaries_[static_cast<std::size_t>(patch)] = std::move(condition);
    }
    BoundaryCondition<T>& boundary(Index patch) {
        requireCellBoundary(patch);
        return boundaries_[static_cast<std::size_t>(patch)];
    }
    const BoundaryCondition<T>& boundary(Index patch) const {
        requireCellBoundary(patch);
        return boundaries_[static_cast<std::size_t>(patch)];
    }

private:
    static std::size_t entityCount(const Mesh& mesh, FieldLocation location) {
        switch (location) {
            case FieldLocation::Cell:
                return static_cast<std::size_t>(mesh.cellCount());
            case FieldLocation::Face:
                return static_cast<std::size_t>(mesh.faceCount());
            case FieldLocation::Vertex:
                return static_cast<std::size_t>(mesh.vertexCount());
        }
        throw std::invalid_argument("unknown field location");
    }

    void requireCellBoundary(Index patch) const {
        if (location_ != FieldLocation::Cell || patch < 0 ||
            static_cast<std::size_t>(patch) >= boundaries_.size()) {
            throw std::out_of_range("field boundary patch is invalid");
        }
    }

    const Mesh* mesh_;
    FieldLocation location_;
    std::string name_;
    std::vector<T> values_;
    std::vector<BoundaryCondition<T>> boundaries_;
};

using ScalarField = Field<double>;
using VectorField = Field<Vec3>;
using TensorField = Field<Tensor3>;

inline double boundaryNormalDistance(const Mesh& mesh, Index face) {
    const auto f = static_cast<std::size_t>(face);
    const Index owner = mesh.face_owner[f];
    return dot(
        mesh.face_centres[f] - mesh.cell_centres[static_cast<std::size_t>(owner)],
        mesh.faceNormal(face));
}

inline double boundaryFaceValue(
    const ScalarField& field,
    Index face,
    double outward_flux = 0.0)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = mesh.face_owner[f];
    const Index patch = mesh.face_patch[f];
    const auto& condition = field.boundary(patch);
    const double owner_value = field[owner];
    const double distance = boundaryNormalDistance(mesh, face);
    switch (condition.type) {
        case BoundaryType::FixedValue:
            return condition.value;
        case BoundaryType::FixedGradient:
            return owner_value + condition.value * distance;
        case BoundaryType::ZeroGradient:
        case BoundaryType::Symmetry:
            return owner_value;
        case BoundaryType::InletOutlet:
            return outward_flux >= 0.0 ? owner_value : condition.value;
    }
    throw std::invalid_argument("unknown scalar boundary condition");
}

inline Vec3 boundaryFaceValue(
    const VectorField& field,
    Index face,
    double outward_flux = 0.0)
{
    const Mesh& mesh = field.mesh();
    const auto f = static_cast<std::size_t>(face);
    const Index owner = mesh.face_owner[f];
    const Index patch = mesh.face_patch[f];
    const auto& condition = field.boundary(patch);
    const Vec3 owner_value = field[owner];
    const double distance = boundaryNormalDistance(mesh, face);
    switch (condition.type) {
        case BoundaryType::FixedValue:
            return condition.value;
        case BoundaryType::FixedGradient:
            return owner_value + condition.value * distance;
        case BoundaryType::ZeroGradient:
            return owner_value;
        case BoundaryType::InletOutlet:
            return outward_flux >= 0.0 ? owner_value : condition.value;
        case BoundaryType::Symmetry: {
            const Vec3 normal = mesh.faceNormal(face);
            return owner_value - dot(owner_value, normal) * normal;
        }
    }
    throw std::invalid_argument("unknown vector boundary condition");
}

}  // babelsim 命名空间
