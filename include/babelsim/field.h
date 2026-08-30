#pragma once

#include "babelsim/mesh.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <string_view>
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
        : m_mesh(&mesh),
          m_location(location),
          m_name(std::move(name)),
          m_values(entityCount(mesh, location), std::move(initial))
    {
        if (location == FieldLocation::Cell) {
            m_boundaries.resize(mesh.patches.size());
        }
    }

    Field(const Mesh&&, FieldLocation, std::string = {}, T initial = T{}) = delete;

    const Mesh& mesh() const { return *m_mesh; }
    FieldLocation location() const { return m_location; }
    const std::string& name() const { return m_name; }
    std::size_t size() const { return m_values.size(); }

    // 计算热路径只暴露连续的只读视图；需要原地交换 halo 时使用明确命名的
    // mutableData()，避免调用者绕过存储不变量替换底层容器。
    T* mutableData() { return m_values.data(); }
    // 保留旧 API 的元素级写入兼容性；容量仍只能由构造函数决定。
    T* data() { return mutableData(); }
    const T* data() const { return m_values.data(); }
    T& operator[](Index index) { return m_values[static_cast<std::size_t>(index)]; }
    const T& operator[](Index index) const {
        return m_values[static_cast<std::size_t>(index)];
    }
    const std::vector<T>& values() const { return m_values; }
    T& at(Index index) {
        return m_values.at(checkedIndex(index));
    }
    const T& at(Index index) const {
        return m_values.at(checkedIndex(index));
    }
    // Field 的长度由 (Mesh, FieldLocation) 唯一决定，禁止外部 resize。
    // 该检查在 MPI halo、算子和输出入口调用，尽早捕获生命周期/越界错误。
    void validateStorage() const {
        if (m_mesh == nullptr || m_values.size() != entityCount(*m_mesh, m_location) ||
            (m_location == FieldLocation::Cell &&
             m_boundaries.size() != m_mesh->patches.size())) {
            throw std::logic_error("field storage invariant is violated");
        }
    }

    Field(const Field&) = default;
    Field(Field&&) noexcept = default;
    Field& operator=(const Field&) = delete;
    Field& operator=(Field&&) = delete;
    void fill(const T& value) { std::fill(m_values.begin(), m_values.end(), value); }

    void setBoundary(Index patch, BoundaryCondition<T> condition) {
        requireCellBoundary(patch);
        m_boundaries[static_cast<std::size_t>(patch)] = std::move(condition);
    }
    BoundaryCondition<T>& boundary(Index patch) {
        requireCellBoundary(patch);
        return m_boundaries[static_cast<std::size_t>(patch)];
    }
    const BoundaryCondition<T>& boundary(Index patch) const {
        requireCellBoundary(patch);
        return m_boundaries[static_cast<std::size_t>(patch)];
    }
    BoundaryCondition<T>& boundary(std::string_view patch_name) {
        return boundary(findPatch(patch_name));
    }
    const BoundaryCondition<T>& boundary(std::string_view patch_name) const {
        return boundary(findPatch(patch_name));
    }

private:
    std::size_t checkedIndex(Index index) const {
        if (index < 0 || static_cast<std::size_t>(index) >= m_values.size()) {
            throw std::out_of_range("field index is outside storage");
        }
        return static_cast<std::size_t>(index);
    }

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
        if (m_location != FieldLocation::Cell || patch < 0 ||
            static_cast<std::size_t>(patch) >= m_boundaries.size()) {
            throw std::out_of_range("field boundary patch is invalid");
        }
    }

    Index findPatch(std::string_view patch_name) const {
        if (m_location != FieldLocation::Cell) {
            throw std::logic_error("only cell fields have boundary conditions");
        }
        for (Index patch = 0; patch < static_cast<Index>(m_mesh->patches.size()); ++patch) {
            if (m_mesh->patches[static_cast<std::size_t>(patch)].name == patch_name) {
                return patch;
            }
        }
        throw std::out_of_range("field boundary patch name is unknown");
    }

    const Mesh* m_mesh;
    FieldLocation m_location;
    std::string m_name;
    std::vector<T> m_values;
    std::vector<BoundaryCondition<T>> m_boundaries;
};

using ScalarField = Field<double>;
using VectorField = Field<Vec3>;
using TensorField = Field<Tensor3>;

template <typename T>
inline BoundaryCondition<T> fixedValue(T value) {
    return BoundaryCondition<T>::fixedValue(std::move(value));
}

template <typename T>
inline BoundaryCondition<T> fixedGradient(T gradient) {
    return BoundaryCondition<T>::fixedGradient(std::move(gradient));
}

struct ZeroGradientBoundary {
    template <typename T>
    operator BoundaryCondition<T>() const {
        return BoundaryCondition<T>::zeroGradient();
    }
};

struct SymmetryBoundary {
    template <typename T>
    operator BoundaryCondition<T>() const {
        return BoundaryCondition<T>::symmetry();
    }
};

inline ZeroGradientBoundary zeroGradient() { return {}; }
inline SymmetryBoundary symmetry() { return {}; }

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
