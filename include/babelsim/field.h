#pragma once

#include "babelsim/mesh.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace babelsim {
namespace detail { struct FieldAccess; }

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
            m_boundaries.resize(mesh.patchCount());
        }
    }

    Field(const Mesh&&, FieldLocation, std::string = {}, T initial = T{}) = delete;

    const Mesh& mesh() const { return *m_mesh; }
    FieldLocation location() const { return m_location; }
    const std::string& name() const { return m_name; }
    std::size_t size() const { return m_values.size(); }

private:
    friend struct detail::FieldAccess;
    template <typename> friend class Field;
public:
    // Field 的长度由 (Mesh, FieldLocation) 唯一决定，禁止外部 resize。
    // 该检查在 MPI halo、算子和输出入口调用，尽早捕获生命周期/越界错误。
    void validateStorage() const {
        if (m_mesh == nullptr || m_values.size() != entityCount(*m_mesh, m_location) ||
            (m_location == FieldLocation::Cell &&
             m_boundaries.size() != static_cast<std::size_t>(m_mesh->patchCount()))) {
            throw std::logic_error("field storage invariant is violated");
        }
    }

    Field(const Field&) = default;
    Field(Field&&) noexcept = default;
    Field& operator=(const Field&) = delete;
    Field& operator=(Field&&) = delete;
    void fill(const T& value) { std::fill(m_values.begin(), m_values.end(), value); }

    // 按空间位置定义已知场（初值、物性或源）。函数应只依赖位置和捕获的物理参数，
    // 不依赖调用次数或分区；框架遍历正确的数据位置，Solver 不接触本地索引。
    template <typename Function>
    void evaluate(Function function) {
        validateStorage();
        for (Index index = 0; index < static_cast<Index>(m_values.size()); ++index) {
            const Vec3& position = m_location == FieldLocation::Cell
                ? m_mesh->cellCentre(index)
                : m_location == FieldLocation::Face ? m_mesh->faceCentre(index)
                                                   : m_mesh->vertex(index);
            m_values[static_cast<std::size_t>(index)] = function(position);
        }
    }

    // 点值物性/源关系，例如 k(T) 或动能(U)。输入输出可有不同值类型，布局必须相同；
    // 不创建临时 Field，且不把数据指针传给用户函数。
    template <typename U, typename Function>
    void evaluate(const Field<U>& source, Function function) {
        validateStorage();
        source.validateStorage();
        if (m_mesh != source.m_mesh || m_location != source.m_location)
            throw std::invalid_argument("field evaluation requires the same mesh and location");
        for (std::size_t index = 0; index < m_values.size(); ++index)
            m_values[index] = function(source.m_values[index]);
    }

    // 显式场赋值保留 Mesh、位置、名称和边界定义，只复制数值。它用于算法历史场和
    // 已知物性变换，避免 Solver 接触底层连续存储或重新分配容器。
    void assign(const Field& source) {
        requireCompatible(source, "field assignment");
        std::copy(source.m_values.begin(), source.m_values.end(), m_values.begin());
    }

    void assignScaled(double factor, const Field& source) {
        if (!std::isfinite(factor)) {
            throw std::invalid_argument("field scale factor must be finite");
        }
        requireCompatible(source, "field scaling");
        std::transform(
            source.m_values.begin(), source.m_values.end(), m_values.begin(),
            [factor](const T& value) { return factor * value; });
    }

    void addScaled(double factor, const Field& source) {
        if (!std::isfinite(factor)) {
            throw std::invalid_argument("field scale factor must be finite");
        }
        requireCompatible(source, "field scaled addition");
        for (std::size_t index = 0; index < m_values.size(); ++index) {
            m_values[index] += factor * source.m_values[index];
        }
    }

    // 通用逐点乘积。该操作覆盖 owned+ghost 的连续本地存储，使 Physics 不需要
    // 手写 cell 循环；分布式输入同步仍由调用它的 math/Algorithm 步骤负责。
    void assignProduct(const Field<double>& coefficient, const Field& source) {
        requireCompatible(source, "field product");
        coefficient.validateStorage();
        if (m_mesh != &coefficient.mesh() || m_location != coefficient.location()) {
            throw std::invalid_argument(
                "field product requires fields on the same mesh and location");
        }
        for (std::size_t index = 0; index < m_values.size(); ++index) {
            m_values[index] = coefficient.m_values[index] *
                source.m_values[index];
        }
    }

    void addProduct(
        double factor,
        const Field<double>& coefficient,
        const Field& source)
    {
        if (!std::isfinite(factor)) {
            throw std::invalid_argument("field product factor must be finite");
        }
        requireCompatible(source, "field product addition");
        coefficient.validateStorage();
        if (m_mesh != &coefficient.mesh() || m_location != coefficient.location()) {
            throw std::invalid_argument(
                "field product requires fields on the same mesh and location");
        }
        for (std::size_t index = 0; index < m_values.size(); ++index) {
            m_values[index] += factor * coefficient.m_values[index] *
                source.m_values[index];
        }
    }

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
    void requireCompatible(const Field& source, const char* operation) const {
        validateStorage();
        source.validateStorage();
        if (m_mesh != source.m_mesh || m_location != source.m_location) {
            throw std::invalid_argument(std::string(operation) +
                                        " requires fields on the same mesh and location");
        }
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
        for (Index patch = 0; patch < static_cast<Index>(m_mesh->patchCount()); ++patch) {
            if (m_mesh->patchName(patch) == patch_name) {
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

// 根据原 Field 的边界生成增量/修正 Field 的齐次边界：固定值变为零固定值，
// symmetry 保持 symmetry，其余边界对应零法向梯度。返回值用于判断方程是否需要参考点。
template <typename T>
bool setHomogeneousCorrectionBoundaries(
    Field<T>& correction,
    const Field<T>& reference)
{
    correction.validateStorage();
    reference.validateStorage();
    if (&correction.mesh() != &reference.mesh() ||
        correction.location() != FieldLocation::Cell ||
        reference.location() != FieldLocation::Cell) {
        throw std::invalid_argument(
            "correction boundaries require cell fields on the same mesh");
    }
    bool has_fixed_value = false;
    const Mesh& mesh = reference.mesh();
    for (Index patch = 0; patch < static_cast<Index>(mesh.patchCount()); ++patch) {
        const BoundaryType type = reference.boundary(patch).type;
        if (type == BoundaryType::FixedValue) {
            correction.setBoundary(patch, BoundaryCondition<T>::fixedValue(T{}));
            has_fixed_value = true;
        } else if (type == BoundaryType::Symmetry) {
            correction.setBoundary(patch, BoundaryCondition<T>::symmetry());
        } else {
            correction.setBoundary(patch, BoundaryCondition<T>::zeroGradient());
        }
    }
    return has_fixed_value;
}

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

}  // babelsim 命名空间
