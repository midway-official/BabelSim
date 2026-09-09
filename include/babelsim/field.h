#pragma once

#include "babelsim/mesh.h"

#include <algorithm>
#include <cmath>
#include <limits>
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
          m_values(entityCount(mesh, location), std::move(initial)),
          m_halo_valid(true)
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
            (m_calculated_boundary && m_boundary_trace.size() != static_cast<std::size_t>(m_mesh->faceCount())) ||
            (!m_boundary_flux.empty() && m_boundary_flux.size() != static_cast<std::size_t>(m_mesh->faceCount())) ||
            (m_location == FieldLocation::Cell &&
             m_boundaries.size() != static_cast<std::size_t>(m_mesh->patchCount()))) {
            throw std::logic_error("field storage invariant is violated");
        }
    }

    Field(const Field&) = default;
    Field(Field&&) noexcept = default;
    Field& operator=(const Field&) = delete;
    Field& operator=(Field&&) = delete;
    void fill(const T& value) {
        updateBoundaryTrace([&](Index) { return value; });
        std::fill(m_values.begin(), m_values.end(), value);
        // owned 与 ghost 同时被同一常量覆盖，不需要再进行 halo 交换。
        m_halo_valid = true;
    }

    // 按空间位置定义已知场（初值、物性或源）。函数应只依赖位置和捕获的物理参数，
    // 不依赖调用次数或分区；框架遍历正确的数据位置，Solver 不接触本地索引。
    template <typename Function>
    void evaluate(Function function) {
        validateStorage();
        updateBoundaryTrace([&](Index face) { return function(m_mesh->faceCentre(face)); });
        for (Index index = 0; index < static_cast<Index>(m_values.size()); ++index) {
            const Vec3& position = m_location == FieldLocation::Cell
                ? m_mesh->cellCentre(index)
                : m_location == FieldLocation::Face ? m_mesh->faceCentre(index)
                                                   : m_mesh->vertex(index);
            m_values[static_cast<std::size_t>(index)] = function(position);
        }
        // 坐标函数在本地完整布局上求值，分区两侧的重复实体天然一致。
        m_halo_valid = true;
    }

    // 点值物性/源关系，例如 k(T) 或动能(U)。输入输出可有不同值类型，布局必须相同；
    // 不创建临时 Field，且不把数据指针传给用户函数。
    template <typename U, typename Function>
    void evaluate(const Field<U>& source, Function function) {
        validateStorage();
        source.validateStorage();
        if (m_mesh != source.m_mesh || m_location != source.m_location)
            throw std::invalid_argument("field evaluation requires the same mesh and location");
        updateBoundaryTrace([&](Index face) { return function(source.boundaryTrace(face)); });
        for (std::size_t index = 0; index < m_values.size(); ++index)
            m_values[index] = function(source.m_values[index]);
        m_halo_valid = source.m_halo_valid;
    }

    // 显式场赋值保留 Mesh、位置、名称和边界定义，只复制数值。它用于算法历史场和
    // 已知物性变换，避免 Solver 接触底层连续存储或重新分配容器。
    void assign(const Field& source) {
        requireCompatible(source, "field assignment");
        updateBoundaryTrace([&](Index face) { return source.boundaryTrace(face); });
        std::copy(source.m_values.begin(), source.m_values.end(), m_values.begin());
        m_halo_valid = source.m_halo_valid;
    }

    void assignScaled(double factor, const Field& source) {
        if (!std::isfinite(factor)) {
            throw std::invalid_argument("field scale factor must be finite");
        }
        requireCompatible(source, "field scaling");
        updateBoundaryTrace([&](Index face) { return factor * source.boundaryTrace(face); });
        std::transform(
            source.m_values.begin(), source.m_values.end(), m_values.begin(),
            [factor](const T& value) { return factor * value; });
        m_halo_valid = source.m_halo_valid;
    }

    void addScaled(double factor, const Field& source) {
        if (!std::isfinite(factor)) {
            throw std::invalid_argument("field scale factor must be finite");
        }
        requireCompatible(source, "field scaled addition");
        updateBoundaryTrace([&](Index face) {
            return boundaryTrace(face) + factor * source.boundaryTrace(face);
        });
        for (std::size_t index = 0; index < m_values.size(); ++index) {
            m_values[index] += factor * source.m_values[index];
        }
        m_halo_valid = m_halo_valid && source.m_halo_valid;
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
        updateBoundaryTrace([&](Index face) {
            return coefficient.boundaryTrace(face) * source.boundaryTrace(face);
        });
        for (std::size_t index = 0; index < m_values.size(); ++index) {
            m_values[index] = coefficient.m_values[index] *
                source.m_values[index];
        }
        m_halo_valid = coefficient.m_halo_valid && source.m_halo_valid;
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
        updateBoundaryTrace([&](Index face) {
            return boundaryTrace(face) + factor * coefficient.boundaryTrace(face) * source.boundaryTrace(face);
        });
        for (std::size_t index = 0; index < m_values.size(); ++index) {
            m_values[index] += factor * coefficient.m_values[index] *
                source.m_values[index];
        }
        m_halo_valid = m_halo_valid && coefficient.m_halo_valid && source.m_halo_valid;
    }

    void setBoundary(Index patch, BoundaryCondition<T> condition) {
        requireCellBoundary(patch);
        if (m_calculated_boundary) throw std::logic_error("calculated field has boundary traces, not equation constraints");
        m_boundaries[static_cast<std::size_t>(patch)] = std::move(condition);
    }
    BoundaryCondition<T>& boundary(Index patch) {
        requireCellBoundary(patch);
        if (m_calculated_boundary) throw std::logic_error("calculated field has per-face traces, not patch constraints");
        return m_boundaries[static_cast<std::size_t>(patch)];
    }
    const BoundaryCondition<T>& boundary(Index patch) const {
        requireCellBoundary(patch);
        if (m_calculated_boundary) throw std::logic_error("calculated field has per-face traces, not patch constraints");
        return m_boundaries[static_cast<std::size_t>(patch)];
    }
    BoundaryCondition<T>& boundary(std::string_view patch_name) {
        return boundary(findPatch(patch_name));
    }
    const BoundaryCondition<T>& boundary(std::string_view patch_name) const {
        return boundary(findPatch(patch_name));
    }

    // 未知场的值更新保留其边界约束。已知/派生场可显式选择 calculated 模式：
    // 所有逐点代数同时作用于单元值和边界迹，不推断字段名或物理意义。
    // 迹是当前操作的快照，不保留源字段或用户函数的悬空引用。
    void useCalculatedBoundary() {
        if (m_location != FieldLocation::Cell || m_calculated_boundary) return;
        std::vector<T> values(m_mesh->faceCount());
        for (Index f = 0; f < m_mesh->faceCount(); ++f)
            if (m_mesh->boundaryFace(f)) values[f] = boundaryTrace(f);
        m_boundary_trace = std::move(values);
        m_calculated_boundary = true;
    }
    bool calculatedBoundary() const { return m_calculated_boundary; }

    // inletOutlet 是依赖面通量符号的混合数学约束。保存求值上下文的快照；
    // eqn 的对流项自动提供它，独立数学操作也可由调用者显式设置。
    void setBoundaryFlux(const Field<double>& flux) {
        if (m_location != FieldLocation::Cell || &flux.mesh() != m_mesh ||
            flux.location() != FieldLocation::Face)
            throw std::invalid_argument("boundary flux must be a face field on the same mesh");
        m_boundary_flux = flux.m_values;
    }

private:
    BoundaryCondition<T> faceCondition(Index face,
        double outward_flux = std::numeric_limits<double>::quiet_NaN()) const {
        if (!m_mesh->boundaryFace(face)) throw std::invalid_argument("expected boundary face");
        if (m_calculated_boundary) return BoundaryCondition<T>::fixedValue(m_boundary_trace.at(face));
        auto condition = boundary(m_mesh->boundaryPatch(face));
        if (condition.type == BoundaryType::InletOutlet) {
            if (!std::isfinite(outward_flux)) {
                if (m_boundary_flux.empty())
                    throw std::invalid_argument("inletOutlet requires a boundary flux context");
                outward_flux = m_boundary_flux.at(face);
            }
            if (!std::isfinite(outward_flux)) throw std::invalid_argument("nonfinite boundary flux");
            return outward_flux < 0.0 ? BoundaryCondition<T>::fixedValue(condition.value)
                                     : BoundaryCondition<T>::zeroGradient();
        }
        return condition;
    }
    T boundaryTrace(Index face,
        double outward_flux = std::numeric_limits<double>::quiet_NaN()) const {
        const auto condition = faceCondition(face, outward_flux);
        const T& owner = m_values.at(m_mesh->owner(face));
        if (condition.type == BoundaryType::FixedValue) return condition.value;
        if (condition.type == BoundaryType::FixedGradient)
            return owner + dot(m_mesh->faceCentre(face) - m_mesh->cellCentre(m_mesh->owner(face)),
                               m_mesh->faceNormal(face)) * condition.value;
        if (condition.type == BoundaryType::Symmetry)
            return symmetricBoundaryValue(owner, m_mesh->faceNormal(face));
        return owner;
    }
    template <typename Function> void updateBoundaryTrace(Function function) {
        if (!m_calculated_boundary) return;
        // 每个面的计算先于单元值更新，因此支持 f.assignProduct(a,f) 等原位代数。
        for (Index face = 0; face < m_mesh->faceCount(); ++face)
            if (m_mesh->boundaryFace(face)) m_boundary_trace[face] = function(face);
    }
public:

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
    bool m_calculated_boundary = false;
    std::vector<T> m_boundary_trace;
    std::vector<double> m_boundary_flux;
    // 仅由 Field 与计算后端维护。Solver 看见的仍是完整数学场，不接触 ghost 状态。
    bool m_halo_valid = false;
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
