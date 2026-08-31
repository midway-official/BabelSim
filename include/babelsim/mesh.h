#pragma once

#include "babelsim/vector.h"

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace babelsim {
namespace detail { struct MeshAccess; }

using Index = std::int32_t;
constexpr Index invalid_index = -1;

// Mesh 的数组一旦构造完成，其长度决定所有整数索引和 Field 布局。该轻量容器
// 保留连续 vector 存储和索引性能，但把会改变容量的操作限制为 Mesh 的成员函数。
// 外部仍可在受控的构造阶段修改元素值，不能通过公开 API resize/clear/push_back。
template <typename T>
class MeshStorage {
public:
    MeshStorage() = default;
    MeshStorage(const MeshStorage&) = default;
    MeshStorage(MeshStorage&&) noexcept = default;

    std::size_t size() const { return m_data.size(); }
    bool empty() const { return m_data.empty(); }
    T& operator[](std::size_t index) { return m_data[index]; }
    const T& operator[](std::size_t index) const { return m_data[index]; }
    T& at(std::size_t index) { return m_data.at(index); }
    const T& at(std::size_t index) const { return m_data.at(index); }
    T* data() { return m_data.data(); }
    const T* data() const { return m_data.data(); }
    T& front() { return m_data.front(); }
    const T& front() const { return m_data.front(); }
    auto begin() { return m_data.begin(); }
    auto end() { return m_data.end(); }
    auto begin() const { return m_data.begin(); }
    auto end() const { return m_data.end(); }

private:
    friend struct Mesh;
    void reserve(std::size_t count) { m_data.reserve(count); }
    void resize(std::size_t count) { m_data.resize(count); }
    void assign(std::size_t count, const T& value) { m_data.assign(count, value); }
    void assign(std::vector<T>&& values) { m_data = std::move(values); }
    void clear() { m_data.clear(); }
    void push_back(const T& value) { m_data.push_back(value); }
    void push_back(T&& value) { m_data.push_back(std::move(value)); }
    MeshStorage& operator=(const MeshStorage&) = default;
    MeshStorage& operator=(MeshStorage&&) noexcept = default;

    std::vector<T> m_data;
};

enum class Side : std::size_t {
    XMin,
    XMax,
    YMin,
    YMax,
    ZMin,
    ZMax,
};

enum class PatchKind {
    Generic,
    Wall,
    Inlet,
    Outlet,
    Symmetry,
    Processor,
};

struct PatchSpec {
    std::string name;
    PatchKind kind = PatchKind::Generic;
};

struct BoundaryPatch {
    std::string name;
    PatchKind kind = PatchKind::Generic;
    MeshStorage<Index> faces;
};

std::array<PatchSpec, 6> defaultPatches();

struct Mesh {
private:
    friend struct detail::MeshAccess;
    struct Storage {
        std::array<Index, 3> dimensions{};
        std::array<Index, 3> global_dimensions{};
        Index global_i_offset = 0;
        Index owned_i_begin = 0;
        Index owned_i_end = 0;
        Index ghost_layers = 0;
        // 网格构造阶段识别的正交几何标志。正交网格无需重复执行偏斜面重构。
        bool orthogonal_geometry = true;

        // 结构化数组形式的拓扑与几何。相同数组同时服务 nx*ny*1 和完整三维网格。
        MeshStorage<Vec3> vertices;
        MeshStorage<Vec3> cell_centres;
        MeshStorage<double> cell_volumes;
        // 常用几何量的倒数/归一化缓存，避免算子热路径重复做除法。
        MeshStorage<double> cell_inverse_volumes;
        MeshStorage<std::array<Index, 6>> cell_faces;
        MeshStorage<std::array<Index, 6>> cell_neighbours;

        MeshStorage<std::array<Index, 4>> face_vertices;
        MeshStorage<Index> face_owner;
        MeshStorage<Index> face_neighbour;
        MeshStorage<Index> face_patch;
        MeshStorage<Vec3> face_centres;
        MeshStorage<Vec3> face_area_vectors;
        MeshStorage<Vec3> face_normals;
        MeshStorage<Vec3> face_non_orthogonal;
        MeshStorage<Vec3> face_skewness;
        MeshStorage<double> face_areas;
        MeshStorage<double> face_orthogonal_coefficients;
        MeshStorage<double> face_owner_weights;
        MeshStorage<BoundaryPatch> patches;

        // 局部结构化 cell ID 仍是存储索引。这些紧凑映射选择分布式代数和输出使用的 owned 子集。
        MeshStorage<Index> owned_cells;
        // 至少连接一个 owned cell 的面。分区后代数装配和物理更新无需扫描 ghost-only 面。
        MeshStorage<Index> owned_faces;
        MeshStorage<Index> cell_owned_indices;
        MeshStorage<Index> cell_global_ids;
    } m_storage;

public:
    Mesh() = default;
    Mesh(const Mesh&) = default;
    Mesh(Mesh&&) noexcept = default;

    static Mesh structured(
        std::array<Index, 3> cells,
        std::vector<Vec3> points,
        const std::array<PatchSpec, 6>& boundary = defaultPatches());

    static Mesh cartesian(
        std::array<Index, 3> cells,
        Vec3 minimum,
        Vec3 maximum,
        const std::array<PatchSpec, 6>& boundary = defaultPatches());

    Index cellCount() const;
    Index faceCount() const { return static_cast<Index>(m_storage.face_owner.size()); }
    Index vertexCount() const { return static_cast<Index>(m_storage.vertices.size()); }
    Index owner(Index face) const { return m_storage.face_owner[static_cast<std::size_t>(face)]; }
    Index neighbour(Index face) const {
        return m_storage.face_neighbour[static_cast<std::size_t>(face)];
    }
    const Vec3& faceCentre(Index face) const {
        return m_storage.face_centres[static_cast<std::size_t>(face)];
    }
    const Vec3& cellCentre(Index cell) const {
        return m_storage.cell_centres[static_cast<std::size_t>(cell)];
    }
    const Vec3& faceAreaVector(Index face) const {
        return m_storage.face_area_vectors[static_cast<std::size_t>(face)];
    }
    double faceArea(Index face) const { return m_storage.face_areas[static_cast<std::size_t>(face)]; }
    double faceOrthogonalCoefficient(Index face) const {
        return m_storage.face_orthogonal_coefficients[static_cast<std::size_t>(face)];
    }
    const Vec3& faceNonOrthogonal(Index face) const {
        return m_storage.face_non_orthogonal[static_cast<std::size_t>(face)];
    }
    double faceOwnerWeight(Index face) const {
        return m_storage.face_owner_weights[static_cast<std::size_t>(face)];
    }
    bool orthogonalGeometry() const { return m_storage.orthogonal_geometry; }
    bool boundaryFace(Index face) const {
        return m_storage.face_neighbour[static_cast<std::size_t>(face)] == invalid_index;
    }
    Vec3 faceNormal(Index face) const {
        return m_storage.face_normals[static_cast<std::size_t>(face)];
    }

    Index cellId(Index i, Index j, Index k) const;
    Index vertexId(Index i, Index j, Index k) const;
    Index patchCount() const { return static_cast<Index>(m_storage.patches.size()); }
    const std::string& patchName(Index patch) const { return m_storage.patches.at(patch).name; }
    PatchKind patchKind(Index patch) const { return m_storage.patches.at(patch).kind; }
    double cellVolume(Index cell) const { return m_storage.cell_volumes.at(cell); }
    const Vec3& vertex(Index index) const { return m_storage.vertices.at(index); }
    void validate() const;

private:
    Index ownedCellCount() const { return static_cast<Index>(m_storage.owned_cells.size()); }
    bool isOwned(Index cell) const;
    Index ownedIndex(Index cell) const;
    Index globalCellId(Index cell) const;
    Mesh& operator=(const Mesh&) = default;
    Mesh& operator=(Mesh&&) noexcept = default;
    void setOwnership(
        std::array<Index, 3> global_cells,
        Index global_x_offset,
        Index owned_x_begin,
        Index owned_x_end,
        Index halo_layers);
    // 仅用于构造/分区阶段一次性替换 patch 拓扑；调用后应立即 validate()。
    void setPatches(std::vector<BoundaryPatch> patches);
    // 受控地登记一个边界面，避免调用者直接改变 patch 面列表容量。
    void addPatchFace(Index patch, Index face);
};

}  // babelsim 命名空间
