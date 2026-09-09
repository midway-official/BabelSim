#pragma once

#include "babelsim/vector.h"

#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace babelsim {
namespace detail { struct MeshAccess; }

using Index = std::int32_t;
constexpr Index invalid_index = -1;

// Mesh 的数组一旦构造完成，其长度决定所有整数索引和 Field 布局。该轻量容器
// 保留连续 vector 存储和索引性能，但把会改变容量的操作限制为 Mesh 的成员函数。
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

// 每个边界四边形通过全局于该 Mesh 的顶点编号及所属 patch 显式给出。
// 顶点环绕方向可以任选；Mesh 会按 owner 单元的外法向修正面方向。
struct BoundaryFaceSpec {
    std::array<Index, 4> vertices{};
    Index patch = invalid_index;
};

// 只表示显式连接的非结构六面体网格。Hex 顶点顺序采用 VTK_HEXAHEDRON：
// (0,1,2,3) 为一侧环，(4,5,6,7) 为对侧对应环。网格没有逻辑坐标、维度
// 或规则编号；单元、面所有权和 ghost 信息只服务局部并行分区。
struct Mesh {
private:
    friend struct detail::MeshAccess;
    struct Storage {
        Index global_cell_count = 0;
        Index ghost_layers = 0;
        bool orthogonal_geometry = true;

        MeshStorage<Vec3> vertices;
        MeshStorage<std::array<Index, 8>> cell_vertices;
        MeshStorage<Vec3> cell_centres;
        MeshStorage<double> cell_volumes;
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

        MeshStorage<Index> owned_cells;
        MeshStorage<Index> owned_faces;
        MeshStorage<Index> cell_owned_indices;
        MeshStorage<Index> cell_global_ids;
        MeshStorage<Index> cell_owner_ranks;
        MeshStorage<Index> cell_ghost_depths;
        MeshStorage<Index> face_global_ids;
        MeshStorage<Index> face_owner_ranks;
    } m_storage;

public:
    Mesh() = default;
    Mesh(const Mesh&) = default;
    Mesh(Mesh&&) noexcept = default;

    static Mesh unstructured(
        std::vector<Vec3> vertices,
        std::vector<std::array<Index, 8>> cells,
        std::vector<PatchSpec> patches,
        std::vector<BoundaryFaceSpec> boundary_faces);

    Index cellCount() const { return static_cast<Index>(m_storage.cell_vertices.size()); }
    Index globalCellCount() const { return m_storage.global_cell_count; }
    Index faceCount() const { return static_cast<Index>(m_storage.face_owner.size()); }
    Index vertexCount() const { return static_cast<Index>(m_storage.vertices.size()); }
    Index owner(Index face) const { return m_storage.face_owner.at(face); }
    Index neighbour(Index face) const { return m_storage.face_neighbour.at(face); }
    Index boundaryPatch(Index face) const { return m_storage.face_patch.at(face); }
    const Vec3& faceCentre(Index face) const { return m_storage.face_centres.at(face); }
    const Vec3& cellCentre(Index cell) const { return m_storage.cell_centres.at(cell); }
    const Vec3& faceAreaVector(Index face) const { return m_storage.face_area_vectors.at(face); }
    const std::array<Index, 8>& cellVertices(Index cell) const {
        return m_storage.cell_vertices.at(cell);
    }
    const std::array<Index, 4>& faceVertices(Index face) const {
        return m_storage.face_vertices.at(face);
    }
    double faceArea(Index face) const { return m_storage.face_areas.at(face); }
    double faceOrthogonalCoefficient(Index face) const {
        return m_storage.face_orthogonal_coefficients.at(face);
    }
    const Vec3& faceNonOrthogonal(Index face) const {
        return m_storage.face_non_orthogonal.at(face);
    }
    double faceOwnerWeight(Index face) const { return m_storage.face_owner_weights.at(face); }
    bool orthogonalGeometry() const { return m_storage.orthogonal_geometry; }
    bool boundaryFace(Index face) const { return m_storage.face_neighbour.at(face) == invalid_index; }
    Vec3 faceNormal(Index face) const { return m_storage.face_normals.at(face); }

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
    Index globalFaceId(Index face) const;
    Index cellOwnerRank(Index cell) const;
    Index faceOwnerRank(Index face) const;
    Mesh& operator=(const Mesh&) = default;
    Mesh& operator=(Mesh&&) noexcept = default;
    void setPartition(
        Index global_cells,
        Index layers,
        std::vector<Index> cell_ids,
        std::vector<Index> cell_owners,
        std::vector<Index> cell_depths,
        std::vector<Index> face_ids,
        std::vector<Index> face_owners,
        Index local_rank);
};

}  // babelsim 命名空间
