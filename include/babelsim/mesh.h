#pragma once

#include "babelsim/vector.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace babelsim {
namespace detail { struct MeshAccess; }

using Index = std::int32_t;
constexpr Index invalid_index = -1;

// A zero-allocation view over one variable-length connectivity row.  Mesh
// finalization stores all rows in contiguous CSR-style arrays; this view keeps
// face/cell traversal out of the solver's hot path free of temporary vectors.
class IndexRange {
public:
    using const_iterator = const Index*;

    IndexRange() = default;
    IndexRange(const Index* first, const Index* last) : m_first(first) {
        m_size = first != nullptr && last != nullptr
            ? static_cast<std::size_t>(last - first) : 0U;
    }

    const_iterator begin() const { return m_first; }
    const_iterator end() const { return m_first == nullptr ? nullptr : m_first + m_size; }
    std::size_t size() const { return m_size; }
    bool empty() const { return m_size == 0U; }
    const Index& operator[](std::size_t i) const { return m_first[i]; }

private:
    const Index* m_first = nullptr;
    std::size_t m_size = 0U;
};

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

// A boundary face declaration used by mesh producers before owner/neighbour
// resolution.  The ring is variable length; Mesh::polyhedral resolves it to a
// face with one owner and no neighbour.
struct BoundaryFaceSpec {
    std::vector<Index> vertices;
    Index patch = invalid_index;
};

// General face input.  The vertex ring may contain any number of vertices >= 3.
// The ring orientation is normalized by Mesh::polyhedral so that the face area
// vector points out of its owner cell.  A boundary face uses invalid_index for
// neighbour and must name a physical patch.
struct PolyhedralFaceSpec {
    std::vector<Index> vertices;
    Index owner = invalid_index;
    Index neighbour = invalid_index;
    Index patch = invalid_index;
};

// Read-only partition facts used by generic performance/reporting tools.  The
// values describe storage owned by one rank; they do not expose mesh arrays or
// alter the geometry/field DSL.
struct MeshPartitionInfo {
    Index global_cells = 0;
    Index local_cells = 0;
    Index owned_cells = 0;
    Index ghost_cells = 0;
    Index local_faces = 0;
    Index owned_faces = 0;
    Index communication_faces = 0;
    Index neighbour_ranks = 0;
};

// Explicit face based polyhedral mesh.  Every face has exactly one owner and
// either one neighbour or one boundary patch.  Cell/face connectivity is
// stored in contiguous CSR arrays, so traversal remains allocation-free.
struct Mesh {
private:
    friend struct detail::MeshAccess;
    struct Storage {
        Index global_cell_count = 0;
        Index ghost_layers = 0;
        bool orthogonal_geometry = true;

        MeshStorage<Vec3> vertices;
        MeshStorage<Vec3> cell_centres;
        MeshStorage<double> cell_volumes;
        MeshStorage<double> cell_inverse_volumes;
        MeshStorage<Index> face_point_offsets;
        MeshStorage<Index> face_point_ids;
        MeshStorage<Index> cell_face_offsets;
        MeshStorage<Index> cell_face_ids;
        MeshStorage<Index> cell_face_signs;
        MeshStorage<Index> cell_vertex_offsets;
        MeshStorage<Index> cell_vertex_ids;
        MeshStorage<Index> cell_neighbour_offsets;
        MeshStorage<Index> cell_neighbour_ids;
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

    static Mesh polyhedral(
        std::vector<Vec3> vertices,
        std::vector<PolyhedralFaceSpec> faces,
        std::vector<PatchSpec> patches);

    static Mesh polyhedral(
        std::vector<Vec3> vertices,
        std::vector<std::vector<Index>> face_vertices,
        std::vector<Index> owners,
        std::vector<Index> neighbours,
        std::vector<Index> face_patches,
        std::vector<PatchSpec> patches);

    Index cellCount() const {
        return static_cast<Index>(m_storage.cell_face_offsets.size() - 1U);
    }
    Index globalCellCount() const { return m_storage.global_cell_count; }
    Index faceCount() const { return static_cast<Index>(m_storage.face_owner.size()); }
    MeshPartitionInfo partitionInfo() const;
    Index vertexCount() const { return static_cast<Index>(m_storage.vertices.size()); }
    Index owner(Index face) const { return m_storage.face_owner.at(face); }
    Index neighbour(Index face) const { return m_storage.face_neighbour.at(face); }
    Index boundaryPatch(Index face) const { return m_storage.face_patch.at(face); }
    const Vec3& faceCentre(Index face) const { return m_storage.face_centres.at(face); }
    const Vec3& cellCentre(Index cell) const { return m_storage.cell_centres.at(cell); }
    const Vec3& faceAreaVector(Index face) const { return m_storage.face_area_vectors.at(face); }
    IndexRange cellVertices(Index cell) const {
        const std::size_t c = static_cast<std::size_t>(cell);
        const Index first = m_storage.cell_vertex_offsets.at(c);
        const Index last = m_storage.cell_vertex_offsets.at(c + 1U);
        return {m_storage.cell_vertex_ids.data() + first,
                m_storage.cell_vertex_ids.data() + last};
    }
    IndexRange faceVertices(Index face) const { return facePoints(face); }
    IndexRange facePoints(Index face) const {
        const std::size_t f = static_cast<std::size_t>(face);
        const Index first = m_storage.face_point_offsets.at(f);
        const Index last = m_storage.face_point_offsets.at(f + 1U);
        return {m_storage.face_point_ids.data() + first,
                m_storage.face_point_ids.data() + last};
    }
    IndexRange cellFaces(Index cell) const {
        const std::size_t c = static_cast<std::size_t>(cell);
        const Index first = m_storage.cell_face_offsets.at(c);
        const Index last = m_storage.cell_face_offsets.at(c + 1U);
        return {m_storage.cell_face_ids.data() + first,
                m_storage.cell_face_ids.data() + last};
    }
    IndexRange cellFaceSigns(Index cell) const {
        const std::size_t c = static_cast<std::size_t>(cell);
        const Index first = m_storage.cell_face_offsets.at(c);
        const Index last = m_storage.cell_face_offsets.at(c + 1U);
        return {m_storage.cell_face_signs.data() + first,
                m_storage.cell_face_signs.data() + last};
    }
    IndexRange cellNeighbours(Index cell) const {
        const std::size_t c = static_cast<std::size_t>(cell);
        const Index first = m_storage.cell_neighbour_offsets.at(c);
        const Index last = m_storage.cell_neighbour_offsets.at(c + 1U);
        return {m_storage.cell_neighbour_ids.data() + first,
                m_storage.cell_neighbour_ids.data() + last};
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
