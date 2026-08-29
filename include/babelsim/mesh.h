#pragma once

#include "babelsim/vector.h"

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace babelsim {

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

    std::size_t size() const { return data_.size(); }
    bool empty() const { return data_.empty(); }
    T& operator[](std::size_t index) { return data_[index]; }
    const T& operator[](std::size_t index) const { return data_[index]; }
    T& at(std::size_t index) { return data_.at(index); }
    const T& at(std::size_t index) const { return data_.at(index); }
    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }
    T& front() { return data_.front(); }
    const T& front() const { return data_.front(); }
    auto begin() { return data_.begin(); }
    auto end() { return data_.end(); }
    auto begin() const { return data_.begin(); }
    auto end() const { return data_.end(); }

private:
    friend struct Mesh;
    void reserve(std::size_t count) { data_.reserve(count); }
    void resize(std::size_t count) { data_.resize(count); }
    void assign(std::size_t count, const T& value) { data_.assign(count, value); }
    void assign(std::vector<T>&& values) { data_ = std::move(values); }
    void clear() { data_.clear(); }
    void push_back(const T& value) { data_.push_back(value); }
    void push_back(T&& value) { data_.push_back(std::move(value)); }
    MeshStorage& operator=(const MeshStorage&) = default;
    MeshStorage& operator=(MeshStorage&&) noexcept = default;

    std::vector<T> data_;
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
    std::array<Index, 3> dimensions{};
    std::array<Index, 3> global_dimensions{};
    Index global_i_offset = 0;
    Index owned_i_begin = 0;
    Index owned_i_end = 0;
    Index ghost_layers = 0;

    // 结构化数组形式的拓扑与几何。相同数组同时服务 nx*ny*1 和完整三维网格。
    MeshStorage<Vec3> vertices;
    MeshStorage<Vec3> cell_centres;
    MeshStorage<double> cell_volumes;
    MeshStorage<std::array<Index, 6>> cell_faces;
    MeshStorage<std::array<Index, 6>> cell_neighbours;

    MeshStorage<std::array<Index, 4>> face_vertices;
    MeshStorage<Index> face_owner;
    MeshStorage<Index> face_neighbour;
    MeshStorage<Index> face_patch;
    MeshStorage<Vec3> face_centres;
    MeshStorage<Vec3> face_area_vectors;
    MeshStorage<Vec3> face_non_orthogonal;
    MeshStorage<Vec3> face_skewness;
    MeshStorage<double> face_areas;
    MeshStorage<double> face_orthogonal_coefficients;
    MeshStorage<double> face_owner_weights;
    MeshStorage<BoundaryPatch> patches;

    // 局部结构化 cell ID 仍是存储索引。这些紧凑映射选择分布式代数和输出使用的 owned 子集。
    MeshStorage<Index> owned_cells;
    MeshStorage<Index> cell_owned_indices;
    MeshStorage<Index> cell_global_ids;

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
    Index ownedCellCount() const {
        return static_cast<Index>(owned_cells.size());
    }
    Index faceCount() const { return static_cast<Index>(face_owner.size()); }
    Index vertexCount() const { return static_cast<Index>(vertices.size()); }
    Index owner(Index face) const { return face_owner[static_cast<std::size_t>(face)]; }
    Index neighbour(Index face) const {
        return face_neighbour[static_cast<std::size_t>(face)];
    }
    const Vec3& faceCentre(Index face) const {
        return face_centres[static_cast<std::size_t>(face)];
    }
    const Vec3& cellCentre(Index cell) const {
        return cell_centres[static_cast<std::size_t>(cell)];
    }
    const Vec3& faceAreaVector(Index face) const {
        return face_area_vectors[static_cast<std::size_t>(face)];
    }
    double faceArea(Index face) const { return face_areas[static_cast<std::size_t>(face)]; }
    double faceOrthogonalCoefficient(Index face) const {
        return face_orthogonal_coefficients[static_cast<std::size_t>(face)];
    }
    const Vec3& faceNonOrthogonal(Index face) const {
        return face_non_orthogonal[static_cast<std::size_t>(face)];
    }
    double faceOwnerWeight(Index face) const {
        return face_owner_weights[static_cast<std::size_t>(face)];
    }
    bool boundaryFace(Index face) const {
        return face_neighbour[static_cast<std::size_t>(face)] == invalid_index;
    }
    Vec3 faceNormal(Index face) const {
        const auto n = static_cast<std::size_t>(face);
        return face_area_vectors[n] / face_areas[n];
    }

    Index cellId(Index i, Index j, Index k) const;
    Index vertexId(Index i, Index j, Index k) const;
    bool isOwned(Index cell) const;
    Index ownedIndex(Index cell) const;
    Index globalCellId(Index cell) const;
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
    void validate() const;
};

}  // babelsim 命名空间
