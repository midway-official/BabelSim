#pragma once

#include "babelsim/vector.h"

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace babelsim {

using Index = std::int32_t;
constexpr Index invalid_index = -1;

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
    std::vector<Index> faces;
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
    std::vector<Vec3> vertices;
    std::vector<Vec3> cell_centres;
    std::vector<double> cell_volumes;
    std::vector<std::array<Index, 6>> cell_faces;
    std::vector<std::array<Index, 6>> cell_neighbours;

    std::vector<std::array<Index, 4>> face_vertices;
    std::vector<Index> face_owner;
    std::vector<Index> face_neighbour;
    std::vector<Index> face_patch;
    std::vector<Vec3> face_centres;
    std::vector<Vec3> face_area_vectors;
    std::vector<Vec3> face_non_orthogonal;
    std::vector<Vec3> face_skewness;
    std::vector<double> face_areas;
    std::vector<double> face_orthogonal_coefficients;
    std::vector<double> face_owner_weights;
    std::vector<BoundaryPatch> patches;

    // 局部结构化 cell ID 仍是存储索引。这些紧凑映射选择分布式代数和输出使用的 owned 子集。
    std::vector<Index> owned_cells;
    std::vector<Index> cell_owned_indices;
    std::vector<Index> cell_global_ids;

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
    void validate() const;
};

}  // babelsim 命名空间
