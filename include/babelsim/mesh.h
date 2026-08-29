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

    // Structure-of-arrays topology and geometry. The same arrays are used for
    // nx*ny*1 and fully three-dimensional meshes.
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

    // Local structured cell ids remain the storage indices. These compact
    // maps select the owned subset used by distributed algebra and output.
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

}  // namespace babelsim
