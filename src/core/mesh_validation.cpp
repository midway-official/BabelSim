#include "babelsim/mesh.h"
#include "internal/mesh_access.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace babelsim {
namespace {

void validateStorageShape(const Mesh& mesh) {
    const auto& storage = detail::meshData(mesh);
    const Index cells = mesh.cellCount();
    const Index faces = mesh.faceCount();
    if (cells <= 0 || faces <= 0 || mesh.vertexCount() <= 0 || mesh.patchCount() <= 0 ||
        storage.global_cell_count <= 0 ||
        storage.face_point_offsets.size() != static_cast<std::size_t>(faces) + 1U ||
        storage.cell_face_offsets.size() != static_cast<std::size_t>(cells) + 1U ||
        storage.cell_face_signs.size() != storage.cell_face_ids.size() ||
        storage.cell_vertex_offsets.size() != static_cast<std::size_t>(cells) + 1U ||
        storage.cell_neighbour_offsets.size() != static_cast<std::size_t>(cells) + 1U ||
        storage.cell_centres.size() != static_cast<std::size_t>(cells) ||
        storage.cell_volumes.size() != static_cast<std::size_t>(cells) ||
        storage.cell_inverse_volumes.size() != static_cast<std::size_t>(cells) ||
        storage.cell_global_ids.size() != static_cast<std::size_t>(cells) ||
        storage.cell_owner_ranks.size() != static_cast<std::size_t>(cells) ||
        storage.cell_ghost_depths.size() != static_cast<std::size_t>(cells) ||
        storage.cell_owned_indices.size() != static_cast<std::size_t>(cells) ||
        storage.face_global_ids.size() != static_cast<std::size_t>(faces) ||
        storage.face_owner_ranks.size() != static_cast<std::size_t>(faces) ||
        storage.face_neighbour.size() != static_cast<std::size_t>(faces) ||
        storage.face_patch.size() != static_cast<std::size_t>(faces) ||
        storage.face_centres.size() != static_cast<std::size_t>(faces) ||
        storage.face_area_vectors.size() != static_cast<std::size_t>(faces) ||
        storage.face_normals.size() != static_cast<std::size_t>(faces) ||
        storage.face_areas.size() != static_cast<std::size_t>(faces) ||
        storage.face_orthogonal_coefficients.size() != static_cast<std::size_t>(faces) ||
        storage.face_owner_weights.size() != static_cast<std::size_t>(faces) ||
        storage.face_non_orthogonal.size() != static_cast<std::size_t>(faces) ||
        storage.face_skewness.size() != static_cast<std::size_t>(faces)) {
        throw std::runtime_error("polyhedral mesh storage is inconsistent");
    }
}

void validateCellConnectivity(const Mesh& mesh, std::set<Index>& cell_ids) {
    const auto& storage = detail::meshData(mesh);
    const Index cells = mesh.cellCount();
    const Index faces = mesh.faceCount();
    for (Index cell = 0; cell < cells; ++cell) {
        const std::size_t c = static_cast<std::size_t>(cell);
        if (storage.cell_global_ids[c] < 0 ||
            storage.cell_global_ids[c] >= storage.global_cell_count ||
            !cell_ids.insert(storage.cell_global_ids[c]).second ||
            storage.cell_owner_ranks[c] < 0 || storage.cell_ghost_depths[c] < 0 ||
            !(storage.cell_volumes[c] > 0.0) || !std::isfinite(storage.cell_volumes[c])) {
            throw std::runtime_error("polyhedral cell metadata is invalid");
        }
        const Index face_first = storage.cell_face_offsets[c];
        const Index face_last = storage.cell_face_offsets[c + 1U];
        const Index vertex_first = storage.cell_vertex_offsets[c];
        const Index vertex_last = storage.cell_vertex_offsets[c + 1U];
        const Index neighbour_first = storage.cell_neighbour_offsets[c];
        const Index neighbour_last = storage.cell_neighbour_offsets[c + 1U];
        if (face_first < 0 || face_last <= face_first ||
            face_last > static_cast<Index>(storage.cell_face_ids.size()) || vertex_first < 0 ||
            vertex_last <= vertex_first ||
            vertex_last > static_cast<Index>(storage.cell_vertex_ids.size()) ||
            neighbour_first < 0 || neighbour_last < neighbour_first ||
            neighbour_last > static_cast<Index>(storage.cell_neighbour_ids.size())) {
            throw std::runtime_error("polyhedral cell connectivity is invalid");
        }
        for (Index i = face_first; i < face_last; ++i) {
            const Index face = storage.cell_face_ids[static_cast<std::size_t>(i)];
            const Index sign = storage.cell_face_signs[static_cast<std::size_t>(i)];
            if (face < 0 || face >= faces || (sign != 1 && sign != -1)) {
                throw std::runtime_error("polyhedral cell-face sign is invalid");
            }
        }
        for (Index i = vertex_first; i < vertex_last; ++i) {
            const Index vertex = storage.cell_vertex_ids[static_cast<std::size_t>(i)];
            if (vertex < 0 || vertex >= mesh.vertexCount()) {
                throw std::runtime_error("polyhedral cell vertex is invalid");
            }
        }
        for (Index i = neighbour_first; i < neighbour_last; ++i) {
            const Index neighbour = storage.cell_neighbour_ids[static_cast<std::size_t>(i)];
            if (neighbour < 0 || neighbour >= cells || neighbour == cell) {
                throw std::runtime_error("polyhedral cell neighbour is invalid");
            }
        }
    }
}

void validateFaceConnectivity(const Mesh& mesh, std::set<Index>& face_ids) {
    const auto& storage = detail::meshData(mesh);
    const Index cells = mesh.cellCount();
    const Index faces = mesh.faceCount();
    for (Index face = 0; face < faces; ++face) {
        const std::size_t f = static_cast<std::size_t>(face);
        const Index first = storage.face_point_offsets[f];
        const Index last = storage.face_point_offsets[f + 1U];
        const Index owner = storage.face_owner[f];
        const Index neighbour = storage.face_neighbour[f];
        const Index patch = storage.face_patch[f];
        if (first < 0 || last - first < 3 ||
            last > static_cast<Index>(storage.face_point_ids.size()) || owner < 0 ||
            owner >= cells ||
            (neighbour != invalid_index &&
             (neighbour < 0 || neighbour >= cells || neighbour == owner)) ||
            (neighbour == invalid_index) != (patch != invalid_index) ||
            (patch != invalid_index && (patch < 0 || patch >= mesh.patchCount())) ||
            !face_ids.insert(storage.face_global_ids[f]).second || storage.face_global_ids[f] < 0 ||
            storage.face_owner_ranks[f] < 0 || !(storage.face_areas[f] > 0.0) ||
            !std::isfinite(storage.face_areas[f])) {
            throw std::runtime_error("polyhedral face metadata is invalid");
        }
        std::set<Index> unique;
        for (Index i = first; i < last; ++i) {
            const Index vertex = storage.face_point_ids[static_cast<std::size_t>(i)];
            if (vertex < 0 || vertex >= mesh.vertexCount() || !unique.insert(vertex).second) {
                throw std::runtime_error("polyhedral face vertex connectivity is invalid");
            }
        }
    }
}

void validateCellClosure(const Mesh& mesh) {
    const auto& storage = detail::meshData(mesh);
    const Index cells = mesh.cellCount();
    // A face based mesh is only usable for finite-volume assembly when every
    // cell is a closed oriented surface.  Check both the vector area closure
    // and the edge manifold here, once at mesh setup time.  The edge test is
    // deliberately performed on each cell's outward orientation, so explicit
    // coarse/fine subfaces are valid as long as their shared edge appears once
    // in each adjacent face with opposite direction.
    for (Index cell = 0; cell < cells; ++cell) {
        const std::size_t c = static_cast<std::size_t>(cell);
        Vec3 closure{};
        double total_area = 0.0;
        std::map<std::pair<Index, Index>, std::pair<Index, Index>> edges;
        const Index first = storage.cell_face_offsets[c];
        const Index last = storage.cell_face_offsets[c + 1U];
        for (Index position = first; position < last; ++position) {
            const Index face = storage.cell_face_ids[static_cast<std::size_t>(position)];
            const Index sign = storage.cell_face_signs[static_cast<std::size_t>(position)];
            const std::size_t f = static_cast<std::size_t>(face);
            closure += static_cast<double>(sign) * storage.face_area_vectors[f];
            total_area += storage.face_areas[f];
            const Index point_first = storage.face_point_offsets[f];
            const Index point_last = storage.face_point_offsets[f + 1U];
            std::vector<Index> ring;
            ring.reserve(static_cast<std::size_t>(point_last - point_first));
            for (Index point = point_first; point < point_last; ++point) {
                ring.push_back(storage.face_point_ids[static_cast<std::size_t>(point)]);
            }
            if (sign < 0)
                std::reverse(ring.begin(), ring.end());
            for (std::size_t edge = 0; edge < ring.size(); ++edge) {
                const Index from = ring[edge];
                const Index to = ring[(edge + 1U) % ring.size()];
                const auto key = std::minmax(from, to);
                auto& count_and_orientation = edges[key];
                ++count_and_orientation.first;
                count_and_orientation.second += from == key.first ? 1 : -1;
            }
        }
        const double closure_tolerance = 1e-9 * std::max(total_area, 1.0);
        if (!(std::isfinite(total_area) && norm(closure) <= closure_tolerance)) {
            throw std::runtime_error("polyhedral cell face area vectors do not close");
        }
        for (const auto& [edge, count_and_orientation] : edges) {
            if (count_and_orientation.first != 2 || count_and_orientation.second != 0) {
                throw std::runtime_error(
                    "polyhedral cell surface has a non-manifold or inconsistently oriented edge " +
                    std::to_string(edge.first) + "-" + std::to_string(edge.second));
            }
        }
    }
}

void validatePatchAssignments(const Mesh& mesh, std::vector<bool>& patch_seen) {
    const auto& storage = detail::meshData(mesh);
    const Index faces = mesh.faceCount();
    for (Index patch = 0; patch < mesh.patchCount(); ++patch) {
        const BoundaryPatch& boundary = storage.patches[static_cast<std::size_t>(patch)];
        if (boundary.name.empty())
            throw std::runtime_error("mesh has an unnamed patch");
        for (Index face : boundary.faces) {
            if (face < 0 || face >= faces ||
                storage.face_patch[static_cast<std::size_t>(face)] != patch ||
                patch_seen[static_cast<std::size_t>(face)]) {
                throw std::runtime_error("polyhedral patch connectivity is invalid");
            }
            patch_seen[static_cast<std::size_t>(face)] = true;
        }
    }
    for (Index face = 0; face < faces; ++face) {
        if (storage.face_neighbour[static_cast<std::size_t>(face)] == invalid_index &&
            !patch_seen[static_cast<std::size_t>(face)]) {
            throw std::runtime_error("polyhedral boundary face is not assigned to a patch");
        }
    }
}

void validateOwnership(const Mesh& mesh) {
    const auto& storage = detail::meshData(mesh);
    const Index cells = mesh.cellCount();
    Index owned_count = 0;
    for (Index cell = 0; cell < cells; ++cell) {
        const Index owned = storage.cell_owned_indices[static_cast<std::size_t>(cell)];
        if (owned != invalid_index) {
            if (owned != owned_count++ ||
                storage.cell_ghost_depths[static_cast<std::size_t>(cell)] != 0) {
                throw std::runtime_error("owned cell mapping is invalid");
            }
        }
    }
    if (owned_count != detail::ownedCellCount(mesh) || owned_count == 0) {
        throw std::runtime_error("owned cell count is invalid");
    }
}

} // namespace

void Mesh::validate() const {
    // Validate lengths before dereferencing connectivity, then metadata before
    // closure and patch membership. Failure order is part of the diagnostics.
    validateStorageShape(*this);
    // Keep the scratch containers alive until all checks finish, as in the
    // original validator. Splitting checks must not reshuffle allocator reuse
    // before fields and PETSc work buffers are subsequently constructed.
    std::set<Index> cell_ids;
    validateCellConnectivity(*this, cell_ids);
    std::set<Index> face_ids;
    std::vector<bool> patch_seen(static_cast<std::size_t>(faceCount()), false);
    validateFaceConnectivity(*this, face_ids);
    validateCellClosure(*this);
    validatePatchAssignments(*this, patch_seen);
    validateOwnership(*this);
}

} // namespace babelsim
