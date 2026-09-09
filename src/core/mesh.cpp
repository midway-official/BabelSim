#include "babelsim/mesh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace babelsim {
namespace {

using Quad = std::array<Index, 4>;
using Hex = std::array<Index, 8>;

constexpr std::array<std::array<int, 4>, 6> hex_faces{{
    {{0, 4, 7, 3}}, {{1, 2, 6, 5}}, {{0, 1, 5, 4}},
    {{3, 7, 6, 2}}, {{0, 3, 2, 1}}, {{4, 5, 6, 7}},
}};

Quad canonical(Quad vertices) {
    std::sort(vertices.begin(), vertices.end());
    return vertices;
}

void requireIndex(Index value, Index count, const char* what) {
    if (value < 0 || value >= count) {
        throw std::invalid_argument(std::string(what) + " references an invalid vertex");
    }
}

std::pair<Vec3, Vec3> quadGeometry(const Quad& ids, const MeshStorage<Vec3>& vertices) {
    const Vec3& a = vertices[static_cast<std::size_t>(ids[0])];
    const Vec3& b = vertices[static_cast<std::size_t>(ids[1])];
    const Vec3& c = vertices[static_cast<std::size_t>(ids[2])];
    const Vec3& d = vertices[static_cast<std::size_t>(ids[3])];
    const Vec3 area_abc = 0.5 * cross(b - a, c - a);
    const Vec3 area_acd = 0.5 * cross(c - a, d - a);
    const double magnitude_abc = norm(area_abc);
    const double magnitude_acd = norm(area_acd);
    const double magnitude = magnitude_abc + magnitude_acd;
    if (!(magnitude > 0.0) || !std::isfinite(magnitude)) {
        throw std::runtime_error("unstructured mesh contains a degenerate quadrilateral");
    }
    const Vec3 centre =
        (magnitude_abc * ((a + b + c) / 3.0) +
         magnitude_acd * ((a + c + d) / 3.0)) / magnitude;
    return {centre, area_abc + area_acd};
}

void addTetrahedron(
    const Vec3& reference,
    const Vec3& a,
    const Vec3& b,
    const Vec3& c,
    double& volume,
    Vec3& first_moment)
{
    const double tetra_volume =
        std::abs(dot(a - reference, cross(b - reference, c - reference))) / 6.0;
    volume += tetra_volume;
    first_moment += tetra_volume * ((reference + a + b + c) / 4.0);
}

double signedHexVolume(const Hex& cell, const MeshStorage<Vec3>& vertices) {
    double volume = 0.0;
    for (const auto& face : hex_faces) {
        const Vec3& a = vertices[static_cast<std::size_t>(cell[static_cast<std::size_t>(face[0])])];
        const Vec3& b = vertices[static_cast<std::size_t>(cell[static_cast<std::size_t>(face[1])])];
        const Vec3& c = vertices[static_cast<std::size_t>(cell[static_cast<std::size_t>(face[2])])];
        const Vec3& d = vertices[static_cast<std::size_t>(cell[static_cast<std::size_t>(face[3])])];
        volume += dot(a, cross(b, c)) / 6.0;
        volume += dot(a, cross(c, d)) / 6.0;
    }
    return volume;
}

struct FaceCandidate {
    Quad vertices{};
    Index owner = invalid_index;
    Index owner_slot = invalid_index;
    Index neighbour = invalid_index;
    Index neighbour_slot = invalid_index;
};

}  // namespace

bool Mesh::isOwned(Index cell) const {
    return ownedIndex(cell) != invalid_index;
}

Index Mesh::ownedIndex(Index cell) const {
    if (cell < 0 || static_cast<std::size_t>(cell) >= m_storage.cell_owned_indices.size()) {
        throw std::out_of_range("cell index is outside the local mesh");
    }
    return m_storage.cell_owned_indices[static_cast<std::size_t>(cell)];
}

Index Mesh::globalCellId(Index cell) const {
    if (cell < 0 || static_cast<std::size_t>(cell) >= m_storage.cell_global_ids.size()) {
        throw std::out_of_range("cell index is outside the local mesh");
    }
    return m_storage.cell_global_ids[static_cast<std::size_t>(cell)];
}

Index Mesh::globalFaceId(Index face) const {
    if (face < 0 || static_cast<std::size_t>(face) >= m_storage.face_global_ids.size()) {
        throw std::out_of_range("face index is outside the local mesh");
    }
    return m_storage.face_global_ids[static_cast<std::size_t>(face)];
}

Index Mesh::cellOwnerRank(Index cell) const {
    if (cell < 0 || static_cast<std::size_t>(cell) >= m_storage.cell_owner_ranks.size()) {
        throw std::out_of_range("cell index is outside the local mesh");
    }
    return m_storage.cell_owner_ranks[static_cast<std::size_t>(cell)];
}

Index Mesh::faceOwnerRank(Index face) const {
    if (face < 0 || static_cast<std::size_t>(face) >= m_storage.face_owner_ranks.size()) {
        throw std::out_of_range("face index is outside the local mesh");
    }
    return m_storage.face_owner_ranks[static_cast<std::size_t>(face)];
}

void Mesh::setPartition(
    Index global_cells,
    Index layers,
    std::vector<Index> cell_ids,
    std::vector<Index> cell_owners,
    std::vector<Index> cell_depths,
    std::vector<Index> face_ids,
    std::vector<Index> face_owners,
    Index local_rank)
{
    if (global_cells <= 0 || layers < 0 || local_rank < 0 ||
        cell_ids.size() != static_cast<std::size_t>(cellCount()) ||
        cell_owners.size() != cell_ids.size() || cell_depths.size() != cell_ids.size() ||
        face_ids.size() != static_cast<std::size_t>(faceCount()) ||
        face_owners.size() != face_ids.size()) {
        throw std::invalid_argument("unstructured partition metadata is invalid");
    }
    std::set<Index> seen_cells;
    for (std::size_t cell = 0; cell < cell_ids.size(); ++cell) {
        if (cell_ids[cell] < 0 || cell_ids[cell] >= global_cells ||
            cell_owners[cell] < 0 || cell_depths[cell] < 0 ||
            !seen_cells.insert(cell_ids[cell]).second) {
            throw std::invalid_argument("unstructured partition cell metadata is invalid");
        }
    }
    std::set<Index> seen_faces;
    for (std::size_t face = 0; face < face_ids.size(); ++face) {
        if (face_ids[face] < 0 || face_owners[face] < 0 ||
            !seen_faces.insert(face_ids[face]).second) {
            throw std::invalid_argument("unstructured partition face metadata is invalid");
        }
    }

    m_storage.global_cell_count = global_cells;
    m_storage.ghost_layers = layers;
    m_storage.cell_global_ids.assign(std::move(cell_ids));
    m_storage.cell_owner_ranks.assign(std::move(cell_owners));
    m_storage.cell_ghost_depths.assign(std::move(cell_depths));
    m_storage.face_global_ids.assign(std::move(face_ids));
    m_storage.face_owner_ranks.assign(std::move(face_owners));
    m_storage.owned_cells.clear();
    m_storage.cell_owned_indices.assign(static_cast<std::size_t>(cellCount()), invalid_index);
    for (Index cell = 0; cell < cellCount(); ++cell) {
        const std::size_t c = static_cast<std::size_t>(cell);
        if (m_storage.cell_owner_ranks[c] == local_rank) {
            if (m_storage.cell_ghost_depths[c] != 0) {
                throw std::invalid_argument("owned unstructured cell has ghost depth");
            }
            m_storage.cell_owned_indices[c] = static_cast<Index>(m_storage.owned_cells.size());
            m_storage.owned_cells.push_back(cell);
        }
    }
    if (m_storage.owned_cells.empty()) {
        throw std::invalid_argument("partition owns no cells on this rank");
    }
    m_storage.owned_faces.clear();
    for (Index face = 0; face < faceCount(); ++face) {
        const std::size_t f = static_cast<std::size_t>(face);
        const Index owner = m_storage.face_owner[f];
        const Index neighbour = m_storage.face_neighbour[f];
        if (isOwned(owner) || (neighbour != invalid_index && isOwned(neighbour))) {
            m_storage.owned_faces.push_back(face);
        }
    }
}

Mesh Mesh::unstructured(
    std::vector<Vec3> vertices,
    std::vector<Hex> cells,
    std::vector<PatchSpec> patches,
    std::vector<BoundaryFaceSpec> boundary_faces)
{
    if (vertices.empty() || cells.empty() || patches.empty() ||
        vertices.size() > static_cast<std::size_t>(std::numeric_limits<Index>::max()) ||
        cells.size() > static_cast<std::size_t>(std::numeric_limits<Index>::max())) {
        throw std::invalid_argument("unstructured mesh arrays are empty or exceed 32-bit indexing");
    }
    if (!std::all_of(vertices.begin(), vertices.end(), isFinite)) {
        throw std::invalid_argument("unstructured mesh has a non-finite vertex");
    }
    for (const PatchSpec& patch : patches) {
        if (patch.name.empty()) throw std::invalid_argument("unstructured mesh has an unnamed patch");
    }

    Mesh mesh;
    mesh.m_storage.vertices.assign(std::move(vertices));
    mesh.m_storage.cell_vertices.assign(std::move(cells));
    const Index vertex_count = mesh.vertexCount();
    const Index cell_count = mesh.cellCount();
    mesh.m_storage.cell_centres.resize(static_cast<std::size_t>(cell_count));
    mesh.m_storage.cell_volumes.resize(static_cast<std::size_t>(cell_count));
    mesh.m_storage.cell_inverse_volumes.resize(static_cast<std::size_t>(cell_count));
    mesh.m_storage.cell_faces.resize(static_cast<std::size_t>(cell_count));
    mesh.m_storage.cell_neighbours.resize(static_cast<std::size_t>(cell_count));
    for (auto& values : mesh.m_storage.cell_faces) values.fill(invalid_index);
    for (auto& values : mesh.m_storage.cell_neighbours) values.fill(invalid_index);
    mesh.m_storage.patches.reserve(patches.size());
    for (PatchSpec& patch : patches) {
        mesh.m_storage.patches.push_back({std::move(patch.name), patch.kind, {}});
    }

    for (Index cell = 0; cell < cell_count; ++cell) {
        const Hex& vertices_of_cell = mesh.m_storage.cell_vertices[static_cast<std::size_t>(cell)];
        std::set<Index> unique;
        Vec3 centre{};
        for (Index vertex : vertices_of_cell) {
            requireIndex(vertex, vertex_count, "cell");
            unique.insert(vertex);
            centre += mesh.m_storage.vertices[static_cast<std::size_t>(vertex)];
        }
        if (unique.size() != vertices_of_cell.size()) {
            throw std::invalid_argument("unstructured cell repeats a vertex");
        }
        if (!(signedHexVolume(vertices_of_cell, mesh.m_storage.vertices) > 0.0)) {
            throw std::invalid_argument("unstructured cell has negative or degenerate orientation");
        }
        mesh.m_storage.cell_centres[static_cast<std::size_t>(cell)] = centre / 8.0;
    }

    std::map<Quad, Index> boundary_patch;
    for (const BoundaryFaceSpec& boundary : boundary_faces) {
        if (boundary.patch < 0 || boundary.patch >= mesh.patchCount()) {
            throw std::invalid_argument("boundary face references an invalid patch");
        }
        std::set<Index> unique;
        for (Index vertex : boundary.vertices) {
            requireIndex(vertex, vertex_count, "boundary face");
            unique.insert(vertex);
        }
        if (unique.size() != boundary.vertices.size() ||
            !boundary_patch.emplace(canonical(boundary.vertices), boundary.patch).second) {
            throw std::invalid_argument("unstructured mesh repeats a boundary face");
        }
    }

    std::map<Quad, Index> face_lookup;
    std::vector<FaceCandidate> candidates;
    candidates.reserve(static_cast<std::size_t>(cell_count) * 3U);
    for (Index cell = 0; cell < cell_count; ++cell) {
        const Hex& vertices_of_cell = mesh.m_storage.cell_vertices[static_cast<std::size_t>(cell)];
        for (Index slot = 0; slot < 6; ++slot) {
            Quad vertices_of_face{};
            for (Index local = 0; local < 4; ++local) {
                vertices_of_face[static_cast<std::size_t>(local)] =
                    vertices_of_cell[static_cast<std::size_t>(hex_faces[static_cast<std::size_t>(slot)]
                        [static_cast<std::size_t>(local)])];
            }
            const Quad key = canonical(vertices_of_face);
            const auto [position, inserted] = face_lookup.emplace(
                key, static_cast<Index>(candidates.size()));
            if (inserted) {
                candidates.push_back({vertices_of_face, cell, slot, invalid_index, invalid_index});
            } else {
                FaceCandidate& candidate = candidates[static_cast<std::size_t>(position->second)];
                if (candidate.neighbour != invalid_index) {
                    throw std::invalid_argument("unstructured mesh contains a non-manifold face");
                }
                candidate.neighbour = cell;
                candidate.neighbour_slot = slot;
            }
        }
    }

    for (const auto& entry : boundary_patch) {
        const auto face = face_lookup.find(entry.first);
        if (face == face_lookup.end()) {
            throw std::invalid_argument("boundary declaration does not match a cell face");
        }
        if (candidates[static_cast<std::size_t>(face->second)].neighbour != invalid_index) {
            throw std::invalid_argument("internal face cannot be assigned to a boundary patch");
        }
    }

    for (const FaceCandidate& candidate : candidates) {
        const Quad key = canonical(candidate.vertices);
        const auto boundary = boundary_patch.find(key);
        const Index patch = boundary == boundary_patch.end() ? invalid_index : boundary->second;
        if ((candidate.neighbour == invalid_index) != (patch != invalid_index)) {
            throw std::invalid_argument("every exterior face must be assigned to exactly one patch");
        }

        Quad face_vertices = candidate.vertices;
        auto [centre, area_vector] = quadGeometry(face_vertices, mesh.m_storage.vertices);
        const Vec3 direction = candidate.neighbour == invalid_index
            ? centre - mesh.m_storage.cell_centres[static_cast<std::size_t>(candidate.owner)]
            : mesh.m_storage.cell_centres[static_cast<std::size_t>(candidate.neighbour)] -
                mesh.m_storage.cell_centres[static_cast<std::size_t>(candidate.owner)];
        if (dot(area_vector, direction) < 0.0) {
            std::swap(face_vertices[1], face_vertices[3]);
            area_vector = -area_vector;
        }
        const Index face = mesh.faceCount();
        mesh.m_storage.face_vertices.push_back(face_vertices);
        mesh.m_storage.face_owner.push_back(candidate.owner);
        mesh.m_storage.face_neighbour.push_back(candidate.neighbour);
        mesh.m_storage.face_patch.push_back(patch);
        mesh.m_storage.face_centres.push_back(centre);
        mesh.m_storage.face_area_vectors.push_back(area_vector);
        mesh.m_storage.cell_faces[static_cast<std::size_t>(candidate.owner)]
            [static_cast<std::size_t>(candidate.owner_slot)] = face;
        if (candidate.neighbour != invalid_index) {
            mesh.m_storage.cell_faces[static_cast<std::size_t>(candidate.neighbour)]
                [static_cast<std::size_t>(candidate.neighbour_slot)] = face;
            mesh.m_storage.cell_neighbours[static_cast<std::size_t>(candidate.owner)]
                [static_cast<std::size_t>(candidate.owner_slot)] = candidate.neighbour;
            mesh.m_storage.cell_neighbours[static_cast<std::size_t>(candidate.neighbour)]
                [static_cast<std::size_t>(candidate.neighbour_slot)] = candidate.owner;
        } else {
            mesh.m_storage.patches[static_cast<std::size_t>(patch)].faces.push_back(face);
        }
    }

    for (Index cell = 0; cell < cell_count; ++cell) {
        const Vec3 reference = mesh.m_storage.cell_centres[static_cast<std::size_t>(cell)];
        double volume = 0.0;
        Vec3 first_moment{};
        for (Index face : mesh.m_storage.cell_faces[static_cast<std::size_t>(cell)]) {
            if (face == invalid_index) throw std::logic_error("cell face construction is incomplete");
            const Quad& ids = mesh.m_storage.face_vertices[static_cast<std::size_t>(face)];
            const Vec3& a = mesh.m_storage.vertices[static_cast<std::size_t>(ids[0])];
            const Vec3& b = mesh.m_storage.vertices[static_cast<std::size_t>(ids[1])];
            const Vec3& c = mesh.m_storage.vertices[static_cast<std::size_t>(ids[2])];
            const Vec3& d = mesh.m_storage.vertices[static_cast<std::size_t>(ids[3])];
            addTetrahedron(reference, a, b, c, volume, first_moment);
            addTetrahedron(reference, a, c, d, volume, first_moment);
        }
        if (!(volume > 0.0) || !std::isfinite(volume)) {
            throw std::invalid_argument("unstructured mesh contains a non-positive cell");
        }
        mesh.m_storage.cell_volumes[static_cast<std::size_t>(cell)] = volume;
        mesh.m_storage.cell_inverse_volumes[static_cast<std::size_t>(cell)] = 1.0 / volume;
        mesh.m_storage.cell_centres[static_cast<std::size_t>(cell)] = first_moment / volume;
    }

    const std::size_t face_count = mesh.m_storage.face_owner.size();
    mesh.m_storage.face_areas.resize(face_count);
    mesh.m_storage.face_normals.resize(face_count);
    mesh.m_storage.face_orthogonal_coefficients.resize(face_count);
    mesh.m_storage.face_owner_weights.resize(face_count);
    mesh.m_storage.face_non_orthogonal.resize(face_count);
    mesh.m_storage.face_skewness.resize(face_count);
    mesh.m_storage.orthogonal_geometry = true;
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const std::size_t f = static_cast<std::size_t>(face);
        const Index owner = mesh.m_storage.face_owner[f];
        const Index neighbour = mesh.m_storage.face_neighbour[f];
        const Vec3& centre = mesh.m_storage.face_centres[f];
        const Vec3& area_vector = mesh.m_storage.face_area_vectors[f];
        const Vec3 delta = neighbour == invalid_index
            ? centre - mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)]
            : mesh.m_storage.cell_centres[static_cast<std::size_t>(neighbour)] -
                mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)];
        const double delta_squared = squaredNorm(delta);
        const double area = norm(area_vector);
        const double projected = dot(area_vector, delta);
        if (!(delta_squared > 0.0) || !(area > 0.0) || !(projected > 0.0)) {
            throw std::invalid_argument("face orientation or geometry is invalid");
        }
        const double orthogonal = projected / delta_squared;
        mesh.m_storage.face_areas[f] = area;
        mesh.m_storage.face_normals[f] = area_vector / area;
        mesh.m_storage.face_orthogonal_coefficients[f] = orthogonal;
        mesh.m_storage.face_non_orthogonal[f] = area_vector - orthogonal * delta;
        if (norm(mesh.m_storage.face_non_orthogonal[f]) > 1e-12 * area) {
            mesh.m_storage.orthogonal_geometry = false;
        }
        if (neighbour == invalid_index) {
            mesh.m_storage.face_owner_weights[f] = 1.0;
            mesh.m_storage.face_skewness[f] = {};
        } else {
            const Vec3 normal = mesh.m_storage.face_normals[f];
            const double owner_distance = std::abs(dot(
                centre - mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)], normal));
            const double neighbour_distance = std::abs(dot(
                mesh.m_storage.cell_centres[static_cast<std::size_t>(neighbour)] - centre, normal));
            const double distance_sum = owner_distance + neighbour_distance;
            if (!(distance_sum > 0.0)) {
                throw std::invalid_argument("face interpolation distance is invalid");
            }
            mesh.m_storage.face_owner_weights[f] = neighbour_distance / distance_sum;
            const double denominator = dot(delta, normal);
            const double fraction = dot(centre -
                mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)], normal) / denominator;
            mesh.m_storage.face_skewness[f] = centre -
                (mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)] + fraction * delta);
            if (norm(mesh.m_storage.face_skewness[f]) > 1e-12 * std::sqrt(delta_squared)) {
                mesh.m_storage.orthogonal_geometry = false;
            }
        }
    }

    std::vector<Index> cell_ids(static_cast<std::size_t>(cell_count));
    std::vector<Index> cell_owners(static_cast<std::size_t>(cell_count), 0);
    std::vector<Index> cell_depths(static_cast<std::size_t>(cell_count), 0);
    std::vector<Index> face_ids(static_cast<std::size_t>(mesh.faceCount()));
    std::vector<Index> face_owners(static_cast<std::size_t>(mesh.faceCount()), 0);
    for (Index cell = 0; cell < cell_count; ++cell) cell_ids[static_cast<std::size_t>(cell)] = cell;
    for (Index face = 0; face < mesh.faceCount(); ++face) face_ids[static_cast<std::size_t>(face)] = face;
    mesh.setPartition(cell_count, 0, std::move(cell_ids), std::move(cell_owners),
                      std::move(cell_depths), std::move(face_ids), std::move(face_owners), 0);
    mesh.validate();
    return mesh;
}

void Mesh::validate() const {
    const Index cells = cellCount();
    const Index faces = faceCount();
    if (cells <= 0 || faces <= 0 || vertexCount() <= 0 || m_storage.global_cell_count <= 0 ||
        m_storage.cell_centres.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_volumes.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_inverse_volumes.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_faces.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_neighbours.size() != static_cast<std::size_t>(cells) ||
        m_storage.face_vertices.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_neighbour.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_patch.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_centres.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_area_vectors.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_normals.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_areas.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_orthogonal_coefficients.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_owner_weights.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_non_orthogonal.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_skewness.size() != static_cast<std::size_t>(faces) ||
        m_storage.cell_global_ids.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_owner_ranks.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_ghost_depths.size() != static_cast<std::size_t>(cells) ||
        m_storage.face_global_ids.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_owner_ranks.size() != static_cast<std::size_t>(faces) ||
        m_storage.cell_owned_indices.size() != static_cast<std::size_t>(cells) ||
        m_storage.patches.empty()) {
        throw std::runtime_error("unstructured mesh storage is inconsistent");
    }
    std::set<Index> cell_ids;
    std::set<Index> face_ids;
    std::vector<bool> patch_seen(static_cast<std::size_t>(faces), false);
    for (Index cell = 0; cell < cells; ++cell) {
        const std::size_t c = static_cast<std::size_t>(cell);
        const Hex& vertices = m_storage.cell_vertices[c];
        std::set<Index> unique(vertices.begin(), vertices.end());
        if (unique.size() != vertices.size() || !cell_ids.insert(m_storage.cell_global_ids[c]).second ||
            m_storage.cell_global_ids[c] < 0 || m_storage.cell_global_ids[c] >= m_storage.global_cell_count ||
            m_storage.cell_owner_ranks[c] < 0 || m_storage.cell_ghost_depths[c] < 0 ||
            !(m_storage.cell_volumes[c] > 0.0) || !std::isfinite(m_storage.cell_volumes[c])) {
            throw std::runtime_error("unstructured cell metadata is invalid");
        }
        for (Index vertex : vertices) {
            if (vertex < 0 || vertex >= vertexCount()) throw std::runtime_error("cell vertex is invalid");
        }
        for (Index slot = 0; slot < 6; ++slot) {
            const Index face = m_storage.cell_faces[c][static_cast<std::size_t>(slot)];
            const Index neighbour = m_storage.cell_neighbours[c][static_cast<std::size_t>(slot)];
            if (face < 0 || face >= faces || (neighbour != invalid_index &&
                (neighbour < 0 || neighbour >= cells))) {
                throw std::runtime_error("cell connectivity is invalid");
            }
        }
    }
    for (Index face = 0; face < faces; ++face) {
        const std::size_t f = static_cast<std::size_t>(face);
        const Index owner = m_storage.face_owner[f];
        const Index neighbour = m_storage.face_neighbour[f];
        const Index patch = m_storage.face_patch[f];
        if (owner < 0 || owner >= cells || (neighbour != invalid_index &&
            (neighbour < 0 || neighbour >= cells || neighbour == owner)) ||
            m_storage.face_owner_ranks[f] < 0 || !face_ids.insert(m_storage.face_global_ids[f]).second ||
            m_storage.face_global_ids[f] < 0 || (neighbour == invalid_index) != (patch != invalid_index) ||
            (patch != invalid_index && (patch < 0 || patch >= patchCount())) ||
            !(m_storage.face_areas[f] > 0.0) || !std::isfinite(m_storage.face_areas[f])) {
            throw std::runtime_error("face metadata is invalid");
        }
    }
    for (Index patch = 0; patch < patchCount(); ++patch) {
        const BoundaryPatch& boundary = m_storage.patches[static_cast<std::size_t>(patch)];
        if (boundary.name.empty()) throw std::runtime_error("mesh has an unnamed patch");
        for (Index face : boundary.faces) {
            if (face < 0 || face >= faces || m_storage.face_patch[static_cast<std::size_t>(face)] != patch ||
                patch_seen[static_cast<std::size_t>(face)]) {
                throw std::runtime_error("boundary patch connectivity is invalid");
            }
            patch_seen[static_cast<std::size_t>(face)] = true;
        }
    }
    for (Index face = 0; face < faces; ++face) {
        const std::size_t f = static_cast<std::size_t>(face);
        if ((m_storage.face_patch[f] != invalid_index) != patch_seen[f]) {
            throw std::runtime_error("boundary face is missing from its patch");
        }
    }
    Index owned_count = 0;
    for (Index cell = 0; cell < cells; ++cell) {
        const Index owned = m_storage.cell_owned_indices[static_cast<std::size_t>(cell)];
        if (owned != invalid_index) {
            if (owned != owned_count++ || m_storage.cell_ghost_depths[static_cast<std::size_t>(cell)] != 0) {
                throw std::runtime_error("owned cell mapping is invalid");
            }
        }
    }
    if (owned_count != ownedCellCount()) throw std::runtime_error("owned cell count is invalid");
}

}  // namespace babelsim
