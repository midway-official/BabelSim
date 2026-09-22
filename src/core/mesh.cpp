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

void requireIndex(Index value, Index count, const char* what) {
    if (value < 0 || value >= count) {
        throw std::invalid_argument(std::string(what) + " references an invalid vertex");
    }
}

}  // namespace

namespace {

struct PolygonTriangle {
    Index a;
    Index b;
    Index c;
};

struct ProjectedPoint {
    double x;
    double y;
};

ProjectedPoint project(const Vec3& point, int dropped_axis) {
    if (dropped_axis == 0) return {point.y, point.z};
    if (dropped_axis == 1) return {point.z, point.x};
    return {point.x, point.y};
}

double cross2d(const ProjectedPoint& a, const ProjectedPoint& b, const ProjectedPoint& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool insideOrOnTriangle(
    const ProjectedPoint& point,
    const ProjectedPoint& a,
    const ProjectedPoint& b,
    const ProjectedPoint& c,
    double orientation,
    double tolerance)
{
    const double ab = orientation * cross2d(a, b, point);
    const double bc = orientation * cross2d(b, c, point);
    const double ca = orientation * cross2d(c, a, point);
    return ab >= -tolerance && bc >= -tolerance && ca >= -tolerance;
}

std::vector<PolygonTriangle> polygonTriangles(
    const std::vector<Index>& ids, const MeshStorage<Vec3>& vertices)
{
    if (ids.size() < 3) throw std::invalid_argument("polyhedral face has fewer than three vertices");

    Vec3 normal{};
    double minimum[3]{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    double maximum[3]{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()};
    for (std::size_t i = 0; i < ids.size(); ++i) {
        const Vec3& current = vertices[static_cast<std::size_t>(ids[i])];
        const Vec3& next = vertices[static_cast<std::size_t>(ids[(i + 1U) % ids.size()])];
        normal += 0.5 * cross(current, next);
        minimum[0] = std::min(minimum[0], current.x);
        minimum[1] = std::min(minimum[1], current.y);
        minimum[2] = std::min(minimum[2], current.z);
        maximum[0] = std::max(maximum[0], current.x);
        maximum[1] = std::max(maximum[1], current.y);
        maximum[2] = std::max(maximum[2], current.z);
    }
    const double normal_components[3]{std::abs(normal.x), std::abs(normal.y), std::abs(normal.z)};
    const int dropped_axis = normal_components[1] > normal_components[0]
        ? (normal_components[2] > normal_components[1] ? 2 : 1)
        : (normal_components[2] > normal_components[0] ? 2 : 0);
    const double scale = std::max({
        maximum[0] - minimum[0], maximum[1] - minimum[1], maximum[2] - minimum[2], 1.0});
    const double tolerance = 1e-13 * scale * scale;
    std::vector<ProjectedPoint> projected;
    projected.reserve(ids.size());
    for (Index id : ids) projected.push_back(project(vertices[static_cast<std::size_t>(id)], dropped_axis));

    double signed_area = 0.0;
    for (std::size_t i = 0; i < projected.size(); ++i) {
        const auto& current = projected[i];
        const auto& next = projected[(i + 1U) % projected.size()];
        signed_area += current.x * next.y - next.x * current.y;
    }
    if (std::abs(signed_area) <= tolerance) {
        throw std::invalid_argument("polyhedral face projects to zero area");
    }
    const double orientation = signed_area > 0.0 ? 1.0 : -1.0;

    bool convex = true;
    bool saw_turn = false;
    for (std::size_t i = 0; i < projected.size(); ++i) {
        const double turn = orientation * cross2d(
            projected[(i + projected.size() - 1U) % projected.size()],
            projected[i], projected[(i + 1U) % projected.size()]);
        saw_turn = saw_turn || std::abs(turn) > tolerance;
        if (turn < -tolerance) {
            convex = false;
            break;
        }
    }
    if (!saw_turn) throw std::invalid_argument("polyhedral face is collinear");
    if (convex) {
        std::vector<PolygonTriangle> triangles;
        triangles.reserve(ids.size() - 2U);
        for (Index i = 1; i + 1 < static_cast<Index>(ids.size()); ++i)
            triangles.push_back({ids[0], ids[static_cast<std::size_t>(i)], ids[static_cast<std::size_t>(i + 1)]});
        return triangles;
    }

    // Ear clipping is only needed for concave rings.  Convex rings retain the
    // historical fan triangulation so old warped-hex geometry remains stable.
    std::vector<std::size_t> remaining(ids.size());
    for (std::size_t i = 0; i < remaining.size(); ++i) remaining[i] = i;
    std::vector<PolygonTriangle> triangles;
    triangles.reserve(ids.size() - 2U);
    while (remaining.size() > 3U) {
        bool clipped = false;
        for (std::size_t position = 0; position < remaining.size(); ++position) {
            const std::size_t previous = remaining[(position + remaining.size() - 1U) % remaining.size()];
            const std::size_t current = remaining[position];
            const std::size_t next = remaining[(position + 1U) % remaining.size()];
            if (orientation * cross2d(projected[previous], projected[current], projected[next]) <= tolerance)
                continue;
            bool contains_vertex = false;
            for (std::size_t candidate : remaining) {
                if (candidate == previous || candidate == current || candidate == next) continue;
                if (insideOrOnTriangle(projected[candidate], projected[previous], projected[current],
                                       projected[next], orientation, tolerance)) {
                    contains_vertex = true;
                    break;
                }
            }
            if (contains_vertex) continue;
            triangles.push_back({ids[previous], ids[current], ids[next]});
            remaining.erase(remaining.begin() + static_cast<std::ptrdiff_t>(position));
            clipped = true;
            break;
        }
        if (!clipped) throw std::invalid_argument("polyhedral face is self-intersecting or not simple");
    }
    triangles.push_back({ids[remaining[0]], ids[remaining[1]], ids[remaining[2]]});
    return triangles;
}

Vec3 polygonAreaVector(const std::vector<Index>& ids, const MeshStorage<Vec3>& vertices) {
    Vec3 area{};
    for (const PolygonTriangle& triangle : polygonTriangles(ids, vertices)) {
        const Vec3& a = vertices[static_cast<std::size_t>(triangle.a)];
        const Vec3& b = vertices[static_cast<std::size_t>(triangle.b)];
        const Vec3& c = vertices[static_cast<std::size_t>(triangle.c)];
        area += 0.5 * cross(b - a, c - a);
    }
    return area;
}

std::pair<Vec3, double> polygonGeometry(
    const std::vector<Index>& ids, const MeshStorage<Vec3>& vertices)
{
    Vec3 centre{};
    double total = 0.0;
    for (const PolygonTriangle& triangle : polygonTriangles(ids, vertices)) {
        const Vec3& a = vertices[static_cast<std::size_t>(triangle.a)];
        const Vec3& b = vertices[static_cast<std::size_t>(triangle.b)];
        const Vec3& c = vertices[static_cast<std::size_t>(triangle.c)];
        const Vec3 area = 0.5 * cross(b - a, c - a);
        const double weight = norm(area);
        if (!(weight > 0.0) || !std::isfinite(weight)) {
            throw std::invalid_argument("polyhedral mesh contains a degenerate face triangle");
        }
        centre += weight * ((a + b + c) / 3.0);
        total += weight;
    }
    if (!(total > 0.0) || !std::isfinite(total)) {
        throw std::invalid_argument("polyhedral mesh contains a degenerate face");
    }
    return {centre / total, total};
}

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

MeshPartitionInfo Mesh::partitionInfo() const {
    MeshPartitionInfo result;
    result.global_cells = m_storage.global_cell_count;
    result.local_cells = cellCount();
    result.owned_cells = static_cast<Index>(m_storage.owned_cells.size());
    result.ghost_cells = result.local_cells - result.owned_cells;
    result.local_faces = faceCount();
    result.owned_faces = static_cast<Index>(m_storage.owned_faces.size());
    for (Index patch = 0; patch < patchCount(); ++patch) {
        const auto& boundary = m_storage.patches[static_cast<std::size_t>(patch)];
        if (boundary.kind != PatchKind::Processor || boundary.faces.empty()) continue;
        ++result.neighbour_ranks;
        result.communication_faces += static_cast<Index>(boundary.faces.size());
    }
    return result;
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

Mesh Mesh::polyhedral(
    std::vector<Vec3> vertices,
    std::vector<PolyhedralFaceSpec> faces,
    std::vector<PatchSpec> patches)
{
    if (vertices.empty() || faces.empty() || patches.empty() ||
        vertices.size() > static_cast<std::size_t>(std::numeric_limits<Index>::max()) ||
        faces.size() > static_cast<std::size_t>(std::numeric_limits<Index>::max())) {
        throw std::invalid_argument("polyhedral mesh arrays are empty or exceed 32-bit indexing");
    }
    if (!std::all_of(vertices.begin(), vertices.end(), isFinite)) {
        throw std::invalid_argument("polyhedral mesh has a non-finite vertex");
    }
    for (const PatchSpec& patch : patches) {
        if (patch.name.empty()) throw std::invalid_argument("polyhedral mesh has an unnamed patch");
    }

    Index cell_count = 0;
    for (const PolyhedralFaceSpec& face : faces) {
        if (face.owner < 0 || face.neighbour < invalid_index ||
            (face.neighbour != invalid_index && face.neighbour == face.owner)) {
            throw std::invalid_argument("polyhedral face has invalid owner/neighbour");
        }
        cell_count = std::max(cell_count, face.owner + 1);
        if (face.neighbour != invalid_index) cell_count = std::max(cell_count, face.neighbour + 1);
    }
    if (cell_count <= 0) throw std::invalid_argument("polyhedral mesh has no cells");

    Mesh mesh;
    mesh.m_storage.vertices.assign(std::move(vertices));
    mesh.m_storage.patches.reserve(patches.size());
    for (PatchSpec& patch : patches) {
        mesh.m_storage.patches.push_back({std::move(patch.name), patch.kind, {}});
    }

    std::vector<Vec3> cell_point_sums(static_cast<std::size_t>(cell_count));
    std::vector<Index> cell_point_counts(static_cast<std::size_t>(cell_count), 0);
    std::vector<std::vector<Index>> cell_faces(static_cast<std::size_t>(cell_count));
    std::vector<std::vector<Index>> cell_signs(static_cast<std::size_t>(cell_count));
    std::vector<std::set<Index>> cell_neighbours(static_cast<std::size_t>(cell_count));

    // Establish a stable orientation reference before processing any face.  A
    // running average would make the orientation depend on input face order.
    for (const PolyhedralFaceSpec& face : faces) {
        Vec3 sum{};
        for (Index vertex : face.vertices) {
            requireIndex(vertex, mesh.vertexCount(), "face");
            sum += mesh.m_storage.vertices[static_cast<std::size_t>(vertex)];
        }
        cell_point_sums[static_cast<std::size_t>(face.owner)] += sum;
        cell_point_counts[static_cast<std::size_t>(face.owner)] +=
            static_cast<Index>(face.vertices.size());
        if (face.neighbour != invalid_index) {
            cell_point_sums[static_cast<std::size_t>(face.neighbour)] += sum;
            cell_point_counts[static_cast<std::size_t>(face.neighbour)] +=
                static_cast<Index>(face.vertices.size());
        }
    }

    mesh.m_storage.face_point_offsets.push_back(0);
    mesh.m_storage.face_neighbour.reserve(faces.size());
    mesh.m_storage.face_owner.reserve(faces.size());
    mesh.m_storage.face_patch.reserve(faces.size());
    mesh.m_storage.face_centres.reserve(faces.size());
    mesh.m_storage.face_area_vectors.reserve(faces.size());

    for (std::size_t input_face = 0; input_face < faces.size(); ++input_face) {
        PolyhedralFaceSpec face = std::move(faces[input_face]);
        if (face.vertices.size() < 3) {
            throw std::invalid_argument("polyhedral face has fewer than three vertices");
        }
        std::set<Index> unique;
        Vec3 face_point_sum{};
        for (Index vertex : face.vertices) {
            requireIndex(vertex, mesh.vertexCount(), "face");
            if (!unique.insert(vertex).second) {
                throw std::invalid_argument("polyhedral face repeats a vertex");
            }
            face_point_sum += mesh.m_storage.vertices[static_cast<std::size_t>(vertex)];
        }
        const Vec3 point_average = face_point_sum / static_cast<double>(face.vertices.size());
        cell_faces[static_cast<std::size_t>(face.owner)].push_back(static_cast<Index>(input_face));
        cell_signs[static_cast<std::size_t>(face.owner)].push_back(1);
        if (face.neighbour != invalid_index) {
            cell_faces[static_cast<std::size_t>(face.neighbour)].push_back(static_cast<Index>(input_face));
            cell_signs[static_cast<std::size_t>(face.neighbour)].push_back(-1);
            cell_neighbours[static_cast<std::size_t>(face.owner)].insert(face.neighbour);
            cell_neighbours[static_cast<std::size_t>(face.neighbour)].insert(face.owner);
            if (face.patch != invalid_index) {
                throw std::invalid_argument("internal polyhedral face cannot belong to a boundary patch");
            }
        } else if (face.patch < 0 || face.patch >= mesh.patchCount()) {
            throw std::invalid_argument("boundary polyhedral face " + std::to_string(input_face) +
                                        " has invalid patch " +
                                        std::to_string(face.patch) + " of " +
                                        std::to_string(mesh.patchCount()));
        }

        // Orient the ring from the owner toward the exterior.  The arithmetic
        // average is only used to choose orientation; the final cell centre is
        // computed from signed face tetrahedra below.
        Vec3 area = polygonAreaVector(face.vertices, mesh.m_storage.vertices);
        // A provisional owner centre based on all owner face points gives a
        // deterministic outward direction for non-centred input rings.
        const Vec3 owner_estimate = cell_point_sums[static_cast<std::size_t>(face.owner)] /
            static_cast<double>(std::max<Index>(1, cell_point_counts[static_cast<std::size_t>(face.owner)]));
        if (dot(area, point_average - owner_estimate) < 0.0) {
            std::reverse(face.vertices.begin() + 1, face.vertices.end());
            area = polygonAreaVector(face.vertices, mesh.m_storage.vertices);
        }
        if (!(norm(area) > 0.0) || !isFinite(area)) {
            throw std::invalid_argument("polyhedral face has zero area");
        }
        auto [centre, area_size] = polygonGeometry(face.vertices, mesh.m_storage.vertices);
        for (Index vertex : face.vertices) mesh.m_storage.face_point_ids.push_back(vertex);
        mesh.m_storage.face_point_offsets.push_back(
            static_cast<Index>(mesh.m_storage.face_point_ids.size()));
        mesh.m_storage.face_owner.push_back(face.owner);
        mesh.m_storage.face_neighbour.push_back(face.neighbour);
        mesh.m_storage.face_patch.push_back(face.patch);
        mesh.m_storage.face_centres.push_back(centre);
        mesh.m_storage.face_area_vectors.push_back(area);
        mesh.m_storage.face_areas.push_back(area_size);
        mesh.m_storage.face_normals.push_back(area / area_size);
        mesh.m_storage.face_owner_weights.push_back(1.0);
        mesh.m_storage.face_non_orthogonal.push_back({});
        mesh.m_storage.face_skewness.push_back({});
        mesh.m_storage.face_orthogonal_coefficients.push_back(0.0);
        if (face.patch != invalid_index) {
            mesh.m_storage.patches[static_cast<std::size_t>(face.patch)].faces.push_back(
                static_cast<Index>(input_face));
        }
    }

    mesh.m_storage.cell_face_offsets.push_back(0);
    for (Index cell = 0; cell < cell_count; ++cell) {
        for (std::size_t i = 0; i < cell_faces[static_cast<std::size_t>(cell)].size(); ++i) {
            mesh.m_storage.cell_face_ids.push_back(cell_faces[static_cast<std::size_t>(cell)][i]);
            mesh.m_storage.cell_face_signs.push_back(cell_signs[static_cast<std::size_t>(cell)][i]);
        }
        mesh.m_storage.cell_face_offsets.push_back(
            static_cast<Index>(mesh.m_storage.cell_face_ids.size()));
        for (Index neighbour : cell_neighbours[static_cast<std::size_t>(cell)]) {
            mesh.m_storage.cell_neighbour_ids.push_back(neighbour);
        }
        if (mesh.m_storage.cell_neighbour_offsets.empty()) {
            mesh.m_storage.cell_neighbour_offsets.push_back(0);
        }
        mesh.m_storage.cell_neighbour_offsets.push_back(
            static_cast<Index>(mesh.m_storage.cell_neighbour_ids.size()));
    }

    mesh.m_storage.cell_vertex_offsets.clear();
    mesh.m_storage.cell_vertex_ids.clear();
    mesh.m_storage.cell_vertex_offsets.push_back(0);
    for (Index cell = 0; cell < cell_count; ++cell) {
        std::set<Index> unique;
        for (Index face : mesh.cellFaces(cell)) {
            for (Index vertex : mesh.facePoints(face)) unique.insert(vertex);
        }
        for (Index vertex : unique) mesh.m_storage.cell_vertex_ids.push_back(vertex);
        mesh.m_storage.cell_vertex_offsets.push_back(
            static_cast<Index>(mesh.m_storage.cell_vertex_ids.size()));
    }

    mesh.m_storage.cell_centres.resize(static_cast<std::size_t>(cell_count));
    mesh.m_storage.cell_volumes.resize(static_cast<std::size_t>(cell_count));
    mesh.m_storage.cell_inverse_volumes.resize(static_cast<std::size_t>(cell_count));
    for (Index cell = 0; cell < cell_count; ++cell) {
        const std::size_t c = static_cast<std::size_t>(cell);
        if (cell_faces[c].empty()) throw std::invalid_argument("polyhedral cell has no faces");
        const Vec3 reference = cell_point_sums[c] /
            static_cast<double>(std::max<Index>(1, cell_point_counts[c]));
        double volume = 0.0;
        Vec3 moment{};
        for (Index face : cell_faces[c]) {
            const std::size_t f = static_cast<std::size_t>(face);
            const Index sign = mesh.m_storage.face_owner[f] == cell ? 1 : -1;
            const Index first = mesh.m_storage.face_point_offsets[f];
            const Index last = mesh.m_storage.face_point_offsets[f + 1U];
            std::vector<Index> ids;
            ids.reserve(static_cast<std::size_t>(last - first));
            for (Index i = first; i < last; ++i)
                ids.push_back(mesh.m_storage.face_point_ids[static_cast<std::size_t>(i)]);
            for (const PolygonTriangle& triangle : polygonTriangles(ids, mesh.m_storage.vertices)) {
                const Vec3& a = mesh.m_storage.vertices[static_cast<std::size_t>(triangle.a)];
                const Vec3& b = mesh.m_storage.vertices[static_cast<std::size_t>(triangle.b)];
                const Vec3& d = mesh.m_storage.vertices[static_cast<std::size_t>(triangle.c)];
                const double tetra = sign * dot(a - reference, cross(b - reference, d - reference)) / 6.0;
                volume += tetra;
                moment += tetra * ((reference + a + b + d) / 4.0);
            }
        }
        if (!(volume > 0.0) || !std::isfinite(volume)) {
            throw std::invalid_argument("polyhedral cell has non-positive signed volume");
        }
        mesh.m_storage.cell_volumes[c] = volume;
        mesh.m_storage.cell_inverse_volumes[c] = 1.0 / volume;
        mesh.m_storage.cell_centres[c] = moment / volume;
    }

    mesh.m_storage.orthogonal_geometry = true;
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const std::size_t f = static_cast<std::size_t>(face);
        const Index owner = mesh.m_storage.face_owner[f];
        const Index neighbour = mesh.m_storage.face_neighbour[f];
        const Vec3 delta = neighbour == invalid_index
            ? mesh.m_storage.face_centres[f] - mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)]
            : mesh.m_storage.cell_centres[static_cast<std::size_t>(neighbour)] -
                mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)];
        const double delta_squared = squaredNorm(delta);
        const double projected = dot(mesh.m_storage.face_area_vectors[f], delta);
        if (!(delta_squared > 0.0) || !(projected > 0.0)) {
            throw std::invalid_argument("polyhedral face geometry has invalid owner/neighbour direction");
        }
        const Vec3 normal = mesh.m_storage.face_normals[f];
        mesh.m_storage.face_orthogonal_coefficients[f] = projected / delta_squared;
        mesh.m_storage.face_non_orthogonal[f] =
            mesh.m_storage.face_area_vectors[f] -
            mesh.m_storage.face_orthogonal_coefficients[f] * delta;
        if (norm(mesh.m_storage.face_non_orthogonal[f]) > 1e-12 * mesh.m_storage.face_areas[f]) {
            mesh.m_storage.orthogonal_geometry = false;
        }
        if (neighbour != invalid_index) {
            const double owner_distance = std::abs(dot(
                mesh.m_storage.face_centres[f] - mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)], normal));
            const double neighbour_distance = std::abs(dot(
                mesh.m_storage.cell_centres[static_cast<std::size_t>(neighbour)] - mesh.m_storage.face_centres[f], normal));
            const double distance_sum = owner_distance + neighbour_distance;
            if (!(distance_sum > 0.0)) throw std::invalid_argument("polyhedral interpolation distance is invalid");
            mesh.m_storage.face_owner_weights[f] = neighbour_distance / distance_sum;
            const double denominator = dot(delta, normal);
            if (!(denominator > 0.0) || !std::isfinite(denominator)) {
                throw std::invalid_argument("polyhedral face interpolation direction is invalid");
            }
            const double fraction = dot(
                mesh.m_storage.face_centres[f] - mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)],
                normal) / denominator;
            mesh.m_storage.face_skewness[f] = mesh.m_storage.face_centres[f] -
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

Mesh Mesh::polyhedral(
    std::vector<Vec3> vertices,
    std::vector<std::vector<Index>> face_vertices,
    std::vector<Index> owners,
    std::vector<Index> neighbours,
    std::vector<Index> face_patches,
    std::vector<PatchSpec> patches)
{
    if (face_vertices.size() != owners.size() || owners.size() != neighbours.size() ||
        owners.size() != face_patches.size()) {
        throw std::invalid_argument("polyhedral face arrays have mismatched lengths");
    }
    std::vector<PolyhedralFaceSpec> faces;
    faces.reserve(face_vertices.size());
    for (std::size_t i = 0; i < face_vertices.size(); ++i) {
        faces.push_back({std::move(face_vertices[i]), owners[i], neighbours[i], face_patches[i]});
    }
    return polyhedral(std::move(vertices), std::move(faces), std::move(patches));
}

void Mesh::validate() const
{
    const Index cells = cellCount();
    const Index faces = faceCount();
    if (cells <= 0 || faces <= 0 || vertexCount() <= 0 || patchCount() <= 0 ||
        m_storage.global_cell_count <= 0 ||
        m_storage.face_point_offsets.size() != static_cast<std::size_t>(faces) + 1U ||
        m_storage.cell_face_offsets.size() != static_cast<std::size_t>(cells) + 1U ||
        m_storage.cell_face_signs.size() != m_storage.cell_face_ids.size() ||
        m_storage.cell_vertex_offsets.size() != static_cast<std::size_t>(cells) + 1U ||
        m_storage.cell_neighbour_offsets.size() != static_cast<std::size_t>(cells) + 1U ||
        m_storage.cell_centres.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_volumes.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_inverse_volumes.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_global_ids.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_owner_ranks.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_ghost_depths.size() != static_cast<std::size_t>(cells) ||
        m_storage.cell_owned_indices.size() != static_cast<std::size_t>(cells) ||
        m_storage.face_global_ids.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_owner_ranks.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_neighbour.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_patch.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_centres.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_area_vectors.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_normals.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_areas.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_orthogonal_coefficients.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_owner_weights.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_non_orthogonal.size() != static_cast<std::size_t>(faces) ||
        m_storage.face_skewness.size() != static_cast<std::size_t>(faces)) {
        throw std::runtime_error("polyhedral mesh storage is inconsistent");
    }

    std::set<Index> cell_ids;
    for (Index cell = 0; cell < cells; ++cell) {
        const std::size_t c = static_cast<std::size_t>(cell);
        if (m_storage.cell_global_ids[c] < 0 ||
            m_storage.cell_global_ids[c] >= m_storage.global_cell_count ||
            !cell_ids.insert(m_storage.cell_global_ids[c]).second ||
            m_storage.cell_owner_ranks[c] < 0 || m_storage.cell_ghost_depths[c] < 0 ||
            !(m_storage.cell_volumes[c] > 0.0) || !std::isfinite(m_storage.cell_volumes[c])) {
            throw std::runtime_error("polyhedral cell metadata is invalid");
        }
        const Index face_first = m_storage.cell_face_offsets[c];
        const Index face_last = m_storage.cell_face_offsets[c + 1U];
        const Index vertex_first = m_storage.cell_vertex_offsets[c];
        const Index vertex_last = m_storage.cell_vertex_offsets[c + 1U];
        const Index neighbour_first = m_storage.cell_neighbour_offsets[c];
        const Index neighbour_last = m_storage.cell_neighbour_offsets[c + 1U];
        if (face_first < 0 || face_last <= face_first ||
            face_last > static_cast<Index>(m_storage.cell_face_ids.size()) ||
            vertex_first < 0 || vertex_last <= vertex_first ||
            vertex_last > static_cast<Index>(m_storage.cell_vertex_ids.size()) ||
            neighbour_first < 0 || neighbour_last < neighbour_first ||
            neighbour_last > static_cast<Index>(m_storage.cell_neighbour_ids.size())) {
            throw std::runtime_error("polyhedral cell connectivity is invalid");
        }
        for (Index i = face_first; i < face_last; ++i) {
            const Index face = m_storage.cell_face_ids[static_cast<std::size_t>(i)];
            const Index sign = m_storage.cell_face_signs[static_cast<std::size_t>(i)];
            if (face < 0 || face >= faces || (sign != 1 && sign != -1)) {
                throw std::runtime_error("polyhedral cell-face sign is invalid");
            }
        }
        for (Index i = vertex_first; i < vertex_last; ++i) {
            const Index vertex = m_storage.cell_vertex_ids[static_cast<std::size_t>(i)];
            if (vertex < 0 || vertex >= vertexCount()) {
                throw std::runtime_error("polyhedral cell vertex is invalid");
            }
        }
        for (Index i = neighbour_first; i < neighbour_last; ++i) {
            const Index neighbour = m_storage.cell_neighbour_ids[static_cast<std::size_t>(i)];
            if (neighbour < 0 || neighbour >= cells || neighbour == cell) {
                throw std::runtime_error("polyhedral cell neighbour is invalid");
            }
        }
    }

    std::set<Index> face_ids;
    std::vector<bool> patch_seen(static_cast<std::size_t>(faces), false);
    for (Index face = 0; face < faces; ++face) {
        const std::size_t f = static_cast<std::size_t>(face);
        const Index first = m_storage.face_point_offsets[f];
        const Index last = m_storage.face_point_offsets[f + 1U];
        const Index owner = m_storage.face_owner[f];
        const Index neighbour = m_storage.face_neighbour[f];
        const Index patch = m_storage.face_patch[f];
        if (first < 0 || last - first < 3 ||
            last > static_cast<Index>(m_storage.face_point_ids.size()) ||
            owner < 0 || owner >= cells ||
            (neighbour != invalid_index &&
             (neighbour < 0 || neighbour >= cells || neighbour == owner)) ||
            (neighbour == invalid_index) != (patch != invalid_index) ||
            (patch != invalid_index && (patch < 0 || patch >= patchCount())) ||
            !face_ids.insert(m_storage.face_global_ids[f]).second ||
            m_storage.face_global_ids[f] < 0 ||
            m_storage.face_owner_ranks[f] < 0 ||
            !(m_storage.face_areas[f] > 0.0) || !std::isfinite(m_storage.face_areas[f])) {
            throw std::runtime_error("polyhedral face metadata is invalid");
        }
        std::set<Index> unique;
        for (Index i = first; i < last; ++i) {
            const Index vertex = m_storage.face_point_ids[static_cast<std::size_t>(i)];
            if (vertex < 0 || vertex >= vertexCount() || !unique.insert(vertex).second) {
                throw std::runtime_error("polyhedral face vertex connectivity is invalid");
            }
        }
    }

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
        const Index first = m_storage.cell_face_offsets[c];
        const Index last = m_storage.cell_face_offsets[c + 1U];
        for (Index position = first; position < last; ++position) {
            const Index face = m_storage.cell_face_ids[static_cast<std::size_t>(position)];
            const Index sign = m_storage.cell_face_signs[static_cast<std::size_t>(position)];
            const std::size_t f = static_cast<std::size_t>(face);
            closure += static_cast<double>(sign) * m_storage.face_area_vectors[f];
            total_area += m_storage.face_areas[f];
            const Index point_first = m_storage.face_point_offsets[f];
            const Index point_last = m_storage.face_point_offsets[f + 1U];
            std::vector<Index> ring;
            ring.reserve(static_cast<std::size_t>(point_last - point_first));
            for (Index point = point_first; point < point_last; ++point) {
                ring.push_back(m_storage.face_point_ids[static_cast<std::size_t>(point)]);
            }
            if (sign < 0) std::reverse(ring.begin(), ring.end());
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
    for (Index patch = 0; patch < patchCount(); ++patch) {
        const BoundaryPatch& boundary = m_storage.patches[static_cast<std::size_t>(patch)];
        if (boundary.name.empty()) throw std::runtime_error("mesh has an unnamed patch");
        for (Index face : boundary.faces) {
            if (face < 0 || face >= faces ||
                m_storage.face_patch[static_cast<std::size_t>(face)] != patch ||
                patch_seen[static_cast<std::size_t>(face)]) {
                throw std::runtime_error("polyhedral patch connectivity is invalid");
            }
            patch_seen[static_cast<std::size_t>(face)] = true;
        }
    }
    for (Index face = 0; face < faces; ++face) {
        if (m_storage.face_neighbour[static_cast<std::size_t>(face)] == invalid_index &&
            !patch_seen[static_cast<std::size_t>(face)]) {
            throw std::runtime_error("polyhedral boundary face is not assigned to a patch");
        }
    }
    Index owned_count = 0;
    for (Index cell = 0; cell < cells; ++cell) {
        const Index owned = m_storage.cell_owned_indices[static_cast<std::size_t>(cell)];
        if (owned != invalid_index) {
            if (owned != owned_count++ ||
                m_storage.cell_ghost_depths[static_cast<std::size_t>(cell)] != 0) {
                throw std::runtime_error("owned cell mapping is invalid");
            }
        }
    }
    if (owned_count != ownedCellCount() || owned_count == 0) {
        throw std::runtime_error("owned cell count is invalid");
    }
}
}  // namespace babelsim
