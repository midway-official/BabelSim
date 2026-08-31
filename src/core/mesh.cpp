#include "babelsim/mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <utility>

namespace babelsim {
namespace {

std::size_t checkedProduct(Index a, Index b, Index c, const char* label) {
    if (a <= 0 || b <= 0 || c <= 0) {
        throw std::invalid_argument(std::string(label) + " dimensions must be positive");
    }
    const std::int64_t product =
        static_cast<std::int64_t>(a) * b * c;
    if (product > std::numeric_limits<Index>::max()) {
        throw std::overflow_error(std::string(label) + " exceeds 32-bit indexing");
    }
    return static_cast<std::size_t>(product);
}

Index checkedIncrement(Index value, const char* label) {
    if (value <= 0 || value == std::numeric_limits<Index>::max()) {
        throw std::overflow_error(std::string(label) + " dimension exceeds 32-bit indexing");
    }
    return static_cast<Index>(value + 1);
}

template <typename VertexStorage>
std::pair<Vec3, Vec3> quadGeometry(
    const std::array<Index, 4>& ids,
    const VertexStorage& vertices)
{
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
        throw std::runtime_error("structured mesh contains a degenerate face");
    }
    const Vec3 centre =
        (magnitude_abc * ((a + b + c) / 3.0) +
         magnitude_acd * ((a + c + d) / 3.0)) /
        magnitude;
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

}  // 匿名命名空间

std::array<PatchSpec, 6> defaultPatches() {
    return {{
        {"xmin", PatchKind::Generic},
        {"xmax", PatchKind::Generic},
        {"ymin", PatchKind::Generic},
        {"ymax", PatchKind::Generic},
        {"zmin", PatchKind::Generic},
        {"zmax", PatchKind::Generic},
    }};
}

Index Mesh::cellCount() const {
    return static_cast<Index>(checkedProduct(
        m_storage.dimensions[0], m_storage.dimensions[1], m_storage.dimensions[2], "cell"));
}

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

void Mesh::setOwnership(
    std::array<Index, 3> global_cells,
    Index global_x_offset,
    Index owned_x_begin,
    Index owned_x_end,
    Index halo_layers)
{
    const auto [nx, ny, nz] = m_storage.dimensions;
    if (global_cells[0] < nx || global_cells[1] != ny ||
        global_cells[2] != nz || global_x_offset < 0 ||
        static_cast<std::int64_t>(global_x_offset) + nx > global_cells[0] ||
        owned_x_begin < 0 || owned_x_begin > owned_x_end ||
        owned_x_end > nx || halo_layers < 0) {
        throw std::invalid_argument("structured mesh ownership is invalid");
    }
    m_storage.global_dimensions = global_cells;
    m_storage.global_i_offset = global_x_offset;
    m_storage.owned_i_begin = owned_x_begin;
    m_storage.owned_i_end = owned_x_end;
    m_storage.ghost_layers = halo_layers;

    const auto cells = static_cast<std::size_t>(cellCount());
    m_storage.owned_cells.clear();
    m_storage.owned_cells.reserve(static_cast<std::size_t>(owned_x_end - owned_x_begin) *
                        static_cast<std::size_t>(ny) *
                        static_cast<std::size_t>(nz));
    m_storage.cell_owned_indices.assign(cells, invalid_index);
    m_storage.cell_global_ids.resize(cells);
    for (Index k = 0; k < nz; ++k) {
        for (Index j = 0; j < ny; ++j) {
            for (Index i = 0; i < nx; ++i) {
                const Index local = cellId(i, j, k);
                m_storage.cell_global_ids[static_cast<std::size_t>(local)] =
                    global_x_offset + i + global_cells[0] * (j + ny * k);
                if (i >= owned_x_begin && i < owned_x_end) {
                    m_storage.cell_owned_indices[static_cast<std::size_t>(local)] =
                        static_cast<Index>(m_storage.owned_cells.size());
                    m_storage.owned_cells.push_back(local);
                }
            }
        }
    }
    m_storage.owned_faces.clear();
    m_storage.owned_faces.reserve(static_cast<std::size_t>(faceCount()));
    for (Index face = 0; face < faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = m_storage.face_owner[f];
        const Index neighbour = m_storage.face_neighbour[f];
        if (isOwned(owner) ||
            (neighbour != invalid_index && isOwned(neighbour))) {
            m_storage.owned_faces.push_back(face);
        }
    }
}

void Mesh::setPatches(std::vector<BoundaryPatch> new_patches) {
    if (new_patches.empty()) {
        throw std::invalid_argument("mesh must contain at least one boundary patch");
    }
    m_storage.patches.clear();
    m_storage.patches.reserve(new_patches.size());
    for (BoundaryPatch& patch : new_patches) {
        m_storage.patches.push_back(std::move(patch));
    }
}

void Mesh::addPatchFace(Index patch, Index face) {
    if (patch < 0 || static_cast<std::size_t>(patch) >= m_storage.patches.size() ||
        face < 0 || static_cast<std::size_t>(face) >=
            static_cast<std::size_t>(faceCount()) ||
        m_storage.face_neighbour[static_cast<std::size_t>(face)] != invalid_index) {
        throw std::out_of_range("mesh patch face is invalid");
    }
    m_storage.patches[static_cast<std::size_t>(patch)].faces.push_back(face);
}

Index Mesh::cellId(Index i, Index j, Index k) const {
    const auto [nx, ny, nz] = m_storage.dimensions;
    if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) {
        throw std::out_of_range("cell index is outside the structured mesh");
    }
    return i + nx * (j + ny * k);
}

Index Mesh::vertexId(Index i, Index j, Index k) const {
    const auto [nx, ny, nz] = m_storage.dimensions;
    if (i < 0 || i > nx || j < 0 || j > ny || k < 0 || k > nz) {
        throw std::out_of_range("vertex index is outside the structured mesh");
    }
    return i + (nx + 1) * (j + (ny + 1) * k);
}

Mesh Mesh::cartesian(
    std::array<Index, 3> cells,
    Vec3 minimum,
    Vec3 maximum,
    const std::array<PatchSpec, 6>& boundary)
{
    checkedProduct(cells[0], cells[1], cells[2], "cell");
    const Vec3 extent = maximum - minimum;
    if (!(extent.x > 0.0) || !(extent.y > 0.0) || !(extent.z > 0.0) ||
        !isFinite(minimum) || !isFinite(maximum)) {
        throw std::invalid_argument("Cartesian mesh bounds must be finite and increasing");
    }

    std::vector<Vec3> points;
    points.reserve(checkedProduct(
        checkedIncrement(cells[0], "vertex"),
        checkedIncrement(cells[1], "vertex"),
        checkedIncrement(cells[2], "vertex"), "vertex"));
    for (Index k = 0; k <= cells[2]; ++k) {
        for (Index j = 0; j <= cells[1]; ++j) {
            for (Index i = 0; i <= cells[0]; ++i) {
                points.push_back({
                    minimum.x + extent.x * i / cells[0],
                    minimum.y + extent.y * j / cells[1],
                    minimum.z + extent.z * k / cells[2],
                });
            }
        }
    }
    return structured(cells, std::move(points), boundary);
}

Mesh Mesh::structured(
    std::array<Index, 3> cells,
    std::vector<Vec3> points,
    const std::array<PatchSpec, 6>& boundary)
{
    const std::size_t number_of_cells =
        checkedProduct(cells[0], cells[1], cells[2], "cell");
    const std::size_t number_of_vertices = checkedProduct(
        checkedIncrement(cells[0], "vertex"),
        checkedIncrement(cells[1], "vertex"),
        checkedIncrement(cells[2], "vertex"), "vertex");
    if (points.size() != number_of_vertices) {
        throw std::invalid_argument("structured vertex array has the wrong size");
    }
    if (!std::all_of(points.begin(), points.end(), isFinite)) {
        throw std::invalid_argument("structured vertex array contains non-finite values");
    }

    Mesh mesh;
    mesh.m_storage.dimensions = cells;
    mesh.m_storage.vertices.assign(std::move(points));
    mesh.m_storage.cell_centres.resize(number_of_cells);
    mesh.m_storage.cell_volumes.resize(number_of_cells);
    mesh.m_storage.cell_inverse_volumes.resize(number_of_cells);
    mesh.m_storage.cell_faces.resize(number_of_cells);
    mesh.m_storage.cell_neighbours.resize(number_of_cells);
    for (auto& faces : mesh.m_storage.cell_faces) {
        faces.fill(invalid_index);
    }
    for (auto& neighbours : mesh.m_storage.cell_neighbours) {
        neighbours.fill(invalid_index);
    }
    mesh.m_storage.patches.reserve(boundary.size());
    for (const PatchSpec& patch : boundary) {
        mesh.m_storage.patches.push_back({patch.name, patch.kind, {}});
    }

    // 顶点平均值足以确定面的朝向；精确多面体质心在下方通过四面体分解计算。
    for (Index k = 0; k < cells[2]; ++k) {
        for (Index j = 0; j < cells[1]; ++j) {
            for (Index i = 0; i < cells[0]; ++i) {
                Vec3 centre{};
                for (Index dk = 0; dk <= 1; ++dk) {
                    for (Index dj = 0; dj <= 1; ++dj) {
                        for (Index di = 0; di <= 1; ++di) {
                            centre += mesh.m_storage.vertices[static_cast<std::size_t>(
                                mesh.vertexId(i + di, j + dj, k + dk))];
                        }
                    }
                }
                mesh.m_storage.cell_centres[static_cast<std::size_t>(mesh.cellId(i, j, k))] =
                    centre / 8.0;
            }
        }
    }

    const auto addFace = [&](std::array<Index, 4> vertices,
                             Index owner,
                             Index neighbour,
                             Index patch,
                             Side owner_side,
                             Side neighbour_side) {
        const Index face = static_cast<Index>(mesh.m_storage.face_owner.size());
        auto [centre, area_vector] = quadGeometry(vertices, mesh.m_storage.vertices);
        const Vec3 direction = neighbour == invalid_index
            ? centre - mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)]
            : mesh.m_storage.cell_centres[static_cast<std::size_t>(neighbour)] -
                mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)];
        if (dot(area_vector, direction) < 0.0) {
            std::swap(vertices[1], vertices[3]);
            area_vector = -area_vector;
        }

        mesh.m_storage.face_vertices.push_back(vertices);
        mesh.m_storage.face_owner.push_back(owner);
        mesh.m_storage.face_neighbour.push_back(neighbour);
        mesh.m_storage.face_patch.push_back(patch);
        mesh.m_storage.face_centres.push_back(centre);
        mesh.m_storage.face_area_vectors.push_back(area_vector);
        mesh.m_storage.cell_faces[static_cast<std::size_t>(owner)]
            [static_cast<std::size_t>(owner_side)] = face;
        if (neighbour != invalid_index) {
            mesh.m_storage.cell_faces[static_cast<std::size_t>(neighbour)]
                [static_cast<std::size_t>(neighbour_side)] = face;
            mesh.m_storage.cell_neighbours[static_cast<std::size_t>(owner)]
                [static_cast<std::size_t>(owner_side)] = neighbour;
            mesh.m_storage.cell_neighbours[static_cast<std::size_t>(neighbour)]
                [static_cast<std::size_t>(neighbour_side)] = owner;
        } else {
            mesh.m_storage.patches[static_cast<std::size_t>(patch)].faces.push_back(face);
        }
    };

    // x 法向面；顶点顺序对应正逻辑 x 法向。
    for (Index k = 0; k < cells[2]; ++k) {
        for (Index j = 0; j < cells[1]; ++j) {
            for (Index i = 0; i <= cells[0]; ++i) {
                const Index owner = i == 0
                    ? mesh.cellId(0, j, k)
                    : mesh.cellId(i - 1, j, k);
                const Index neighbour = i == 0 || i == cells[0]
                    ? invalid_index
                    : mesh.cellId(i, j, k);
                const Index patch = i == 0
                    ? static_cast<Index>(Side::XMin)
                    : (i == cells[0] ? static_cast<Index>(Side::XMax) : invalid_index);
                addFace({
                    mesh.vertexId(i, j, k),
                    mesh.vertexId(i, j + 1, k),
                    mesh.vertexId(i, j + 1, k + 1),
                    mesh.vertexId(i, j, k + 1),
                }, owner, neighbour, patch,
                i == 0 ? Side::XMin : Side::XMax, Side::XMin);
            }
        }
    }

    // y 法向面；顶点顺序对应正逻辑 y 法向。
    for (Index k = 0; k < cells[2]; ++k) {
        for (Index j = 0; j <= cells[1]; ++j) {
            for (Index i = 0; i < cells[0]; ++i) {
                const Index owner = j == 0
                    ? mesh.cellId(i, 0, k)
                    : mesh.cellId(i, j - 1, k);
                const Index neighbour = j == 0 || j == cells[1]
                    ? invalid_index
                    : mesh.cellId(i, j, k);
                const Index patch = j == 0
                    ? static_cast<Index>(Side::YMin)
                    : (j == cells[1] ? static_cast<Index>(Side::YMax) : invalid_index);
                addFace({
                    mesh.vertexId(i, j, k),
                    mesh.vertexId(i, j, k + 1),
                    mesh.vertexId(i + 1, j, k + 1),
                    mesh.vertexId(i + 1, j, k),
                }, owner, neighbour, patch,
                j == 0 ? Side::YMin : Side::YMax, Side::YMin);
            }
        }
    }

    // z 法向面；顶点顺序对应正逻辑 z 法向。
    for (Index k = 0; k <= cells[2]; ++k) {
        for (Index j = 0; j < cells[1]; ++j) {
            for (Index i = 0; i < cells[0]; ++i) {
                const Index owner = k == 0
                    ? mesh.cellId(i, j, 0)
                    : mesh.cellId(i, j, k - 1);
                const Index neighbour = k == 0 || k == cells[2]
                    ? invalid_index
                    : mesh.cellId(i, j, k);
                const Index patch = k == 0
                    ? static_cast<Index>(Side::ZMin)
                    : (k == cells[2] ? static_cast<Index>(Side::ZMax) : invalid_index);
                addFace({
                    mesh.vertexId(i, j, k),
                    mesh.vertexId(i + 1, j, k),
                    mesh.vertexId(i + 1, j + 1, k),
                    mesh.vertexId(i, j + 1, k),
                }, owner, neighbour, patch,
                k == 0 ? Side::ZMin : Side::ZMax, Side::ZMin);
            }
        }
    }

    // 对三角化六面体计算精确质心和体积。
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        const Vec3 reference = mesh.m_storage.cell_centres[static_cast<std::size_t>(cell)];
        double volume = 0.0;
        Vec3 first_moment{};
        for (Index face : mesh.m_storage.cell_faces[static_cast<std::size_t>(cell)]) {
            const auto& ids = mesh.m_storage.face_vertices[static_cast<std::size_t>(face)];
            const Vec3& a = mesh.m_storage.vertices[static_cast<std::size_t>(ids[0])];
            const Vec3& b = mesh.m_storage.vertices[static_cast<std::size_t>(ids[1])];
            const Vec3& c = mesh.m_storage.vertices[static_cast<std::size_t>(ids[2])];
            const Vec3& d = mesh.m_storage.vertices[static_cast<std::size_t>(ids[3])];
            addTetrahedron(reference, a, b, c, volume, first_moment);
            addTetrahedron(reference, a, c, d, volume, first_moment);
        }
        if (!(volume > 0.0) || !std::isfinite(volume)) {
            throw std::runtime_error("structured mesh contains a non-positive cell");
        }
        mesh.m_storage.cell_volumes[static_cast<std::size_t>(cell)] = volume;
        mesh.m_storage.cell_inverse_volumes[static_cast<std::size_t>(cell)] = 1.0 / volume;
        mesh.m_storage.cell_centres[static_cast<std::size_t>(cell)] = first_moment / volume;
    }

    const std::size_t number_of_faces = mesh.m_storage.face_owner.size();
    mesh.m_storage.face_areas.resize(number_of_faces);
    mesh.m_storage.face_normals.resize(number_of_faces);
    mesh.m_storage.face_orthogonal_coefficients.resize(number_of_faces);
    mesh.m_storage.face_owner_weights.resize(number_of_faces);
    mesh.m_storage.face_non_orthogonal.resize(number_of_faces);
    mesh.m_storage.face_skewness.resize(number_of_faces);
    mesh.m_storage.orthogonal_geometry = true;
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
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
            throw std::runtime_error(
                "face owner-neighbour orientation or geometry is invalid");
        }
        const double orthogonal = projected / delta_squared;
        mesh.m_storage.face_areas[f] = area;
        mesh.m_storage.face_normals[f] = area_vector / area;
        mesh.m_storage.face_orthogonal_coefficients[f] = orthogonal;
        mesh.m_storage.face_non_orthogonal[f] = area_vector - orthogonal * delta;

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
                throw std::runtime_error("face interpolation distance is invalid");
            }
            mesh.m_storage.face_owner_weights[f] = neighbour_distance / distance_sum;
            const double denominator = dot(delta, normal);
            const double fraction =
                dot(centre - mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)], normal) /
                denominator;
            const Vec3 intersection =
                mesh.m_storage.cell_centres[static_cast<std::size_t>(owner)] + fraction * delta;
            mesh.m_storage.face_skewness[f] = centre - intersection;
        }
        // 两个量的物理量纲不同：非正交向量与面面积同量纲，偏斜量与长度同量纲。
        // 分别归一化后再判断，避免极端长宽比网格被误判为正交网格。
        const double non_orthogonal_ratio =
            norm(mesh.m_storage.face_non_orthogonal[f]) / area;
        const double skew_scale = std::max(norm(delta), std::sqrt(area));
        const double skew_ratio =
            skew_scale > 0.0 ? norm(mesh.m_storage.face_skewness[f]) / skew_scale : 0.0;
        if (non_orthogonal_ratio > 1e-14 || skew_ratio > 1e-14) {
            mesh.m_storage.orthogonal_geometry = false;
        }
    }

    mesh.setOwnership(cells, 0, 0, cells[0], 0);
    mesh.validate();
    return mesh;
}

void Mesh::validate() const {
    const std::size_t cells = checkedProduct(
        m_storage.dimensions[0], m_storage.dimensions[1], m_storage.dimensions[2], "cell");
    const std::size_t expected_vertices = checkedProduct(
        checkedIncrement(m_storage.dimensions[0], "vertex"),
        checkedIncrement(m_storage.dimensions[1], "vertex"),
        checkedIncrement(m_storage.dimensions[2], "vertex"), "vertex");
    if (m_storage.vertices.size() != expected_vertices || m_storage.cell_centres.size() != cells ||
        m_storage.cell_volumes.size() != cells || m_storage.cell_inverse_volumes.size() != cells ||
        m_storage.cell_faces.size() != cells ||
        m_storage.cell_neighbours.size() != cells) {
        throw std::runtime_error("mesh cell or vertex arrays have inconsistent sizes");
    }
    if (m_storage.global_dimensions[0] <= 0 || m_storage.global_dimensions[1] <= 0 ||
        m_storage.global_dimensions[2] <= 0 || m_storage.owned_i_begin < 0 ||
        m_storage.owned_i_end < m_storage.owned_i_begin || m_storage.owned_i_end > m_storage.dimensions[0] ||
        m_storage.global_i_offset < 0 || m_storage.ghost_layers < 0) {
        throw std::runtime_error("mesh ownership arrays are inconsistent");
    }
    checkedProduct(
        m_storage.global_dimensions[0], m_storage.global_dimensions[1], m_storage.global_dimensions[2],
        "global cell");
    const std::int64_t expected_owned =
        static_cast<std::int64_t>(m_storage.owned_i_end - m_storage.owned_i_begin) *
        m_storage.dimensions[1] * m_storage.dimensions[2];
    const std::int64_t local_end =
        static_cast<std::int64_t>(m_storage.global_i_offset) + m_storage.dimensions[0];
    if (m_storage.global_dimensions[0] <= 0 || m_storage.global_dimensions[1] <= 0 ||
        m_storage.global_dimensions[2] <= 0 || m_storage.cell_owned_indices.size() != cells ||
        m_storage.cell_global_ids.size() != cells ||
        m_storage.owned_cells.size() != static_cast<std::size_t>(expected_owned) ||
        m_storage.owned_i_begin < 0 || m_storage.owned_i_end < m_storage.owned_i_begin ||
        m_storage.owned_i_end > m_storage.dimensions[0] || m_storage.global_i_offset < 0 ||
        local_end > m_storage.global_dimensions[0] ||
        m_storage.global_dimensions[1] != m_storage.dimensions[1] ||
        m_storage.global_dimensions[2] != m_storage.dimensions[2] || m_storage.ghost_layers < 0 ||
        (m_storage.owned_i_begin > 0 && m_storage.owned_i_begin < m_storage.ghost_layers) ||
        (m_storage.dimensions[0] - m_storage.owned_i_end > 0 &&
         m_storage.dimensions[0] - m_storage.owned_i_end < m_storage.ghost_layers)) {
        throw std::runtime_error("mesh ownership arrays are inconsistent");
    }
    for (Index k = 0; k < m_storage.dimensions[2]; ++k) {
        for (Index j = 0; j < m_storage.dimensions[1]; ++j) {
            for (Index i = 0; i < m_storage.dimensions[0]; ++i) {
                const Index cell = cellId(i, j, k);
                const auto c = static_cast<std::size_t>(cell);
                const std::int64_t expected_global =
                    static_cast<std::int64_t>(m_storage.global_i_offset) + i +
                    static_cast<std::int64_t>(m_storage.global_dimensions[0]) *
                        (j + static_cast<std::int64_t>(m_storage.global_dimensions[1]) * k);
                const bool owned = i >= m_storage.owned_i_begin && i < m_storage.owned_i_end;
                if (m_storage.cell_global_ids[c] != static_cast<Index>(expected_global) ||
                    (owned != (m_storage.cell_owned_indices[c] != invalid_index))) {
                    throw std::runtime_error(
                        "mesh cell ownership/global-id mapping is inconsistent");
                }
                if (!owned && m_storage.cell_owned_indices[c] != invalid_index) {
                    throw std::runtime_error("ghost cell has an owned index");
                }
            }
        }
    }
    for (std::size_t owned = 0; owned < m_storage.owned_cells.size(); ++owned) {
        const Index cell = m_storage.owned_cells[owned];
        if (cell < 0 || static_cast<std::size_t>(cell) >= cells ||
            m_storage.cell_owned_indices[static_cast<std::size_t>(cell)] !=
                static_cast<Index>(owned)) {
            throw std::runtime_error("mesh owned-cell mapping is inconsistent");
        }
    }
    const std::size_t faces = m_storage.face_owner.size();
    if (m_storage.face_vertices.size() != faces || m_storage.face_neighbour.size() != faces ||
        m_storage.face_patch.size() != faces || m_storage.face_centres.size() != faces ||
        m_storage.face_area_vectors.size() != faces || m_storage.face_normals.size() != faces ||
        m_storage.face_non_orthogonal.size() != faces ||
        m_storage.face_skewness.size() != faces || m_storage.face_areas.size() != faces ||
        m_storage.face_orthogonal_coefficients.size() != faces ||
        m_storage.face_owner_weights.size() != faces) {
        throw std::runtime_error("mesh face arrays have inconsistent sizes");
    }

    std::set<std::string> patch_names;
    std::vector<int> patch_visits(faces, 0);
    for (std::size_t patch = 0; patch < m_storage.patches.size(); ++patch) {
        if (m_storage.patches[patch].name.empty() ||
            static_cast<int>(m_storage.patches[patch].kind) < static_cast<int>(PatchKind::Generic) ||
            static_cast<int>(m_storage.patches[patch].kind) > static_cast<int>(PatchKind::Processor) ||
            !patch_names.insert(m_storage.patches[patch].name).second) {
            throw std::runtime_error("boundary patch names must be non-empty and unique");
        }
        for (Index face : m_storage.patches[patch].faces) {
            if (face < 0 || static_cast<std::size_t>(face) >= faces ||
                m_storage.face_patch[static_cast<std::size_t>(face)] !=
                    static_cast<Index>(patch) ||
                m_storage.face_neighbour[static_cast<std::size_t>(face)] != invalid_index) {
                throw std::runtime_error("boundary patch face mapping is inconsistent");
            }
            ++patch_visits[static_cast<std::size_t>(face)];
        }
    }

    for (std::size_t face = 0; face < faces; ++face) {
        const Index owner = m_storage.face_owner[face];
        const Index neighbour = m_storage.face_neighbour[face];
        const bool owner_valid = owner >= 0 && static_cast<std::size_t>(owner) < cells;
        const bool neighbour_valid =
            neighbour >= 0 && static_cast<std::size_t>(neighbour) < cells;
        if (!owner_valid || (!neighbour_valid && neighbour != invalid_index) ||
            !isFinite(m_storage.face_centres[face]) || !isFinite(m_storage.face_area_vectors[face]) ||
            !isFinite(m_storage.face_normals[face]) ||
            !isFinite(m_storage.face_non_orthogonal[face]) || !isFinite(m_storage.face_skewness[face]) ||
            !(m_storage.face_areas[face] > 0.0) ||
            !(m_storage.face_orthogonal_coefficients[face] > 0.0) ||
            !(m_storage.face_owner_weights[face] >= 0.0 && m_storage.face_owner_weights[face] <= 1.0)) {
            throw std::runtime_error("mesh contains an invalid face");
        }
        if (std::abs(norm(m_storage.face_normals[face]) - 1.0) > 1e-12) {
            throw std::runtime_error("mesh face normal is not normalized");
        }
        for (Index vertex : m_storage.face_vertices[face]) {
            if (vertex < 0 || static_cast<std::size_t>(vertex) >= m_storage.vertices.size()) {
                throw std::runtime_error("face vertex index is outside the mesh");
            }
        }
        if (neighbour == invalid_index) {
            if (m_storage.face_patch[face] < 0 ||
                static_cast<std::size_t>(m_storage.face_patch[face]) >= m_storage.patches.size() ||
                patch_visits[face] != 1) {
                throw std::runtime_error("boundary face has no unique patch");
            }
        } else if (m_storage.face_patch[face] != invalid_index || patch_visits[face] != 0) {
            throw std::runtime_error("internal face is assigned to a boundary patch");
        }
    }

    for (std::size_t cell = 0; cell < cells; ++cell) {
        if (!isFinite(m_storage.cell_centres[cell]) || !(m_storage.cell_volumes[cell] > 0.0) ||
            !std::isfinite(m_storage.cell_volumes[cell]) ||
            !(m_storage.cell_inverse_volumes[cell] > 0.0) ||
            !std::isfinite(m_storage.cell_inverse_volumes[cell])) {
            throw std::runtime_error("mesh contains an invalid cell geometry");
        }
        Vec3 closure{};
        double area_sum = 0.0;
        for (std::size_t side = 0; side < 6; ++side) {
            const Index face = m_storage.cell_faces[cell][side];
            if (face < 0 || static_cast<std::size_t>(face) >= faces) {
                throw std::runtime_error("cell is missing one of its six faces");
            }
            const std::size_t f = static_cast<std::size_t>(face);
            const bool owns = m_storage.face_owner[f] == static_cast<Index>(cell);
            const bool neighbours = m_storage.face_neighbour[f] == static_cast<Index>(cell);
            if (owns == neighbours) {
                throw std::runtime_error("cell-face owner/neighbour mapping is inconsistent");
            }
            const Index expected_neighbour = owns
                ? m_storage.face_neighbour[f] : m_storage.face_owner[f];
            if (m_storage.cell_neighbours[cell][side] != expected_neighbour) {
                throw std::runtime_error("cell-neighbour topology is inconsistent");
            }
            closure += owns ? m_storage.face_area_vectors[f] : -m_storage.face_area_vectors[f];
            area_sum += m_storage.face_areas[f];
        }
        if (norm(closure) > 1e-10 * std::max(area_sum, 1.0)) {
            throw std::runtime_error("cell face area vectors do not form a closed surface");
        }
    }

    std::size_t expected_owned_faces = 0;
    for (std::size_t face = 0; face < faces; ++face) {
        const Index owner = m_storage.face_owner[face];
        const Index neighbour = m_storage.face_neighbour[face];
        if (isOwned(owner) ||
            (neighbour != invalid_index && isOwned(neighbour))) {
            ++expected_owned_faces;
        }
    }
    if (m_storage.owned_faces.size() != expected_owned_faces) {
        throw std::runtime_error("mesh owned-face mapping is inconsistent");
    }
    for (Index face : m_storage.owned_faces) {
        if (face < 0 || static_cast<std::size_t>(face) >= faces) {
            throw std::runtime_error("mesh owned-face index is invalid");
        }
    }
}

}  // babelsim 命名空间
