#pragma once

#include "babelsim/mesh.h"
#include "babelsim/vector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

inline void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

inline bool near(double actual, double expected, double tolerance = 1e-11) {
    return std::abs(actual - expected) <= tolerance *
        std::max({1.0, std::abs(actual), std::abs(expected)});
}

inline bool near(
    const babelsim::Vec3& actual,
    const babelsim::Vec3& expected,
    double tolerance = 1e-11)
{
    return near(actual.x, expected.x, tolerance) &&
        near(actual.y, expected.y, tolerance) &&
        near(actual.z, expected.z, tolerance);
}

// Test-only input producer. It emits explicit Hex cells and boundary quads before
// calling the sole public Mesh construction API; production code has no regular-grid API.
inline std::array<babelsim::PatchSpec, 6> boxPatches() {
    using namespace babelsim;
    return {{{"minus_x", PatchKind::Generic}, {"plus_x", PatchKind::Generic},
             {"minus_y", PatchKind::Generic}, {"plus_y", PatchKind::Generic},
             {"minus_z", PatchKind::Generic}, {"plus_z", PatchKind::Generic}}};
}

inline babelsim::Index hexCellIndex(
    babelsim::Index a, babelsim::Index b, babelsim::Index c,
    babelsim::Index count_a, babelsim::Index count_b)
{
    return a + count_a * (b + count_b * c);
}

inline babelsim::Mesh makeHexFromVertices(
    std::array<babelsim::Index, 3> counts,
    std::vector<babelsim::Vec3> vertices,
    const std::array<babelsim::PatchSpec, 6>& patches = boxPatches())
{
    using namespace babelsim;
    if (counts[0] <= 0 || counts[1] <= 0 || counts[2] <= 0) {
        throw std::invalid_argument("test hex counts must be positive");
    }
    const auto vertex = [&](Index a, Index b, Index c) {
        return a + (counts[0] + 1) * (b + (counts[1] + 1) * c);
    };
    const std::size_t expected = static_cast<std::size_t>(counts[0] + 1) *
        static_cast<std::size_t>(counts[1] + 1) * static_cast<std::size_t>(counts[2] + 1);
    if (vertices.size() != expected) throw std::invalid_argument("test vertex count is invalid");
    std::vector<std::array<Index, 8>> cells;
    std::vector<BoundaryFaceSpec> boundaries;
    cells.reserve(static_cast<std::size_t>(counts[0]) * counts[1] * counts[2]);
    for (Index c = 0; c < counts[2]; ++c) {
        for (Index b = 0; b < counts[1]; ++b) {
            for (Index a = 0; a < counts[0]; ++a) {
                const std::array<Index, 8> cell{{
                    vertex(a, b, c), vertex(a + 1, b, c), vertex(a + 1, b + 1, c), vertex(a, b + 1, c),
                    vertex(a, b, c + 1), vertex(a + 1, b, c + 1),
                    vertex(a + 1, b + 1, c + 1), vertex(a, b + 1, c + 1)}};
                cells.push_back(cell);
                if (a == 0) boundaries.push_back({{{cell[0], cell[4], cell[7], cell[3]}}, 0});
                if (a + 1 == counts[0]) boundaries.push_back({{{cell[1], cell[2], cell[6], cell[5]}}, 1});
                if (b == 0) boundaries.push_back({{{cell[0], cell[1], cell[5], cell[4]}}, 2});
                if (b + 1 == counts[1]) boundaries.push_back({{{cell[3], cell[7], cell[6], cell[2]}}, 3});
                if (c == 0) boundaries.push_back({{{cell[0], cell[3], cell[2], cell[1]}}, 4});
                if (c + 1 == counts[2]) boundaries.push_back({{{cell[4], cell[5], cell[6], cell[7]}}, 5});
            }
        }
    }
    return Mesh::unstructured(
        std::move(vertices), std::move(cells),
        std::vector<PatchSpec>(patches.begin(), patches.end()), std::move(boundaries));
}

inline babelsim::Mesh makeHexBox(
    std::array<babelsim::Index, 3> counts,
    babelsim::Vec3 minimum,
    babelsim::Vec3 maximum,
    const std::array<babelsim::PatchSpec, 6>& patches = boxPatches())
{
    using namespace babelsim;
    const Vec3 span = maximum - minimum;
    if (counts[0] <= 0 || counts[1] <= 0 || counts[2] <= 0 ||
        span.x <= 0.0 || span.y <= 0.0 || span.z <= 0.0) {
        throw std::invalid_argument("test hex box is invalid");
    }
    std::vector<Vec3> vertices;
    vertices.reserve(static_cast<std::size_t>(counts[0] + 1) *
                     static_cast<std::size_t>(counts[1] + 1) *
                     static_cast<std::size_t>(counts[2] + 1));
    for (Index c = 0; c <= counts[2]; ++c) {
        for (Index b = 0; b <= counts[1]; ++b) {
            for (Index a = 0; a <= counts[0]; ++a) {
                vertices.push_back({minimum.x + span.x * a / counts[0],
                                    minimum.y + span.y * b / counts[1],
                                    minimum.z + span.z * c / counts[2]});
            }
        }
    }
    return makeHexFromVertices(counts, std::move(vertices), patches);
}
