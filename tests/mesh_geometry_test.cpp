#include "babelsim/mesh.h"

#include "test_util.h"

#include <algorithm>
#include <iostream>

using namespace babelsim;

int main() {
    const Mesh planar = Mesh::cartesian({2, 3, 1}, {0, 0, 0}, {2, 3, 1});
    require(planar.orthogonalGeometry(), "Cartesian mesh was not recognized as orthogonal");
    require(planar.cellCount() == 6, "nz=1 cell count is incorrect");
    require(planar.faceCount() == 29, "nz=1 face count is incorrect");
    require(planar.vertexCount() == 24, "nz=1 vertex count is incorrect");
    const std::array<std::size_t, 6> patch_sizes{{3, 3, 2, 2, 6, 6}};
    for (std::size_t patch = 0; patch < patch_sizes.size(); ++patch) {
        require(
            planar.patches[patch].faces.size() == patch_sizes[patch],
            "structured patch face count is incorrect");
    }

    Index internal_faces = 0;
    for (Index face = 0; face < planar.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        if (!planar.boundaryFace(face)) {
            ++internal_faces;
        }
        require(near(planar.face_areas[f], 1.0), "Cartesian face area is incorrect");
        require(
            norm(planar.face_non_orthogonal[f]) < 1e-12,
            "Cartesian face has a non-orthogonal correction");
        require(
            norm(planar.face_skewness[f]) < 1e-12,
            "Cartesian face has skewness");
    }
    require(internal_faces == 7, "structured internal-face count is incorrect");
    for (double volume : planar.cell_volumes) {
        require(near(volume, 1.0), "Cartesian cell volume is incorrect");
    }

    const std::array<Index, 3> dimensions{{3, 3, 2}};
    std::vector<Vec3> points;
    points.reserve(static_cast<std::size_t>(4 * 4 * 3));
    for (Index k = 0; k <= dimensions[2]; ++k) {
        for (Index j = 0; j <= dimensions[1]; ++j) {
            for (Index i = 0; i <= dimensions[0]; ++i) {
                points.push_back({
                    i + 0.22 * j + 0.025 * i * j + 0.015 * j * k,
                    j + 0.13 * k + 0.012 * i * k,
                    k + 0.018 * i * j,
                });
            }
        }
    }
    const Mesh skewed = Mesh::structured(dimensions, std::move(points));
    double maximum_non_orthogonal = 0.0;
    double maximum_skewness = 0.0;
    double total_volume = 0.0;
    for (Index face = 0; face < skewed.faceCount(); ++face) {
        maximum_non_orthogonal = std::max(
            maximum_non_orthogonal,
            norm(skewed.face_non_orthogonal[static_cast<std::size_t>(face)]));
        maximum_skewness = std::max(
            maximum_skewness,
            norm(skewed.face_skewness[static_cast<std::size_t>(face)]));
    }
    for (double volume : skewed.cell_volumes) {
        require(volume > 0.0, "skewed cell volume is not positive");
        total_volume += volume;
    }
    require(maximum_non_orthogonal > 1e-3, "skewed mesh was treated as orthogonal");
    require(maximum_skewness > 1e-5, "skewness metric was not detected");

    std::cout << "mesh_geometry_test: cells=" << skewed.cellCount()
              << " volume=" << total_volume
              << " max_nonorth=" << maximum_non_orthogonal
              << " max_skew=" << maximum_skewness << '\n';
}
