#include "babelsim/simple.h"
#include "babelsim/simple_control.h"

#include "test_util.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace babelsim;

namespace {

Mesh warpedCavity(Index n) {
    constexpr double pi = 3.14159265358979323846;
    std::vector<Vec3> points;
    points.reserve(static_cast<std::size_t>((n + 1) * (n + 1) * (n + 1)));
    for (Index k = 0; k <= n; ++k) {
        for (Index j = 0; j <= n; ++j) {
            for (Index i = 0; i <= n; ++i) {
                const double xi = static_cast<double>(i) / n;
                const double eta = static_cast<double>(j) / n;
                const double zeta = static_cast<double>(k) / n;
                const double envelope =
                    std::sin(pi * xi) * std::sin(pi * eta) * std::sin(pi * zeta);
                points.push_back({
                    xi + 0.10 * envelope,
                    eta + 0.045 * std::sin(2.0 * pi * xi) *
                        std::sin(pi * eta) * std::sin(pi * zeta),
                    zeta + 0.035 * std::sin(pi * xi) *
                        std::sin(pi * eta) * std::sin(2.0 * pi * zeta),
                });
            }
        }
    }
    auto patches = defaultPatches();
    for (auto& patch : patches) {
        patch.kind = PatchKind::Wall;
    }
    return Mesh::structured({n, n, n}, std::move(points), patches);
}

}  // 匿名命名空间

int main() {
    constexpr Index n = 5;
    const Mesh mesh = warpedCavity(n);
    double maximum_nonorthogonal_ratio = 0.0;
    double maximum_skewness_ratio = 0.0;
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        maximum_nonorthogonal_ratio = std::max(
            maximum_nonorthogonal_ratio,
            norm(mesh.face_non_orthogonal[f]) / mesh.face_areas[f]);
        maximum_skewness_ratio = std::max(
            maximum_skewness_ratio,
            norm(mesh.face_skewness[f]) /
                std::cbrt(mesh.cell_volumes[
                    static_cast<std::size_t>(mesh.face_owner[f])]));
    }
    require(
        maximum_nonorthogonal_ratio > 0.05 && maximum_skewness_ratio > 1e-3,
        "3D warped cavity does not exercise non-orthogonality and skewness");

    IncompressibleFields fields(mesh);
    for (Index patch = 0;
         patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        fields.velocity.setBoundary(
            patch, BoundaryCondition<Vec3>::fixedValue({}));
    }
    fields.velocity.setBoundary(
        static_cast<Index>(Side::YMax),
        BoundaryCondition<Vec3>::fixedValue({1.0, 0.0, 0.0}));

    Methods methods;
    methods.interpolation = InterpolationMethod::Corrected;
    methods.gradient = GradientMethod::LeastSquares;
    methods.convection = ConvectionMethod::Upwind;
    methods.diffusion = DiffusionMethod::Corrected;
    SimpleControl control;
    control.max_iterations = 3000;
    control.non_orthogonal_corrections = 2;
    control.velocity_relaxation = 0.4;
    control.pressure_relaxation = 0.2;
    control.continuity_tolerance = 2e-8;
    control.velocity_tolerance = 2e-6;
    control.velocity_solver.absolute_tolerance = 1e-15;
    control.velocity_solver.relative_tolerance = 1e-9;
    control.pressure_solver.absolute_tolerance = 1e-15;
    control.pressure_solver.relative_tolerance = 1e-9;

    RunTime run_time = RunTime::forMesh(mesh, simpleRunTimeControl(methods, control));
    SimpleSolver solver(run_time, fields, {1.0, 0.01}, control);
    SimpleIterationResult result;
    int iterations = 0;
    for (int iteration = 1; iteration <= control.max_iterations; ++iteration) {
        result = solver.iterate();
        iterations = iteration;
        require(result.healthy, "3D non-orthogonal cavity became unhealthy");
        if (result.converged) break;
    }
    require(result.converged, "3D non-orthogonal cavity did not converge");

    const Index centre = mesh.cellId(n / 2, n / 2, n / 2);
    require(
        fields.velocity[centre].x < -0.015,
        "3D non-orthogonal cavity primary vortex is missing");
    double maximum_spanwise_velocity = 0.0;
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        maximum_spanwise_velocity = std::max(
            maximum_spanwise_velocity, std::abs(fields.velocity[cell].z));
    }
    require(
        maximum_spanwise_velocity > 1e-4,
        "3D non-orthogonal cavity collapsed to a duplicated 2D solution");
    require(
        result.continuity.relative <= control.continuity_tolerance,
        "3D non-orthogonal cavity is not mass conservative");

    std::cout << "nonorthogonal_cavity_3d_test: cells=" << mesh.cellCount()
              << " iterations=" << iterations
              << " mass=" << result.continuity.relative
              << " centreU=" << fields.velocity[centre].x
              << " maxW=" << maximum_spanwise_velocity
              << " maxKbySf=" << maximum_nonorthogonal_ratio
              << " maxSkew=" << maximum_skewness_ratio << '\n';
}
