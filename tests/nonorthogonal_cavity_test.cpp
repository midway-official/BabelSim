#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "physics/simple/algorithm.h"
#include "babelsim/runtime.h"

#include "test_util.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

using namespace babelsim;

namespace {

Mesh warpedCavity(Index n) {
    constexpr double pi = 3.14159265358979323846;
    const std::array<Index, 3> dimensions{{n, n, 1}};
    std::vector<Vec3> points;
    points.reserve(static_cast<std::size_t>((n + 1) * (n + 1) * 2));
    for (Index k = 0; k <= 1; ++k) {
        for (Index j = 0; j <= n; ++j) {
            for (Index i = 0; i <= n; ++i) {
                const double xi = static_cast<double>(i) / n;
                const double eta = static_cast<double>(j) / n;
                const double envelope = std::sin(pi * xi) * std::sin(pi * eta);
                points.push_back({
                    xi + 0.12 * envelope,
                    eta + 0.04 * std::sin(2.0 * pi * xi) * std::sin(pi * eta),
                    static_cast<double>(k),
                });
            }
        }
    }
    auto patches = defaultPatches();
    for (Side side : {Side::XMin, Side::XMax, Side::YMin, Side::YMax}) {
        patches[static_cast<std::size_t>(side)].kind = PatchKind::Wall;
    }
    patches[static_cast<std::size_t>(Side::ZMin)].kind = PatchKind::Symmetry;
    patches[static_cast<std::size_t>(Side::ZMax)].kind = PatchKind::Symmetry;
    return Mesh::structured(dimensions, std::move(points), patches);
}

}  // 匿名命名空间

int main() {
    constexpr Index n = 10;
    const Mesh mesh = warpedCavity(n);
    double maximum_nonorthogonal_ratio = 0.0;
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        maximum_nonorthogonal_ratio = std::max(
            maximum_nonorthogonal_ratio,
            norm(detail::meshData(mesh).face_non_orthogonal[f]) / detail::meshData(mesh).face_areas[f]);
    }
    require(
        maximum_nonorthogonal_ratio > 0.05,
        "warped cavity is not sufficiently non-orthogonal");

    IncompressibleFields fields(mesh);
    for (Side side : {Side::XMin, Side::XMax, Side::YMin}) {
        fields.velocity.setBoundary(
            static_cast<Index>(side),
            BoundaryCondition<Vec3>::fixedValue({}));
    }
    fields.velocity.setBoundary(
        static_cast<Index>(Side::YMax),
        BoundaryCondition<Vec3>::fixedValue({1.0, 0.0, 0.0}));
    for (Side side : {Side::ZMin, Side::ZMax}) {
        fields.velocity.setBoundary(
            static_cast<Index>(side), BoundaryCondition<Vec3>::symmetry());
        fields.pressure.setBoundary(
            static_cast<Index>(side), BoundaryCondition<double>::symmetry());
    }

    RuntimeControl run_control;
    run_control.scalar_solver.solver = LinearSolverType::ConjugateGradient;
    run_control.scalar_solver.preconditioner = PreconditionerType::IncompleteCholesky;
    Methods& methods = run_control.methods;
    methods.gradient = GradientMethod::LeastSquares;
    methods.convection = ConvectionMethod::Upwind;
    methods.diffusion = DiffusionMethod::Corrected;
    SimpleControl control;
    control.max_iterations = 3000;
    control.velocity_relaxation = 0.4;
    control.pressure_relaxation = 0.2;
    control.continuity_tolerance = 1e-8;
    control.velocity_tolerance = 1e-6;
    run_control.vector_solver.absolute_tolerance = 1e-15;
    run_control.vector_solver.relative_tolerance = 1e-9;
    run_control.scalar_solver.absolute_tolerance = 1e-15;
    run_control.scalar_solver.relative_tolerance = 1e-9;

    RunTime run_time = RunTime::forMesh(mesh, run_control);
    SteadySimpleAlgorithm solver(fields, {1.0, 0.01}, control);
    SimpleIterationResult result;
    int iterations = 0;
    for (int iteration = 1; iteration <= control.max_iterations; ++iteration) {
        result = solver.iterate();
        iterations = iteration;
        require(result.healthy, "non-orthogonal cavity SIMPLE became unhealthy");
        if (result.converged) {
            break;
        }
    }
    require(result.converged, "non-orthogonal cavity SIMPLE did not converge");
    const Index centre = mesh.cellId(n / 2 - 1, n / 2 - 1, 0);
    require(
        detail::fieldData(fields.velocity)[centre].x < -0.04,
        "non-orthogonal cavity primary vortex is missing");
    require(
        result.continuity.relative <= control.continuity_tolerance,
        "non-orthogonal cavity is not mass conservative");

    std::cout << "nonorthogonal_cavity_test: cells=" << mesh.cellCount()
              << " iterations=" << iterations
              << " mass=" << result.continuity.relative
              << " centreU=" << detail::fieldData(fields.velocity)[centre].x
              << " maxKbySf=" << maximum_nonorthogonal_ratio << '\n';
}
