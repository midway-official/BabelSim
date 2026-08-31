#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/simple.h"
#include "babelsim/simple_control.h"
#include "babelsim/runtime.h"
#include "babelsim/parallel.h"
#include "babelsim/parallel_writer.h"

#include "test_util.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

using namespace babelsim;

int main(int argc, char* argv[]) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 1;
    const ParallelContext parallel = ParallelContext::world();
    try {
        require(
            parallel.size == 1 || parallel.size == 2 || parallel.size == 4,
            "parallel_simple_test supports one, two, or four MPI ranks");
        constexpr Index n = 12;
        auto patches = defaultPatches();
        for (Side side : {Side::XMin, Side::XMax, Side::YMin, Side::YMax}) {
            patches[static_cast<std::size_t>(side)].kind = PatchKind::Wall;
        }
        patches[static_cast<std::size_t>(Side::ZMin)].kind = PatchKind::Symmetry;
        patches[static_cast<std::size_t>(Side::ZMax)].kind = PatchKind::Symmetry;
        const Mesh global = Mesh::cartesian(
            {n, n, 1}, {0, 0, 0}, {1, 1, 1}, patches);
        const Mesh mesh = decompose(global, parallel);
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
                static_cast<Index>(side),
                BoundaryCondition<Vec3>::symmetry());
            fields.pressure.setBoundary(
                static_cast<Index>(side),
                BoundaryCondition<double>::symmetry());
        }

        RuntimeControl run_control;
        run_control.scalar_solver.solver = LinearSolverType::ConjugateGradient;
        run_control.scalar_solver.preconditioner = PreconditionerType::IncompleteCholesky;
        Methods& methods = run_control.methods;
        methods.gradient = GradientMethod::GreenGauss;
        methods.convection = ConvectionMethod::Upwind;
        methods.diffusion = DiffusionMethod::Orthogonal;
        SimpleControl control;
        control.max_iterations = 3000;
        control.velocity_relaxation = 0.5;
        control.pressure_relaxation = 0.3;
        control.continuity_tolerance = 1e-8;
        control.velocity_tolerance = 1e-6;
        run_control.vector_solver.absolute_tolerance = 1e-14;
        run_control.vector_solver.relative_tolerance = 1e-10;
        run_control.scalar_solver.absolute_tolerance = 1e-14;
        run_control.scalar_solver.relative_tolerance = 1e-10;

        RunTime run_time = RunTime::forMesh(mesh, run_control);
        SimpleSolver solver(run_time, fields, {1.0, 0.01}, control);
        SimpleIterationResult result;
        int iterations = 0;
        for (int iteration = 1;
             iteration <= control.max_iterations; ++iteration) {
            result = solver.iterate();
            iterations = iteration;
            require(result.healthy, "distributed cavity SIMPLE became unhealthy");
            if (result.converged) {
                break;
            }
        }
        require(result.converged, "distributed cavity SIMPLE did not converge");
        // 收敛决定必须是 collective；若某个 rank 使用了本地状态提前/延后
        // 跳出，这里的最小和最大外迭代次数会立即暴露控制流分叉。
        const int maximum_iterations = parallel.maximum(iterations);
        const int minimum_iterations = -parallel.maximum(-iterations);
        require(
            minimum_iterations == maximum_iterations,
            "MPI ranks stopped SIMPLE at different outer iterations");

        const Index centre_global = global.cellId(n / 2 - 1, n / 2 - 1, 0);
        double local_centre_u = 0.0;
        double local_maximum_z = 0.0;
        for (Index cell : detail::meshData(mesh).owned_cells) {
            if (detail::globalCellId(mesh, cell) == centre_global) {
                local_centre_u = detail::fieldData(fields.velocity)[cell].x;
            }
            local_maximum_z = std::max(
                local_maximum_z, std::abs(detail::fieldData(fields.velocity)[cell].z));
        }
        double centre_u = 0.0;
        double maximum_z = 0.0;
        parallel.sum(&local_centre_u, &centre_u, 1);
        parallel.maximum(&local_maximum_z, &maximum_z, 1);
        require(
            std::abs(centre_u - (-0.149236)) < 5e-4,
            "distributed cavity differs excessively from the serial reference");
        require(maximum_z < 1e-12, "distributed nz=1 cavity generated z velocity");

        const std::filesystem::path output = argc > 1
            ? argv[1] : "build/mpi-output";
        const auto time = output / ("cavity" + std::to_string(parallel.size));
        writeOwnedFieldCsv(time, fields.velocity, parallel);
        writeOwnedFieldCsv(time, fields.pressure, parallel);
        writeOwnedResultMetadata(time, mesh, parallel, "cavity", {
            {"U", "vector", FieldLocation::Cell}, {"p", "scalar", FieldLocation::Cell}});
        if (parallel.rank == 0) {
            std::cout << "parallel_simple_test: ranks=" << parallel.size
                      << " iterations=" << iterations
                      << " mass=" << result.continuity.relative
                      << " dU=" << result.relative_velocity_change
                      << " centreU=" << centre_u << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "rank " << parallel.rank << ": " << error.what() << '\n';
        const int abort_status = MPI_Abort(parallel.communicator, 1);
        if (abort_status != MPI_SUCCESS) return 1;
    }
    return MPI_Finalize() == MPI_SUCCESS ? 0 : 1;
}
