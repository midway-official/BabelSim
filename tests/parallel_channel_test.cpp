#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/simple.h"
#include "babelsim/simple_control.h"
#include "babelsim/runtime.h"
#include "babelsim/mesh_io.h"
#include "babelsim/parallel.h"
#include "babelsim/parallel_writer.h"

#include "test_util.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace babelsim;

namespace {

void configureChannelBoundaries(IncompressibleFields& fields) {
    const Mesh& mesh = fields.velocity.mesh();
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        switch (detail::meshData(mesh).patches[static_cast<std::size_t>(patch)].kind) {
            case PatchKind::Inlet:
                fields.velocity.setBoundary(
                    patch, BoundaryCondition<Vec3>::fixedValue({1.0, 0.0, 0.0}));
                fields.pressure.setBoundary(
                    patch, BoundaryCondition<double>::zeroGradient());
                break;
            case PatchKind::Outlet:
                fields.velocity.setBoundary(
                    patch, BoundaryCondition<Vec3>::zeroGradient());
                fields.pressure.setBoundary(
                    patch, BoundaryCondition<double>::fixedValue(0.0));
                break;
            case PatchKind::Wall:
                fields.velocity.setBoundary(
                    patch, BoundaryCondition<Vec3>::fixedValue({}));
                fields.pressure.setBoundary(
                    patch, BoundaryCondition<double>::zeroGradient());
                break;
            case PatchKind::Symmetry:
                fields.velocity.setBoundary(
                    patch, BoundaryCondition<Vec3>::symmetry());
                fields.pressure.setBoundary(
                    patch, BoundaryCondition<double>::symmetry());
                break;
            case PatchKind::Generic:
            case PatchKind::Processor:
                fields.velocity.setBoundary(
                    patch, BoundaryCondition<Vec3>::zeroGradient());
                fields.pressure.setBoundary(
                    patch, BoundaryCondition<double>::zeroGradient());
                break;
        }
    }
}

}  // 匿名命名空间

int main(int argc, char* argv[]) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 1;
    const ParallelContext parallel = ParallelContext::world();
    try {
        if (argc < 2 || argc > 3) {
            throw std::invalid_argument(
                "usage: parallel_channel_test <BabelSim-mesh> [output-dir]");
        }
        require(
            parallel.size == 1 || parallel.size == 2 || parallel.size == 4,
            "parallel_channel_test supports one, two, or four MPI ranks");
        // 验证原生网格的根 rank 读取与局部分发路径；测试进程不保留全局 Mesh。
        const Mesh mesh = readDistributedMesh(argv[1], parallel);
        IncompressibleFields fields(mesh);
        configureChannelBoundaries(fields);

        RuntimeControl run_control;
        run_control.scalar_solver.solver = LinearSolverType::ConjugateGradient;
        run_control.scalar_solver.preconditioner = PreconditionerType::IncompleteCholesky;
        Methods& methods = run_control.methods;
        methods.gradient = GradientMethod::GreenGauss;
        methods.convection = ConvectionMethod::Upwind;
        methods.diffusion = DiffusionMethod::Orthogonal;
        SimpleControl control;
        control.max_iterations = 1200;
        control.velocity_relaxation = 0.5;
        control.pressure_relaxation = 0.3;
        control.continuity_tolerance = 1e-7;
        control.velocity_tolerance = 1e-6;
        run_control.vector_solver.absolute_tolerance = 1e-14;
        run_control.vector_solver.relative_tolerance = 1e-7;
        run_control.vector_solver.max_iterations = 300;
        run_control.scalar_solver.absolute_tolerance = 1e-14;
        run_control.scalar_solver.relative_tolerance = 1e-7;
        run_control.scalar_solver.max_iterations = 1200;

        RunTime run_time = RunTime::forMesh(mesh, run_control);
        SimpleSolver solver(fields, {1.0, 0.01}, control);
        SimpleIterationResult result;
        int iterations = 0;
        for (int iteration = 1; iteration <= control.max_iterations; ++iteration) {
            result = solver.iterate();
            iterations = iteration;
            require(result.healthy, "distributed channel SIMPLE became unhealthy");
            if (result.converged) {
                break;
            }
        }
        require(result.converged, "distributed channel SIMPLE did not converge");

        Vec3 local_sum{};
        double local_max = 0.0;
        for (Index cell : detail::meshData(mesh).owned_cells) {
            local_sum += detail::fieldData(fields.velocity)[cell];
            local_max = std::max(local_max, norm(detail::fieldData(fields.velocity)[cell]));
        }
        const double local_sum_values[3] = {
            local_sum.x, local_sum.y, local_sum.z};
        double global_sum_values[3] = {};
        parallel.sum(local_sum_values, global_sum_values, 3);
        double global_max = 0.0;
        parallel.maximum(&local_max, &global_max, 1);
        require(global_max > 1.0, "native Poiseuille flow is empty");
        if (parallel.rank == 0) {
            std::cout << "parallel_channel_test: ranks=" << parallel.size
                      << " iterations=" << iterations
                      << " mass=" << result.continuity.relative
                      << " dU=" << result.relative_velocity_change
                      << " sumU=(" << global_sum_values[0] << ','
                      << global_sum_values[1] << ',' << global_sum_values[2]
                      << ") maxU=" << global_max << '\n';
        }
        const std::filesystem::path output = argc == 3
            ? argv[2] : "build/mpi-output";
        const auto time = output / ("poiseuille" + std::to_string(parallel.size));
        writeOwnedFieldCsv(time, fields.velocity, parallel);
        writeOwnedFieldCsv(time, fields.pressure, parallel);
        writeOwnedResultMetadata(time, mesh, parallel, "poiseuille", {
            {"U", "vector", FieldLocation::Cell}, {"p", "scalar", FieldLocation::Cell}});
    } catch (const std::exception& error) {
        std::cerr << "rank " << parallel.rank << ": " << error.what() << '\n';
        const int abort_status = MPI_Abort(parallel.communicator, 1);
        if (abort_status != MPI_SUCCESS) return 1;
    }
    return MPI_Finalize() == MPI_SUCCESS ? 0 : 1;
}
