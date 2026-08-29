#include "babelsim/incompressible.h"
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
    for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
        switch (mesh.patches[static_cast<std::size_t>(patch)].kind) {
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
    MPI_Init(&argc, &argv);
    const ParallelContext parallel = ParallelContext::world();
    try {
        if (argc < 2 || argc > 3) {
            throw std::invalid_argument(
                "usage: parallel_channel_test <BabelSim-mesh> [output-dir]");
        }
        require(
            parallel.size == 1 || parallel.size == 2 || parallel.size == 4,
            "parallel_channel_test supports one, two, or four MPI ranks");
        const Mesh global = readMeshFile(argv[1]);
        const Mesh mesh = decompose(global, parallel);
        IncompressibleFields fields(mesh);
        configureChannelBoundaries(fields);

        Methods methods;
        methods.gradient = GradientMethod::GreenGauss;
        methods.convection = ConvectionMethod::Upwind;
        methods.diffusion = DiffusionMethod::Orthogonal;
        SimpleControl control;
        control.max_iterations = 1200;
        control.velocity_relaxation = 0.5;
        control.pressure_relaxation = 0.3;
        control.continuity_tolerance = 1e-7;
        control.velocity_tolerance = 1e-6;
        control.velocity_solver.absolute_tolerance = 1e-14;
        control.velocity_solver.relative_tolerance = 1e-7;
        control.velocity_solver.max_iterations = 300;
        control.pressure_solver.absolute_tolerance = 1e-14;
        control.pressure_solver.relative_tolerance = 1e-7;
        control.pressure_solver.max_iterations = 1200;

        SimpleSolver solver(
            fields, {1.0, 0.01}, methods, control, parallel);
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
        for (Index cell : mesh.owned_cells) {
            local_sum += fields.velocity[cell];
            local_max = std::max(local_max, norm(fields.velocity[cell]));
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
        MPI_Abort(parallel.communicator, 1);
    }
    MPI_Finalize();
    return 0;
}
