#include "babelsim/incompressible.h"
#include "babelsim/legacy_taiho.h"
#include "babelsim/parallel.h"

#include "test_util.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

using namespace babelsim;

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    const ParallelContext parallel = ParallelContext::world();
    try {
        if (argc < 2 || argc > 3) {
            throw std::invalid_argument(
                "usage: parallel_imported_test <Taiho-mesh> [output-dir]");
        }
        require(
            parallel.size == 1 || parallel.size == 2 || parallel.size == 4,
            "parallel_imported_test supports one, two, or four MPI ranks");
        ImportedTaihoMesh imported = readTaihoMesh(argv[1]);
        IncompressibleFields global_fields(imported.mesh);
        applyImportedBoundaryConditions(
            imported, global_fields.velocity, global_fields.pressure);
        const Mesh mesh = decompose(imported.mesh, parallel);
        IncompressibleFields fields(mesh);
        copyBoundaryConditions(global_fields.velocity, fields.velocity);
        copyBoundaryConditions(global_fields.pressure, fields.pressure);

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
            require(result.healthy, "distributed imported SIMPLE became unhealthy");
            if (result.converged) {
                break;
            }
        }
        require(result.converged, "distributed imported SIMPLE did not converge");

        Vec3 local_sum{};
        double local_max = 0.0;
        for (Index cell : mesh.owned_cells) {
            local_sum += fields.velocity[cell];
            local_max = std::max(local_max, norm(fields.velocity[cell]));
        }
        double global_sum_values[3] = {};
        const double local_sum_values[3] = {
            local_sum.x, local_sum.y, local_sum.z};
        parallel.sum(local_sum_values, global_sum_values, 3);
        double global_max = 0.0;
        parallel.maximum(&local_max, &global_max, 1);
        require(global_max > 1.0, "imported Poiseuille flow is empty");
        if (parallel.rank == 0) {
            std::cout << "parallel_imported_test: ranks=" << parallel.size
                      << " iterations=" << iterations
                      << " mass=" << result.continuity.relative
                      << " dU=" << result.relative_velocity_change
                      << " sumU=(" << global_sum_values[0] << ','
                      << global_sum_values[1] << ',' << global_sum_values[2]
                      << ") maxU=" << global_max << '\n';
        }
        const std::filesystem::path output = argc == 3
            ? argv[2] : "build/mpi-output";
        writeOwnedCsv(
            output, fields.velocity, fields.pressure, parallel,
            "poiseuille" + std::to_string(parallel.size));
    } catch (const std::exception& error) {
        std::cerr << "rank " << parallel.rank << ": " << error.what() << '\n';
        MPI_Abort(parallel.communicator, 1);
    }
    MPI_Finalize();
    return 0;
}
