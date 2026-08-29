#include "babelsim/incompressible.h"
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

        Methods methods;
        methods.gradient = GradientMethod::GreenGauss;
        methods.convection = ConvectionMethod::Upwind;
        methods.diffusion = DiffusionMethod::Orthogonal;
        SimpleControl control;
        control.max_iterations = 3000;
        control.velocity_relaxation = 0.5;
        control.pressure_relaxation = 0.3;
        control.continuity_tolerance = 1e-8;
        control.velocity_tolerance = 1e-6;
        control.velocity_solver.absolute_tolerance = 1e-14;
        control.velocity_solver.relative_tolerance = 1e-10;
        control.pressure_solver.absolute_tolerance = 1e-14;
        control.pressure_solver.relative_tolerance = 1e-10;

        SimpleSolver solver(
            fields, {1.0, 0.01}, methods, control, parallel);
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

        const Index centre_global = global.cellId(n / 2 - 1, n / 2 - 1, 0);
        double local_centre_u = 0.0;
        double local_maximum_z = 0.0;
        for (Index cell : mesh.owned_cells) {
            if (mesh.globalCellId(cell) == centre_global) {
                local_centre_u = fields.velocity[cell].x;
            }
            local_maximum_z = std::max(
                local_maximum_z, std::abs(fields.velocity[cell].z));
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
        writeOwnedCsv(
            output, fields.velocity, fields.pressure, parallel,
            "cavity" + std::to_string(parallel.size));
        if (parallel.rank == 0) {
            std::cout << "parallel_simple_test: ranks=" << parallel.size
                      << " iterations=" << iterations
                      << " mass=" << result.continuity.relative
                      << " dU=" << result.relative_velocity_change
                      << " centreU=" << centre_u << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "rank " << parallel.rank << ": " << error.what() << '\n';
        MPI_Abort(parallel.communicator, 1);
    }
    MPI_Finalize();
    return 0;
}
