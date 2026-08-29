#include "babelsim/incompressible.h"
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
        require(parallel.size == 2, "parallel_cavity_3d_test requires two MPI ranks");
        constexpr Index n = 6;
        auto patches = defaultPatches();
        for (auto& patch : patches) {
            patch.kind = PatchKind::Wall;
        }
        const Mesh global = Mesh::cartesian(
            {n, n, n}, {0, 0, 0}, {1, 1, 1}, patches);
        const Mesh mesh = decompose(global, parallel);
        IncompressibleFields fields(mesh);
        for (Index patch = 0;
             patch < static_cast<Index>(global.patches.size()); ++patch) {
            fields.velocity.setBoundary(
                patch, BoundaryCondition<Vec3>::fixedValue({}));
        }
        fields.velocity.setBoundary(
            static_cast<Index>(Side::YMax),
            BoundaryCondition<Vec3>::fixedValue({1.0, 0.0, 0.0}));

        Methods methods;
        methods.gradient = GradientMethod::GreenGauss;
        methods.convection = ConvectionMethod::Upwind;
        methods.diffusion = DiffusionMethod::Orthogonal;
        SimpleControl control;
        control.max_iterations = 2000;
        control.velocity_relaxation = 0.5;
        control.pressure_relaxation = 0.3;
        control.continuity_tolerance = 1e-8;
        control.velocity_tolerance = 1e-6;
        control.velocity_solver.absolute_tolerance = 1e-15;
        control.velocity_solver.relative_tolerance = 1e-9;
        control.pressure_solver.absolute_tolerance = 1e-15;
        control.pressure_solver.relative_tolerance = 1e-9;

        SimpleSolver solver(fields, {1.0, 0.01}, methods, control, parallel);
        SimpleIterationResult result;
        int iterations = 0;
        for (int iteration = 1; iteration <= control.max_iterations; ++iteration) {
            result = solver.iterate();
            iterations = iteration;
            require(result.healthy, "distributed 3D cavity became unhealthy");
            if (result.converged) {
                break;
            }
        }
        require(result.converged, "distributed 3D cavity did not converge");

        const Index lower_global = global.cellId(n / 2 - 1, n / 2 - 1, n / 2 - 1);
        const Index upper_global = global.cellId(n / 2 - 1, n / 2 - 1, n / 2);
        Vec3 lower{};
        Vec3 upper{};
        double local_maximum_z = 0.0;
        for (Index cell : mesh.owned_cells) {
            const Index id = mesh.globalCellId(cell);
            if (id == lower_global) {
                lower = fields.velocity[cell];
            }
            if (id == upper_global) {
                upper = fields.velocity[cell];
            }
            local_maximum_z = std::max(
                local_maximum_z, std::abs(fields.velocity[cell].z));
        }
        const double local_centres[6] = {
            lower.x, lower.y, lower.z, upper.x, upper.y, upper.z};
        double centres[6]{};
        parallel.sum(local_centres, centres, 6);
        double maximum_z = 0.0;
        parallel.maximum(&local_maximum_z, &maximum_z, 1);
        require(
            0.5 * (centres[0] + centres[3]) < -0.02,
            "distributed 3D cavity primary vortex is missing");
        require(
            near(centres[0], centres[3], 2e-5) &&
                near(centres[1], centres[4], 2e-5) &&
                near(centres[2], -centres[5], 2e-5),
            "distributed 3D cavity violates centre-plane symmetry");
        require(
            maximum_z > 1e-4,
            "distributed 3D cavity collapsed to a duplicated 2D solution");
        require(
            result.continuity.relative <= control.continuity_tolerance,
            "distributed 3D cavity is not mass conservative");

        const std::filesystem::path output = argc > 1
            ? argv[1] : "build/mpi-output";
        const auto time = output / "cavity3d";
        writeOwnedFieldCsv(time, fields.velocity, parallel);
        writeOwnedFieldCsv(time, fields.pressure, parallel);
        writeOwnedResultMetadata(time, mesh, parallel, "cavity3d", {
            {"U", "vector", FieldLocation::Cell}, {"p", "scalar", FieldLocation::Cell}});
        if (parallel.rank == 0) {
            std::cout << "parallel_cavity_3d_test: ranks=" << parallel.size
                      << " iterations=" << iterations
                      << " mass=" << result.continuity.relative
                      << " centreU=" << 0.5 * (centres[0] + centres[3])
                      << " maxW=" << maximum_z << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "rank " << parallel.rank << ": " << error.what() << '\n';
        const int abort_status = MPI_Abort(parallel.communicator, 1);
        if (abort_status != MPI_SUCCESS) return 1;
    }
    return MPI_Finalize() == MPI_SUCCESS ? 0 : 1;
}
