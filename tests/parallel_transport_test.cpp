#include "babelsim/parallel.h"
#include "babelsim/runtime.h"

#include "test_util.h"

#include <mpi.h>

#include <iostream>

using namespace babelsim;

int main(int argc, char* argv[]) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 1;
    const ParallelContext parallel = ParallelContext::world();
    try {
        const Mesh global = Mesh::cartesian({8, 2, 1}, {0, 0, 0}, {1, 1, 1});
        const Mesh mesh = decompose(global, parallel);
        ScalarField concentration(mesh, FieldLocation::Cell, "C", 0.0);
        ScalarField flux(mesh, FieldLocation::Face, "phi", 0.0);
        for (Index patch = 0; patch < static_cast<Index>(mesh.patches.size()); ++patch) {
            concentration.boundary(patch) = zeroGradient();
        }
        RuntimeControl control;
        control.methods.time = TimeMethod::Euler;
        control.time = {0.0, 0.1, 0.1};
        control.scalar_solver.absolute_tolerance = 1e-14;
        RunTime run_time = RunTime::forMesh(mesh, control);
        require(run_time.loop(), "parallel transport did not start its time step");
        const SolveResult result = solve(
            fvm::ddt(2.0, concentration) + fvm::div(flux, concentration) ==
                fvm::laplacian(0.0, concentration) + fvm::source(6.0));
        require(result.converged(), "parallel transport did not converge");
        require(!run_time.loop() && run_time.step() == 1, "parallel transport step count changed");
        double local_error = 0.0;
        for (Index cell : mesh.owned_cells) {
            local_error = std::max(local_error, std::abs(concentration[cell] - 0.3));
        }
        double global_error = 0.0;
        parallel.maximum(&local_error, &global_error, 1);
        require(global_error < 1e-12, "parallel transport changed the solution");
        if (parallel.rank == 0) {
            std::cout << "parallel_transport_test: ranks=" << parallel.size
                      << " max_error=" << global_error << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "rank " << parallel.rank << ": " << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
}
