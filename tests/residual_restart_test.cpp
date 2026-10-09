#include "babelsim/parallel.h"
#include "internal/compute_backend.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"
#include "internal/petsc_session.h"
#include "test_util.h"
#include <mpi.h>
#include <cmath>
#include <iostream>

using namespace babelsim;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    try {
        const auto parallel = ParallelContext::world();
        const auto global = makeHexBox({64, 4, 1}, {0, 0, 0}, {1, 1, 1});
        const auto mesh = decompose(global, parallel, 3);
        auto backend = detail::makeComputeBackend(mesh, parallel);
        ScalarDiscreteEquation equation(mesh);
        equation.numerical_identity = "residualContract";
        const auto exact = [&](Index cell) {
            return 1.25 + std::sin(0.13 * detail::globalCellId(mesh, cell));
        };
        std::fill(equation.diagonal.begin(), equation.diagonal.end(), 0.02);
        for (Index face = 0; face < mesh.faceCount(); ++face) {
            const Index o = mesh.owner(face), n = mesh.neighbour(face);
            if (n == invalid_index) continue;
            equation.diagonal[o] += 1.0; equation.diagonal[n] += 1.0;
            equation.upper[face] = equation.lower[face] = -1.0;
        }
        for (Index cell : detail::meshData(mesh).owned_cells)
            equation.source[cell] = equation.diagonal[cell] * exact(cell);
        for (Index face = 0; face < mesh.faceCount(); ++face) {
            const Index o = mesh.owner(face), n = mesh.neighbour(face);
            if (n == invalid_index) continue;
            if (detail::isOwned(mesh, o)) equation.source[o] -= exact(n);
            if (detail::isOwned(mesh, n)) equation.source[n] -= exact(o);
        }
        LinearSolverConfig config;
        config.solver = LinearSolverType::ConjugateGradient;
        config.preconditioner = PreconditionerType::None;
        config.relative_tolerance = 1e-30;
        config.max_iterations = 1000;
        int repaired = 0, rejected = 0;
        for (double target : {1e-12, 1e-13, 6e-14, 4e-14, 2e-14, 1e-14}) {
            config.absolute_tolerance = target;
            ScalarField x(mesh, FieldLocation::Cell, "manufactured");
            const auto before = backend->performance();
            const auto result = backend->solve(equation, x, config);
            const auto after = backend->performance();
            // Continuation keeps the original acceptance target and shares
            // the initial iteration budget. It must not silently rescale it.
            require(result.iterations <= config.max_iterations, "residual continuation renewed the iteration budget");
            require(result.converged() == residualConverged(result.final_residual, target),
                    "true residual status disagrees with the original absolute target");
            if (target >= 1e-13) require(result.converged(), "attainable manufactured target failed");
            if (result.converged()) {
                if (after.true_residual_checks > before.true_residual_checks + 1) ++repaired;
                double local_error = 0.0, error = 0.0;
                for (Index cell : detail::meshData(mesh).owned_cells)
                    local_error = std::max(local_error, std::abs(detail::fieldData(x)[cell] - exact(cell)));
                parallel.maximum(&local_error, &error, 1);
                require(error < 1e-9, "converged residual lost the independent manufactured solution");
            } else {
                ++rejected;
                require(std::isfinite(result.final_residual), "unattainable residual became nonfinite");
            }
        }
        // A good warm guess skips Krylov without changing the true residual.
        config.absolute_tolerance = 1e-12;
        config.warm_start = true;
        ScalarField exact_field(mesh, FieldLocation::Cell, "exactWarm");
        for (Index cell : detail::meshData(mesh).owned_cells) detail::fieldData(exact_field)[cell] = exact(cell);
        const auto before = backend->performance();
        const auto warm = backend->solve(equation, exact_field, config);
        const auto after = backend->performance();
        require(warm.converged() && warm.iterations == 0 && warm.final_residual == warm.initial_residual,
                "accepted initial true residual changed on the short path");
        require(after.linear_solve_seconds == before.linear_solve_seconds,
                "accepted warm start performed a redundant Krylov solve");
        if (parallel.rank == 0)
            std::cout << "residual_restart_test: ranks=" << parallel.size << " repaired=" << repaired
                      << " honestRejections=" << rejected << " strict targets and warm shortcut passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; MPI_Abort(MPI_COMM_WORLD, 1);
    }
    detail::finalizePetscSession();
    MPI_Finalize();
}
