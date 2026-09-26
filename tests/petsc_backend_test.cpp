#include "internal/compute_backend.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"
#include "internal/petsc_session.h"

#include "babelsim/parallel.h"
#include "babelsim/mpi_support.h"

#include "test_util.h"

#include <mpi.h>

#include <cmath>
#include <iostream>
#include <vector>

using namespace babelsim;

namespace {

double exactValue(const Mesh& mesh, Index cell) {
    return 1.0 + 0.01 * static_cast<double>(detail::globalCellId(mesh, cell));
}

void assembleManufactured(ScalarDiscreteEquation& equation, double scale = 1.0) {
    const Mesh& mesh = *equation.mesh;
    std::fill(equation.diagonal.begin(), equation.diagonal.end(), 6.0);
    std::fill(equation.upper.begin(), equation.upper.end(), -1.0);
    std::fill(equation.lower.begin(), equation.lower.end(), -1.0);
    std::fill(equation.source.begin(), equation.source.end(), 0.0);
    for (Index cell : detail::meshData(mesh).owned_cells)
        equation.source[static_cast<std::size_t>(cell)] = 6.0 * exactValue(mesh, cell);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const Index owner = mesh.owner(face);
        const Index neighbour = mesh.neighbour(face);
        if (neighbour == invalid_index) continue;
        if (detail::isOwned(mesh, owner))
            equation.source[static_cast<std::size_t>(owner)] -= exactValue(mesh, neighbour);
        if (detail::isOwned(mesh, neighbour))
            equation.source[static_cast<std::size_t>(neighbour)] -= exactValue(mesh, owner);
    }
    for (Index cell : detail::meshData(mesh).owned_cells)
        equation.source[static_cast<std::size_t>(cell)] *= scale;
}

void verifyField(const ScalarField& field, double scale, double tolerance) {
    for (Index cell = 0; cell < field.mesh().cellCount(); ++cell)
        require(std::abs(detail::fieldValues(field)[static_cast<std::size_t>(cell)] -
                         scale * exactValue(field.mesh(), cell)) <= tolerance,
                "PETSc field differs from manufactured solution");
}

}  // namespace

int main(int argc, char** argv) {
    int initialized = 0;
    detail::checkMpi(MPI_Initialized(&initialized), "MPI_Initialized");
    const bool owns_mpi = initialized == 0;
    if (owns_mpi) detail::checkMpi(MPI_Init(&argc, &argv), "MPI_Init");
    try {
        const ParallelContext parallel = ParallelContext::world();
        const Mesh global = makeHexBox({4, 2, 2}, {0, 0, 0}, {1, 1, 1});
        const Mesh mesh = parallel.size == 1 ? global : decompose(global, parallel, 3);
        auto backend = detail::makeComputeBackend(mesh, parallel);

        // Only one owner changes its values after all ranks cache a valid halo.
        // Every rank must still participate in the next synchronization.
        ScalarField asymmetric(mesh, FieldLocation::Cell, "asymmetric");
        asymmetric.fill(0.0);
        backend->synchronize(asymmetric);
        if (parallel.rank == 0) {
            double* values = detail::fieldData(asymmetric);
            for (Index cell : detail::meshData(mesh).owned_cells) values[cell] = 7.0;
        }
        backend->synchronize(asymmetric);
        for (Index cell = 0; cell < mesh.cellCount(); ++cell)
            require(detail::fieldValues(asymmetric)[cell] ==
                        (detail::cellOwnerRank(mesh, cell) == 0 ? 7.0 : 0.0),
                    "rank-local update did not reach cached ghost values");
        const auto cached = backend->performance().halo_exchanges;
        backend->synchronize(asymmetric);
        require(backend->performance().halo_exchanges == cached,
                "globally valid halo was exchanged redundantly");

        LinearSolverConfig config;
        config.solver = LinearSolverType::ConjugateGradient;
        config.preconditioner = PreconditionerType::BlockJacobi;
        config.absolute_tolerance = 1e-13;
        config.relative_tolerance = 1e-12;
        config.max_iterations = 200;

        ScalarField scalar(mesh, FieldLocation::Cell, "petscScalar");
        ScalarDiscreteEquation scalar_equation(mesh);
        scalar_equation.numerical_identity = "petscScalarEquation";
        assembleManufactured(scalar_equation);
        SolveResult scalar_first = backend->solve(scalar_equation, scalar, config);
        require(scalar_first.converged(), "PETSc scalar manufactured solve did not converge");
        backend->synchronize(scalar);
        verifyField(scalar, 1.0, 2e-11);

        assembleManufactured(scalar_equation, 2.0);
        SolveResult scalar_rhs_only = backend->solve(scalar_equation, scalar, config);
        require(scalar_rhs_only.converged(), "PETSc scalar RHS-only solve did not converge");
        backend->synchronize(scalar);
        verifyField(scalar, 2.0, 4e-11);

        VectorField vector(mesh, FieldLocation::Cell, "petscVector");
        VectorDiscreteEquation vector_equation(mesh);
        vector_equation.numerical_identity = "petscVectorEquation";
        std::fill(vector_equation.diagonal.begin(), vector_equation.diagonal.end(), 6.0);
        std::fill(vector_equation.upper.begin(), vector_equation.upper.end(), -1.0);
        std::fill(vector_equation.lower.begin(), vector_equation.lower.end(), -1.0);
        for (Index cell : detail::meshData(mesh).owned_cells) {
            const double base = exactValue(mesh, cell);
            vector_equation.source[static_cast<std::size_t>(cell)] = {6.0 * base, -12.0 * base, 3.0 * base};
        }
        for (Index face = 0; face < mesh.faceCount(); ++face) {
            const Index owner = mesh.owner(face);
            const Index neighbour = mesh.neighbour(face);
            if (neighbour == invalid_index) continue;
            for (int endpoint = 0; endpoint < 2; ++endpoint) {
                const Index row = endpoint == 0 ? owner : neighbour;
                const Index column = endpoint == 0 ? neighbour : owner;
                if (!detail::isOwned(mesh, row)) continue;
                const double base = exactValue(mesh, column);
                Vec3& source = vector_equation.source[static_cast<std::size_t>(row)];
                source.x -= base;
                source.y += 2.0 * base;
                source.z -= 0.5 * base;
            }
        }
        const auto vector_results = backend->solve(vector_equation, vector, config);
        for (const SolveResult& result : vector_results) {
            if (!result.converged())
                std::cerr << "vector solve status=" << static_cast<int>(result.status)
                          << " iterations=" << result.iterations
                          << " initial=" << result.initial_residual
                          << " final=" << result.final_residual
                          << " relative=" << result.relative_residual << '\n';
            require(result.converged(), "PETSc vector component did not converge");
        }
        backend->synchronize(vector);
        for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
            const double base = exactValue(mesh, cell);
            const Vec3& actual = detail::fieldValues(vector)[static_cast<std::size_t>(cell)];
            require(std::abs(actual.x - base) < 2e-11, "PETSc x component is inaccurate");
            require(std::abs(actual.y + 2.0 * base) < 4e-11, "PETSc y component is inaccurate");
            require(std::abs(actual.z - 0.5 * base) < 1e-11, "PETSc z component is inaccurate");
        }

        std::unique_ptr<detail::ComputeBackend> split_backend;
        if (parallel.size == 1) {
            // Two triangular subfaces connect the same cell pair. The COO
            // pattern therefore contains repeated coordinates whose values
            // must be summed by PETSc without losing either face contribution.
            const Mesh split_mesh = makeSplitInterfaceMesh();
            split_backend = detail::makeComputeBackend(split_mesh, parallel);
            ScalarField split_field(split_mesh, FieldLocation::Cell, "splitInterface");
            ScalarDiscreteEquation split_equation(split_mesh);
            split_equation.numerical_identity = "splitInterfaceEquation";
            std::fill(split_equation.diagonal.begin(), split_equation.diagonal.end(), 2.0);
            std::fill(split_equation.upper.begin(), split_equation.upper.end(), -0.5);
            std::fill(split_equation.lower.begin(), split_equation.lower.end(), -0.5);
            for (Index cell : detail::meshData(split_mesh).owned_cells)
                split_equation.source[static_cast<std::size_t>(cell)] = 2.0 * exactValue(split_mesh, cell);
            for (Index face = 0; face < split_mesh.faceCount(); ++face) {
                const Index owner = split_mesh.owner(face);
                const Index neighbour = split_mesh.neighbour(face);
                if (neighbour == invalid_index) continue;
                split_equation.source[static_cast<std::size_t>(owner)] -= 0.5 * exactValue(split_mesh, neighbour);
                split_equation.source[static_cast<std::size_t>(neighbour)] -= 0.5 * exactValue(split_mesh, owner);
            }
            require(split_backend->solve(split_equation, split_field, config).converged(),
                    "PETSc repeated-COO interface solve did not converge");
            for (Index cell = 0; cell < split_mesh.cellCount(); ++cell)
                require(std::abs(detail::fieldValues(split_field)[static_cast<std::size_t>(cell)] -
                                 exactValue(split_mesh, cell)) < 1e-12,
                        "PETSc repeated-COO interface lost a face contribution");
        }

        if (parallel.distributed()) {
            // One disconnected cube per rank: a valid partition with zero
            // communication neighbours still participates in collective solves.
            Mesh isolated = makeHexBox({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
            std::vector<Index> face_ids(static_cast<std::size_t>(isolated.faceCount()));
            for (Index face = 0; face < isolated.faceCount(); ++face)
                face_ids[face] = parallel.rank * isolated.faceCount() + face;
            detail::MeshAccess::setPartition(isolated, parallel.size, 3,
                {parallel.rank}, {parallel.rank}, {0}, std::move(face_ids),
                std::vector<Index>(static_cast<std::size_t>(isolated.faceCount()), parallel.rank),
                parallel.rank);
            auto isolated_backend = detail::makeComputeBackend(isolated, parallel);
            ScalarField isolated_field(isolated, FieldLocation::Cell, "isolated");
            ScalarDiscreteEquation isolated_equation(isolated);
            isolated_equation.diagonal[0] = 2.0;
            isolated_equation.source[0] = 2.0 * (parallel.rank + 1.0);
            require(isolated_backend->solve(isolated_equation, isolated_field, config).converged(),
                    "zero-neighbour partition failed to solve");
            require(std::abs(detail::fieldValues(isolated_field)[0] - (parallel.rank + 1.0)) < 1e-12,
                    "zero-neighbour partition solution is incorrect");
            require(isolated_backend->performance().halo_bytes == 0,
                    "zero-neighbour partition transferred halo values");
        }

        const PerformanceCounters counters = backend->performance();
        require(counters.matrix_pattern_builds == 2, "PETSc did not retain one system per equation");
        require(counters.matrix_value_updates == 2, "PETSc matrix update count is incorrect");
        require(counters.rhs_only_solves == 1, "PETSc did not detect the unchanged scalar matrix");
        require(counters.true_residual_checks == 5, "PETSc true residual checks were not recorded");
        require(counters.equation_systems.size() == 2,
                "PETSc per-equation performance snapshots are missing");
        if (parallel.rank == 0)
            std::cout << "petsc_backend_test: ranks=" << parallel.size
                      << " owned=" << detail::ownedCellCount(mesh)
                      << " iterations=" << counters.krylov_iterations << '\n';
        backend.reset();
        split_backend.reset();
        detail::finalizePetscSession();
    } catch (const std::exception& error) {
        std::cerr << "petsc_backend_test: " << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    if (owns_mpi) detail::checkMpi(MPI_Finalize(), "MPI_Finalize");
}
