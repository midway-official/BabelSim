#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/parallel.h"
#include "babelsim/assembly.h"
#include "babelsim/distributed_solver.h"
#include "babelsim/operators.h"
#include "babelsim/parallel_writer.h"

#include "test_util.h"

#include <mpi.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace babelsim;

int main(int argc, char* argv[]) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 1;
    const ParallelContext parallel = ParallelContext::world();
    try {
        require(parallel.size == 2, "parallel_domain_test requires two MPI ranks");
        const Mesh global = Mesh::cartesian(
            {8, 3, 2}, {0, 0, 0}, {2, 1, 0.5});
        const Mesh local = decompose(global, parallel);
        require(
            parallel.sum(detail::ownedCellCount(local)) == global.cellCount(),
            "owned cell counts do not cover the global mesh");
        require(
            detail::ownedCellCount(local) == 24 && detail::meshData(local).ghost_layers == 2,
            "unexpected local ownership or halo width");

        HaloExchange halo(local, parallel);
        ScalarField scalar(local, FieldLocation::Cell, "scalar", -1.0);
        VectorField vector(local, FieldLocation::Cell, "vector", {-1, -1, -1});
        for (Index cell : detail::meshData(local).owned_cells) {
            const double id = detail::globalCellId(local, cell);
            detail::fieldData(scalar)[cell] = id + 0.25;
            detail::fieldData(vector)[cell] = {id, 2.0 * id, -id};
        }
        halo.exchange(scalar);
        halo.exchange(vector);
        for (Index cell = 0; cell < local.cellCount(); ++cell) {
            const double id = detail::globalCellId(local, cell);
            require(
                near(detail::fieldData(scalar)[cell], id + 0.25) &&
                near(detail::fieldData(vector)[cell], {id, 2.0 * id, -id}),
                "halo exchange did not reconstruct global cell values");
        }

        TensorField tensor(local, FieldLocation::Cell, "tensor");
        for (Index cell : detail::meshData(local).owned_cells) {
            const double id = detail::globalCellId(local, cell);
            detail::fieldData(tensor)[cell].rows[0] = {id, id + 1.0, id + 2.0};
            detail::fieldData(tensor)[cell].rows[1] = {id + 3.0, id + 4.0, id + 5.0};
            detail::fieldData(tensor)[cell].rows[2] = {id + 6.0, id + 7.0, id + 8.0};
        }
        halo.exchange(tensor);
        for (Index cell = 0; cell < local.cellCount(); ++cell) {
            const double id = detail::globalCellId(local, cell);
            for (std::size_t row = 0; row < 3; ++row) {
                for (std::size_t column = 0; column < 3; ++column) {
                    require(
                        near(detail::fieldData(tensor)[cell].rows[row][column],
                             id + 3.0 * row + column),
                        "tensor cell halo exchange did not reconstruct values");
                }
            }
        }

        ScalarField face_field(local, FieldLocation::Face, "faceField", 0.0);
        face_field.fill(static_cast<double>(parallel.rank));
        halo.exchange(face_field);
        VectorField face_vector(local, FieldLocation::Face, "faceVector");
        face_vector.fill({
            static_cast<double>(parallel.rank),
            static_cast<double>(parallel.rank + 10),
            static_cast<double>(parallel.rank + 20)});
        halo.exchange(face_vector);
        for (Index cell : detail::meshData(local).owned_cells) {
            const Index i = cell % detail::meshData(local).dimensions[0];
            if (i == detail::meshData(local).owned_i_begin && parallel.rank > 0) {
                const Index face = detail::meshData(local).cell_faces[static_cast<std::size_t>(cell)]
                    [static_cast<std::size_t>(Side::XMin)];
                require(
                    near(detail::fieldData(face_field)[face], 0.0),
                    "face-centred halo exchange failed on the left interface");
                require(
                    near(detail::fieldData(face_vector)[face], {0.0, 10.0, 20.0}),
                    "vector face halo exchange failed on the left interface");
            }
            if (i == detail::meshData(local).owned_i_end - 1 && parallel.rank + 1 < parallel.size) {
                const Index face = detail::meshData(local).cell_faces[static_cast<std::size_t>(cell)]
                    [static_cast<std::size_t>(Side::XMax)];
                require(
                    near(detail::fieldData(face_field)[face], 0.0),
                    "face-centred owner-authoritative exchange changed the owner");
                require(
                    near(detail::fieldData(face_vector)[face], {0.0, 10.0, 20.0}),
                    "vector face owner-authoritative exchange changed the owner");
            }
        }

        ScalarField affine(local, FieldLocation::Cell, "affine");
        affine.setBoundary(
            static_cast<Index>(Side::XMin),
            BoundaryCondition<double>::fixedValue(1.0));
        affine.setBoundary(
            static_cast<Index>(Side::XMax),
            BoundaryCondition<double>::fixedValue(5.0));
        for (Index cell : detail::meshData(local).owned_cells) {
            detail::fieldData(affine)[cell] = 2.0 *
                detail::meshData(local).cell_centres[static_cast<std::size_t>(cell)].x + 1.0;
        }
        halo.exchange(affine);
        VectorField affine_gradient(
            local, FieldLocation::Cell, "gradAffine");
        gradient(affine, affine_gradient, GradientMethod::GreenGauss);
        ScalarField affine_faces(
            local, FieldLocation::Face, "affineFaces");
        interpolate(affine, affine_faces);
        for (Index cell : detail::meshData(local).owned_cells) {
            require(
                near(detail::fieldData(affine_gradient)[cell], {2.0, 0.0, 0.0}, 1e-10),
                "distributed Green-Gauss gradient is incorrect");
        }
        for (Index face = 0; face < local.faceCount(); ++face) {
            const auto f = static_cast<std::size_t>(face);
            if (!detail::isOwned(local, detail::meshData(local).face_owner[f]) &&
                (detail::meshData(local).face_neighbour[f] == invalid_index ||
                 !detail::isOwned(local, detail::meshData(local).face_neighbour[f]))) {
                continue;
            }
            require(
                near(
                    detail::fieldData(affine_faces)[face],
                    2.0 * detail::meshData(local).face_centres[f].x + 1.0,
                    1e-10),
                "distributed interpolation is incorrect");
        }

        VectorField constant_velocity(
            local, FieldLocation::Cell, "constantVelocity", {1.0, 0.2, -0.1});
        for (Index patch = 0;
             patch < static_cast<Index>(detail::meshData(global).patches.size()); ++patch) {
            constant_velocity.setBoundary(
                patch,
                BoundaryCondition<Vec3>::fixedValue({1.0, 0.2, -0.1}));
        }
        halo.exchange(constant_velocity);
        ScalarField face_flux(local, FieldLocation::Face, "faceFlux");
        ScalarField velocity_divergence(
            local, FieldLocation::Cell, "divVelocity");
        flux(constant_velocity, face_flux);
        divergence(face_flux, velocity_divergence);
        for (Index cell : detail::meshData(local).owned_cells) {
            require(
                std::abs(detail::fieldData(velocity_divergence)[cell]) < 1e-11,
                "distributed constant-field divergence is not zero");
        }

        // 每个通用方程算子均在局部 owned 行和已同步的 ghost 输入上测试，
        // 不依赖 SIMPLE。
        ScalarField previous(local, FieldLocation::Cell, "previous", 0.0);
        ScalarField older(local, FieldLocation::Cell, "older", 0.0);
        for (Index cell : detail::meshData(local).owned_cells) {
            detail::fieldData(previous)[cell] = 1.0 + detail::meshData(local).cell_centres[static_cast<std::size_t>(cell)].x;
            detail::fieldData(older)[cell] = detail::fieldData(previous)[cell] - 0.1;
        }
        halo.exchange(previous);
        halo.exchange(older);
        ScalarDiscreteEquation time_equation(local);
        addTimeDerivative(time_equation, previous, 0.2, 2.0, TimeMethod::Euler);
        addTimeDerivative(
            time_equation, previous, 0.2, 2.0, TimeMethod::BDF2, &older);
        ScalarDiscreteEquation convection_equation(local);
        addConvection(convection_equation, face_flux, affine, ConvectionMethod::Upwind);
        ScalarDiscreteEquation corrected_equation(local);
        addDiffusion(
            corrected_equation, 1.0, affine,
            GradientMethod::LeastSquares, DiffusionMethod::Corrected);
        SparseAssembly generic_assembly(local);
        generic_assembly.update(time_equation);
        generic_assembly.update(convection_equation);
        generic_assembly.update(corrected_equation);
        require(
            generic_assembly.matrix().rows() == detail::ownedCellCount(local),
            "generic distributed operators generated ghost matrix rows");

        ScalarDiscreteEquation equation(local);
        std::fill(equation.diagonal.begin(), equation.diagonal.end(), 2.0);
        SparseAssembly assembly(local);
        assembly.update(equation);
        require(
            assembly.matrix().rows() == detail::ownedCellCount(local) &&
            assembly.matrix().cols() == detail::ownedCellCount(local),
            "distributed sparse assembly contains ghost rows");
        Eigen::VectorXd source;
        for (Index cell = 0; cell < local.cellCount(); ++cell) {
            equation.source[static_cast<std::size_t>(cell)] =
                detail::globalCellId(local, cell);
        }
        assembleSource(equation, source);
        for (Index cell : detail::meshData(local).owned_cells) {
            require(
                near(source[detail::ownedIndex(local, cell)], detail::globalCellId(local, cell)),
                "distributed source assembly changed owned ordering");
        }

        ScalarField diffusion(local, FieldLocation::Cell, "diffusion");
        diffusion.setBoundary(
            static_cast<Index>(Side::XMin),
            BoundaryCondition<double>::fixedValue(0.0));
        diffusion.setBoundary(
            static_cast<Index>(Side::XMax),
            BoundaryCondition<double>::fixedValue(2.0));
        ScalarDiscreteEquation diffusion_equation(local);
        addDiffusion(
            diffusion_equation, 1.0, diffusion,
            GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
        SparseAssembly diffusion_assembly(local);
        diffusion_assembly.update(diffusion_equation);
        Eigen::VectorXd diffusion_source;
        assembleSource(diffusion_equation, diffusion_source);
        LinearSolverConfig solver_config;
        solver_config.solver = LinearSolverType::ConjugateGradient;
        solver_config.preconditioner = PreconditionerType::IncompleteCholesky;
        solver_config.absolute_tolerance = 1e-13;
        solver_config.relative_tolerance = 1e-11;
        solver_config.max_iterations = 200;
        DistributedLinearSolver solver(local, parallel, solver_config);
        solver.compute(diffusion_assembly.matrix(), diffusion_equation);
        Eigen::VectorXd diffusion_solution;
        const SolveResult solve_result = solver.solve(
            diffusion_source, diffusion_solution);
        require(solve_result.converged(), "distributed PCG diffusion solve failed");
        double local_error = 0.0;
        for (Index cell : detail::meshData(local).owned_cells) {
            local_error = std::max(
                local_error,
                std::abs(
                    diffusion_solution[detail::ownedIndex(local, cell)] -
                    detail::meshData(local).cell_centres[static_cast<std::size_t>(cell)].x));
        }
        double global_error = 0.0;
        parallel.maximum(&local_error, &global_error, 1);
        require(global_error < 1e-10, "distributed diffusion solution is incorrect");

        const std::array<Index, 3> skew_dimensions{{8, 5, 5}};
        std::vector<Vec3> skew_points;
        for (Index k = 0; k <= skew_dimensions[2]; ++k) {
            for (Index j = 0; j <= skew_dimensions[1]; ++j) {
                for (Index i = 0; i <= skew_dimensions[0]; ++i) {
                    skew_points.push_back({
                        static_cast<double>(i) + 0.25 * j + 0.10 * k,
                        static_cast<double>(j) + 0.15 * k,
                        static_cast<double>(k),
                    });
                }
            }
        }
        const Mesh skew_global = Mesh::structured(
            skew_dimensions, std::move(skew_points));
        const Mesh skew = decompose(skew_global, parallel);
        HaloExchange skew_halo(skew, parallel);
        ScalarField skew_linear(skew, FieldLocation::Cell, "skewLinear");
        for (Index cell : detail::meshData(skew).owned_cells) {
            const Vec3& point =
                detail::meshData(skew).cell_centres[static_cast<std::size_t>(cell)];
            detail::fieldData(skew_linear)[cell] =
                2.0 * point.x - 3.0 * point.y + 0.5 * point.z + 1.0;
        }
        skew_halo.exchange(skew_linear);
        VectorField skew_gradient(skew, FieldLocation::Cell, "skewGradient");
        ScalarField skew_laplacian(skew, FieldLocation::Cell, "skewLaplacian");
        gradient(skew_linear, skew_gradient, GradientMethod::LeastSquares);
        laplacian(
            skew_linear, skew_laplacian,
            GradientMethod::LeastSquares, DiffusionMethod::Corrected);
        for (Index cell : detail::meshData(skew).owned_cells) {
            const Index global_id = detail::globalCellId(skew, cell);
            const Index global_i = global_id % skew_dimensions[0];
            const Index global_j =
                (global_id / skew_dimensions[0]) % skew_dimensions[1];
            const Index global_k = global_id /
                (skew_dimensions[0] * skew_dimensions[1]);
            if ((global_i == 3 || global_i == 4) &&
                global_j == 2 && global_k == 2) {
                require(
                    near(detail::fieldData(skew_gradient)[cell], {2.0, -3.0, 0.5}, 1e-10),
                    "MPI least-squares gradient failed at a partition face");
                require(
                    std::abs(detail::fieldData(skew_laplacian)[cell]) < 1e-10,
                    "MPI corrected Laplacian failed at a partition face");
            }
        }

        // 在分区交界面验证带偏斜修正的通用算子。所有输入只包含本地
        // owned+ghost 数据，算子内部不访问全局场。
        ScalarField skew_faces(skew, FieldLocation::Face, "skewFaces");
        interpolate(
            skew_linear, skew_faces,
            InterpolationMethod::Corrected, GradientMethod::LeastSquares);
        VectorField skew_velocity(skew, FieldLocation::Cell, "skewVelocity");
        for (Index cell : detail::meshData(skew).owned_cells) {
            const Vec3& point =
                detail::meshData(skew).cell_centres[static_cast<std::size_t>(cell)];
            detail::fieldData(skew_velocity)[cell] = {
                point.x + 2.0 * point.y - 0.5 * point.z,
                -point.x + 3.0 * point.z,
                0.25 * point.x - point.y + 2.0 * point.z,
            };
        }
        skew_halo.exchange(skew_velocity);
        ScalarField skew_flux(skew, FieldLocation::Face, "skewFlux");
        ScalarField skew_divergence(
            skew, FieldLocation::Cell, "skewDivergence");
        flux(
            skew_velocity, skew_flux,
            InterpolationMethod::Corrected, GradientMethod::LeastSquares);
        divergence(
            skew_velocity, skew_divergence,
            InterpolationMethod::Corrected, GradientMethod::LeastSquares);
        VectorField advecting_velocity(
            skew, FieldLocation::Cell, "advectingVelocity", {0.7, -0.2, 0.4});
        ScalarField advecting_flux(
            skew, FieldLocation::Face, "advectingFlux");
        flux(
            advecting_velocity, advecting_flux,
            InterpolationMethod::Corrected, GradientMethod::LeastSquares);
        ScalarDiscreteEquation skew_convection(skew);
        addConvection(
            skew_convection, advecting_flux, skew_linear,
            ConvectionMethod::Central, InterpolationMethod::Corrected,
            GradientMethod::LeastSquares);

        double local_operator_error = 0.0;
        for (Index cell : detail::meshData(skew).owned_cells) {
            const Index global_id = detail::globalCellId(skew, cell);
            const Index global_i = global_id % skew_dimensions[0];
            const Index global_j =
                (global_id / skew_dimensions[0]) % skew_dimensions[1];
            const Index global_k = global_id /
                (skew_dimensions[0] * skew_dimensions[1]);
            if ((global_i != 3 && global_i != 4) ||
                global_j != 2 || global_k != 2) {
                continue;
            }
            local_operator_error = std::max(
                local_operator_error,
                std::abs(detail::fieldData(skew_divergence)[cell] - 3.0));
            for (Index face : detail::meshData(skew).cell_faces[static_cast<std::size_t>(cell)]) {
                const auto f = static_cast<std::size_t>(face);
                const Vec3& point = detail::meshData(skew).face_centres[f];
                local_operator_error = std::max(
                    local_operator_error,
                    std::abs(
                        detail::fieldData(skew_faces)[face] -
                        (2.0 * point.x - 3.0 * point.y +
                         0.5 * point.z + 1.0)));
                const Vec3 exact_velocity{
                    point.x + 2.0 * point.y - 0.5 * point.z,
                    -point.x + 3.0 * point.z,
                    0.25 * point.x - point.y + 2.0 * point.z,
                };
                local_operator_error = std::max(
                    local_operator_error,
                    std::abs(
                        detail::fieldData(skew_flux)[face] -
                        dot(exact_velocity, detail::meshData(skew).face_area_vectors[f])));
            }
            double row =
                skew_convection.diagonal[static_cast<std::size_t>(cell)] *
                    detail::fieldData(skew_linear)[cell] -
                skew_convection.source[static_cast<std::size_t>(cell)];
            for (Index face : detail::meshData(skew).cell_faces[static_cast<std::size_t>(cell)]) {
                const auto f = static_cast<std::size_t>(face);
                const Index owner = detail::meshData(skew).face_owner[f];
                const Index neighbour = detail::meshData(skew).face_neighbour[f];
                row += owner == cell
                    ? skew_convection.upper[f] * detail::fieldData(skew_linear)[neighbour]
                    : skew_convection.lower[f] * detail::fieldData(skew_linear)[owner];
            }
            const double expected =
                detail::meshData(skew).cell_volumes[static_cast<std::size_t>(cell)] * 2.2;
            local_operator_error = std::max(
                local_operator_error, std::abs(row - expected));
        }
        double global_operator_error = 0.0;
        parallel.maximum(
            &local_operator_error, &global_operator_error, 1);
        require(
            global_operator_error < 1e-10,
            "MPI corrected operators failed at a partition face");

        ScalarField pressure(local, FieldLocation::Cell, "p");
        const std::filesystem::path directory = argc > 1
            ? argv[1] : "build/mpi-output";
        const auto time = directory / "domain";
        writeOwnedFieldCsv(time, vector, parallel);
        writeOwnedFieldCsv(time, pressure, parallel);
        writeOwnedResultMetadata(time, local, parallel, "domain", {
            {"vector", "vector", FieldLocation::Cell}, {"p", "scalar", FieldLocation::Cell}});
        parallel.barrier();
        std::ostringstream rank_directory;
        rank_directory << "rank-" << std::setfill('0') << std::setw(4) << parallel.rank;
        std::ifstream input(
            time / rank_directory.str() / "vector.csv");
        require(static_cast<bool>(input), "parallel output file is missing");
        std::string line;
        Index lines = 0;
        while (std::getline(input, line)) {
            ++lines;
        }
        require(
            lines == detail::ownedCellCount(local) + 1,
            "parallel output contains ghost cells or misses owned cells");

        if (parallel.rank == 0) {
            std::cout << "parallel_domain_test: global_cells="
                      << global.cellCount() << " ranks=" << parallel.size
                      << " ghost_layers=" << detail::meshData(local).ghost_layers
                      << " pcg_iterations=" << solve_result.iterations << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "rank " << parallel.rank << ": " << error.what() << '\n';
        const int abort_status = MPI_Abort(parallel.communicator, 1);
        if (abort_status != MPI_SUCCESS) return 1;
    }
    return MPI_Finalize() == MPI_SUCCESS ? 0 : 1;
}
