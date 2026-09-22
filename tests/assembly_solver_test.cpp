#include "internal/mesh_access.h"
#include "babelsim/assembly.h"
#include "babelsim/linear_solver.h"
#include "babelsim/operators.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    const Mesh split_mesh = makeSplitInterfaceMesh();
    require(split_mesh.cellCount() == 2 && split_mesh.faceCount() == 12,
            "split-interface polyhedral fixture has the wrong topology");
    ScalarDiscreteEquation repeated_equation(split_mesh);
    repeated_equation.diagonal = {5.0, 7.0};
    const Index first_interface = 1;
    const Index second_interface = 2;
    repeated_equation.upper[static_cast<std::size_t>(first_interface)] = 2.0;
    repeated_equation.upper[static_cast<std::size_t>(second_interface)] = 11.0;
    repeated_equation.lower[static_cast<std::size_t>(first_interface)] = 3.0;
    repeated_equation.lower[static_cast<std::size_t>(second_interface)] = 13.0;
    SparseAssembly repeated_assembly(split_mesh);
    repeated_assembly.update(repeated_equation);
    require(
        near(repeated_assembly.matrix().coeff(0, 0), 5.0) &&
            near(repeated_assembly.matrix().coeff(1, 1), 7.0) &&
            near(repeated_assembly.matrix().coeff(0, 1), 13.0) &&
            near(repeated_assembly.matrix().coeff(1, 0), 16.0),
        "repeated subface coupling was not reduced into one matrix entry");
    Eigen::MatrixXd dense_oracle(2, 2);
    dense_oracle << 5.0, 13.0, 16.0, 7.0;
    require(
        (repeated_assembly.matrix().toDense() - dense_oracle).norm() < 1e-14,
        "polyhedral sparse matrix differs from the independent dense oracle");
    LinearSolverConfig replay_config;
    replay_config.solver = LinearSolverType::BiCGSTAB;
    replay_config.preconditioner = PreconditionerType::None;
    replay_config.absolute_tolerance = 1e-14;
    replay_config.relative_tolerance = 1e-14;
    PreparedLinearSolver replay_solver(replay_config);
    replay_solver.compute(repeated_assembly.matrix());
    Eigen::VectorXd frozen_x(2);
    frozen_x << 1.25, -0.75;
    const Eigen::VectorXd frozen_b = dense_oracle * frozen_x;
    Eigen::VectorXd replay_x;
    const SolveResult replay_result = replay_solver.solve(frozen_b, replay_x);
    require(
        replay_result.converged() && (replay_x - frozen_x).norm() < 1e-12,
        "frozen polyhedral A/b/x0 algebra replay did not reproduce x0");
    repeated_equation.diagonal = {17.0, 19.0};
    repeated_equation.upper[static_cast<std::size_t>(first_interface)] = 1.0;
    repeated_equation.upper[static_cast<std::size_t>(second_interface)] = -2.0;
    repeated_equation.lower[static_cast<std::size_t>(first_interface)] = 4.0;
    repeated_equation.lower[static_cast<std::size_t>(second_interface)] = 8.0;
    repeated_assembly.update(repeated_equation);
    require(
        near(repeated_assembly.matrix().coeff(0, 0), 17.0) &&
            near(repeated_assembly.matrix().coeff(1, 1), 19.0) &&
            near(repeated_assembly.matrix().coeff(0, 1), -1.0) &&
            near(repeated_assembly.matrix().coeff(1, 0), 12.0),
        "repeated subface assembly accumulated stale coefficient values");

    auto patches = boxPatches();
    patches[static_cast<std::size_t>(2)].kind = PatchKind::Symmetry;
    patches[static_cast<std::size_t>(3)].kind = PatchKind::Symmetry;
    patches[static_cast<std::size_t>(4)].kind = PatchKind::Symmetry;
    patches[static_cast<std::size_t>(5)].kind = PatchKind::Symmetry;
    const Mesh mesh = makeHexBox(
        {8, 1, 1}, {0, 0, 0}, {1, 1, 1}, patches);

    ScalarField phi(mesh, FieldLocation::Cell, "phi", 0.0);
    phi.setBoundary(
        static_cast<Index>(0),
        BoundaryCondition<double>::fixedValue(0.0));
    phi.setBoundary(
        static_cast<Index>(1),
        BoundaryCondition<double>::fixedValue(1.0));
    for (Index side : {2, 3, 4, 5}) {
        phi.setBoundary(
            static_cast<Index>(side),
            BoundaryCondition<double>::symmetry());
    }

    ScalarDiscreteEquation equation(mesh);
    addDiffusion(
        equation, 1.0, phi, GradientMethod::GreenGauss,
        DiffusionMethod::Orthogonal);
    SparseAssembly assembly(mesh);
    assembly.update(equation);
    Eigen::VectorXd source;
    assembleSource(equation, source);
    const auto& matrix = assembly.matrix();
    const double symmetry_error =
        (matrix - Eigen::SparseMatrix<double>(matrix.transpose())).norm();
    require(symmetry_error < 1e-13, "diffusion assembly is not symmetric");

    Eigen::VectorXd solution;
    LinearSolverConfig config;
    config.solver = LinearSolverType::ConjugateGradient;
    config.preconditioner = PreconditionerType::IncompleteCholesky;
    config.absolute_tolerance = 1e-14;
    config.relative_tolerance = 1e-12;
    PreparedLinearSolver direct_solver(config);
    direct_solver.compute(matrix);
    const SolveResult result = direct_solver.solve(source, solution);
    require(result.converged(), "CG did not solve the diffusion equation");

    double maximum_error = 0.0;
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        const double exact = detail::meshData(mesh).cell_centres[static_cast<std::size_t>(cell)].x;
        maximum_error = std::max(
            maximum_error,
            std::abs(solution[static_cast<Eigen::Index>(cell)] - exact));
    }
    require(maximum_error < 1e-11, "assembled linear diffusion solution is incorrect");

    LinearSolverConfig no_preconditioner_config = config;
    no_preconditioner_config.preconditioner = PreconditionerType::None;
    Eigen::VectorXd no_preconditioner_solution;
    PreparedLinearSolver no_preconditioner_solver(no_preconditioner_config);
    no_preconditioner_solver.compute(matrix);
    const SolveResult no_preconditioner_result =
        no_preconditioner_solver.solve(source, no_preconditioner_solution);
    require(
        no_preconditioner_result.converged() &&
            (no_preconditioner_solution - solution).norm() < 1e-11 &&
            no_preconditioner_result.performance.preconditioner_applications == 0 &&
            no_preconditioner_result.performance.preconditioner_apply_seconds == 0.0,
        "unpreconditioned CG did not remain an identity preconditioner");

    PreparedLinearSolver prepared(config);
    prepared.compute(matrix);
    Eigen::VectorXd first;
    Eigen::VectorXd second;
    const SolveResult first_result = prepared.solve(source, first);
    const SolveResult second_result = prepared.solve(2.0 * source, second);
    require(
        first_result.converged() && second_result.converged() &&
            (first - solution).norm() < 1e-12 &&
            (second - 2.0 * solution).norm() < 1e-12,
        "prepared linear solver failed to reuse a factorization");

    // AMG 只接收代数系统并作为 Krylov 预条件器；不作为独立线性求解器。
    // 将粗网格阈值压低，确保该小系统也实际建立多层而非退化为一次直接分解。
    LinearSolverConfig amg_config = config;
    amg_config.preconditioner = PreconditionerType::AlgebraicMultigrid;
    amg_config.amg_coarse_size = 2;
    amg_config.amg_smoothing_steps = 2;
    amg_config.max_iterations = 100;
    Eigen::VectorXd amg_cg_solution;
    PreparedLinearSolver amg_cg_solver(amg_config);
    amg_cg_solver.compute(matrix);
    const SolveResult amg_cg_result = amg_cg_solver.solve(source, amg_cg_solution);
    require(
        amg_cg_result.converged() && (amg_cg_solution - solution).norm() < 1e-10,
        "AMG-preconditioned CG did not solve diffusion");

    amg_config.solver = LinearSolverType::BiCGSTAB;
    Eigen::VectorXd bicgstab_solution;
    PreparedLinearSolver bicgstab_solver(amg_config);
    bicgstab_solver.compute(matrix);
    const SolveResult bicgstab_result = bicgstab_solver.solve(source, bicgstab_solution);
    require(
        bicgstab_result.converged() &&
            (bicgstab_solution - solution).norm() < 1e-10,
        "AMG-preconditioned BiCGSTAB did not solve diffusion");
    require(
        bicgstab_result.performance.sparse_matvecs > 0 &&
            bicgstab_result.performance.preconditioner_applications > 0 &&
            bicgstab_result.performance.sparse_matvec_seconds >= 0.0 &&
            bicgstab_result.performance.preconditioner_apply_seconds >= 0.0,
        "serial Krylov performance counters were not collected by the actual kernels");

    // AMG 是预条件器时可短期复用上一轮层级；Krylov matvec 仍使用新矩阵，
    // 因而复用只影响速度和迭代数，不能改变线性系统的解。
    amg_config.amg_refresh_interval = 4;
    PreparedLinearSolver cached_amg(amg_config);
    cached_amg.compute(matrix);
    Eigen::VectorXd cached_amg_solution;
    cached_amg.factorize(2.0 * matrix);
    require(
        cached_amg.solve(2.0 * source, cached_amg_solution).converged() &&
            (cached_amg_solution - solution).norm() < 1e-10,
        "reused AMG preconditioner did not solve the updated system");

    PreparedLinearSolver prepared_amg(amg_config);
    prepared_amg.compute(matrix);
    Eigen::VectorXd prepared_amg_solution;
    require(
        prepared_amg.solve(source, prepared_amg_solution).converged() &&
            (prepared_amg_solution - solution).norm() < 1e-10,
        "prepared AMG did not solve diffusion");
    prepared_amg.factorize(2.0 * matrix);
    Eigen::VectorXd refactorized_amg_solution;
    require(
        prepared_amg.solve(2.0 * source, refactorized_amg_solution).converged() &&
            (refactorized_amg_solution - solution).norm() < 1e-10,
        "AMG factorization did not reuse its hierarchy");

    VectorDiscreteEquation vector_equation(mesh);
    vector_equation.diagonal = equation.diagonal;
    vector_equation.upper = equation.upper;
    vector_equation.lower = equation.lower;
    SparseAssembly vector_assembly(mesh);
    vector_assembly.update(vector_equation);
    require(
        (vector_assembly.matrix() - matrix).norm() < 1e-14,
        "scalar and segregated-vector matrix assembly differ");

    SparseAssembly cached_assembly(mesh);
    cached_assembly.update(equation);
    require(
        (cached_assembly.matrix() - matrix).norm() < 1e-14,
        "precomputed sparse assembly differs from triplet assembly");
    ScalarDiscreteEquation rescaled_equation = equation;
    for (double& value : rescaled_equation.diagonal) {
        value *= 2.0;
    }
    for (double& value : rescaled_equation.upper) {
        value *= 2.0;
    }
    for (double& value : rescaled_equation.lower) {
        value *= 2.0;
    }
    cached_assembly.update(rescaled_equation);
    require(
        (cached_assembly.matrix() - 2.0 * matrix).norm() < 1e-14,
        "precomputed sparse assembly did not update coefficient values");

    ScalarField face_diffusivity(
        mesh, FieldLocation::Face, "faceDiffusivity", 3.0);
    ScalarDiscreteEquation field_diffusion(mesh);
    ScalarDiscreteEquation constant_diffusion(mesh);
    addDiffusion(
        field_diffusion, face_diffusivity, phi,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    addDiffusion(
        constant_diffusion, 3.0, phi,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    SparseAssembly field_diffusion_assembly(mesh);
    SparseAssembly constant_diffusion_assembly(mesh);
    field_diffusion_assembly.update(field_diffusion);
    constant_diffusion_assembly.update(constant_diffusion);
    require(
        (field_diffusion_assembly.matrix() -
         constant_diffusion_assembly.matrix()).norm() < 1e-14 &&
            field_diffusion.source == constant_diffusion.source,
        "face-centred scalar diffusivity differs from a constant coefficient");

    VectorField vector(mesh, FieldLocation::Cell, "vector");
    VectorDiscreteEquation field_vector_diffusion(mesh);
    VectorDiscreteEquation constant_vector_diffusion(mesh);
    addDiffusion(
        field_vector_diffusion, face_diffusivity, vector,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    addDiffusion(
        constant_vector_diffusion, 3.0, vector,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    SparseAssembly field_vector_diffusion_assembly(mesh);
    SparseAssembly constant_vector_diffusion_assembly(mesh);
    field_vector_diffusion_assembly.update(field_vector_diffusion);
    constant_vector_diffusion_assembly.update(constant_vector_diffusion);
    double vector_source_error = 0.0;
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        vector_source_error = std::max(
            vector_source_error,
            norm(
                field_vector_diffusion.source[static_cast<std::size_t>(cell)] -
                constant_vector_diffusion.source[static_cast<std::size_t>(cell)]));
    }
    require(
        (field_vector_diffusion_assembly.matrix() -
         constant_vector_diffusion_assembly.matrix()).norm() < 1e-14 &&
            vector_source_error < 1e-14,
        "face-centred vector diffusivity differs from a constant coefficient");

    ScalarField scalar_neumann(mesh, FieldLocation::Cell, "scalarNeumann");
    scalar_neumann.setBoundary(
        static_cast<Index>(1),
        BoundaryCondition<double>::fixedGradient(2.0));
    ScalarDiscreteEquation scalar_neumann_equation(mesh);
    addDiffusion(
        scalar_neumann_equation, 3.0, scalar_neumann,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    const Index last = hexCellIndex(7, 0, 0, 8, 1);
    require(
        near(scalar_neumann_equation.source[static_cast<std::size_t>(last)], 6.0),
        "scalar outward Neumann flux has the wrong equation sign");

    VectorField vector_neumann(mesh, FieldLocation::Cell, "vectorNeumann");
    vector_neumann.setBoundary(
        static_cast<Index>(1),
        BoundaryCondition<Vec3>::fixedGradient({2.0, -1.0, 0.5}));
    VectorDiscreteEquation vector_neumann_equation(mesh);
    addDiffusion(
        vector_neumann_equation, 3.0, vector_neumann,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    require(
        near(
            vector_neumann_equation.source[static_cast<std::size_t>(last)],
            {6.0, -3.0, 1.5}),
        "vector outward Neumann flux has the wrong equation sign");

    std::cout << "assembly_solver_test: iterations=" << result.iterations
              << " residual=" << result.final_residual
              << " max_error=" << maximum_error << '\n';
}
