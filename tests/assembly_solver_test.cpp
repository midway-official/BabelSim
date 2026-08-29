#include "babelsim/assembly.h"
#include "babelsim/linear_solver.h"
#include "babelsim/operators.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    auto patches = defaultPatches();
    patches[static_cast<std::size_t>(Side::YMin)].kind = PatchKind::Symmetry;
    patches[static_cast<std::size_t>(Side::YMax)].kind = PatchKind::Symmetry;
    patches[static_cast<std::size_t>(Side::ZMin)].kind = PatchKind::Symmetry;
    patches[static_cast<std::size_t>(Side::ZMax)].kind = PatchKind::Symmetry;
    const Mesh mesh = Mesh::cartesian(
        {8, 1, 1}, {0, 0, 0}, {1, 1, 1}, patches);

    ScalarField phi(mesh, FieldLocation::Cell, "phi", 0.0);
    phi.setBoundary(
        static_cast<Index>(Side::XMin),
        BoundaryCondition<double>::fixedValue(0.0));
    phi.setBoundary(
        static_cast<Index>(Side::XMax),
        BoundaryCondition<double>::fixedValue(1.0));
    for (Side side : {Side::YMin, Side::YMax, Side::ZMin, Side::ZMax}) {
        phi.setBoundary(
            static_cast<Index>(side),
            BoundaryCondition<double>::symmetry());
    }

    ScalarEquation equation(mesh);
    addDiffusion(
        equation, 1.0, phi, GradientMethod::GreenGauss,
        DiffusionMethod::Orthogonal);
    const LinearSystem system = assemble(equation);
    const double symmetry_error =
        (system.A - Eigen::SparseMatrix<double>(system.A.transpose())).norm();
    require(symmetry_error < 1e-13, "diffusion assembly is not symmetric");

    Eigen::VectorXd solution;
    LinearSolverConfig config;
    config.solver = LinearSolverType::ConjugateGradient;
    config.preconditioner = PreconditionerType::IncompleteCholesky;
    config.absolute_tolerance = 1e-14;
    config.relative_tolerance = 1e-12;
    const SolveResult result = solve(system.A, system.b, solution, config);
    require(result.converged(), "CG did not solve the diffusion equation");

    double maximum_error = 0.0;
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        const double exact = mesh.cell_centres[static_cast<std::size_t>(cell)].x;
        maximum_error = std::max(
            maximum_error,
            std::abs(solution[static_cast<Eigen::Index>(cell)] - exact));
    }
    require(maximum_error < 1e-11, "assembled linear diffusion solution is incorrect");

    PreparedLinearSolver prepared(config);
    prepared.compute(system.A);
    Eigen::VectorXd first;
    Eigen::VectorXd second;
    const SolveResult first_result = prepared.solve(system.b, first);
    const SolveResult second_result = prepared.solve(2.0 * system.b, second);
    require(
        first_result.converged() && second_result.converged() &&
            (first - solution).norm() < 1e-12 &&
            (second - 2.0 * solution).norm() < 1e-12,
        "prepared linear solver failed to reuse a factorization");

    VectorEquation vector_equation(mesh);
    vector_equation.diagonal = equation.diagonal;
    vector_equation.upper = equation.upper;
    vector_equation.lower = equation.lower;
    require(
        (assembleMatrix(vector_equation) - system.A).norm() < 1e-14,
        "scalar and segregated-vector matrix assembly differ");

    SparseAssembly cached_assembly(mesh);
    cached_assembly.update(equation);
    require(
        (cached_assembly.matrix() - system.A).norm() < 1e-14,
        "precomputed sparse assembly differs from triplet assembly");
    ScalarEquation rescaled_equation = equation;
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
        (cached_assembly.matrix() - 2.0 * system.A).norm() < 1e-14,
        "precomputed sparse assembly did not update coefficient values");

    ScalarField face_diffusivity(
        mesh, FieldLocation::Face, "faceDiffusivity", 3.0);
    ScalarEquation field_diffusion(mesh);
    ScalarEquation constant_diffusion(mesh);
    addDiffusion(
        field_diffusion, face_diffusivity, phi,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    addDiffusion(
        constant_diffusion, 3.0, phi,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    require(
        (assembleMatrix(field_diffusion) -
         assembleMatrix(constant_diffusion)).norm() < 1e-14 &&
            field_diffusion.source == constant_diffusion.source,
        "face-centred scalar diffusivity differs from a constant coefficient");

    VectorField vector(mesh, FieldLocation::Cell, "vector");
    VectorEquation field_vector_diffusion(mesh);
    VectorEquation constant_vector_diffusion(mesh);
    addDiffusion(
        field_vector_diffusion, face_diffusivity, vector,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    addDiffusion(
        constant_vector_diffusion, 3.0, vector,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    double vector_source_error = 0.0;
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        vector_source_error = std::max(
            vector_source_error,
            norm(
                field_vector_diffusion.source[static_cast<std::size_t>(cell)] -
                constant_vector_diffusion.source[static_cast<std::size_t>(cell)]));
    }
    require(
        (assembleMatrix(field_vector_diffusion) -
         assembleMatrix(constant_vector_diffusion)).norm() < 1e-14 &&
            vector_source_error < 1e-14,
        "face-centred vector diffusivity differs from a constant coefficient");

    ScalarField scalar_neumann(mesh, FieldLocation::Cell, "scalarNeumann");
    scalar_neumann.setBoundary(
        static_cast<Index>(Side::XMax),
        BoundaryCondition<double>::fixedGradient(2.0));
    ScalarEquation scalar_neumann_equation(mesh);
    addDiffusion(
        scalar_neumann_equation, 3.0, scalar_neumann,
        GradientMethod::GreenGauss, DiffusionMethod::Orthogonal);
    const Index last = mesh.cellId(7, 0, 0);
    require(
        near(scalar_neumann_equation.source[static_cast<std::size_t>(last)], 6.0),
        "scalar outward Neumann flux has the wrong equation sign");

    VectorField vector_neumann(mesh, FieldLocation::Cell, "vectorNeumann");
    vector_neumann.setBoundary(
        static_cast<Index>(Side::XMax),
        BoundaryCondition<Vec3>::fixedGradient({2.0, -1.0, 0.5}));
    VectorEquation vector_neumann_equation(mesh);
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
