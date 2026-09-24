#include "babelsim/runtime.h"
#include "babelsim/parallel.h"
#include "babelsim/mpi_support.h"
#include "babelsim/equ.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"
#include "internal/petsc_session.h"
#include "test_util.h"
#include <iostream>

using namespace babelsim;

void mixedBoundary() {
    const Mesh mesh = makeHexBox({1,1,1},{0,0,0},{1,1,1});
    RunTime time = RunTime::forMesh(mesh);
    ScalarField c(mesh, FieldLocation::Cell, "c"), phi(mesh, FieldLocation::Face, "phi");
    c.setBoundary(0, BoundaryCondition<double>::inletOutlet(1.0));
    c.setBoundary(1, fixedValue(0.0));
    detail::fieldData(phi)[detail::meshData(mesh).patches[0].faces.front()] = -1;
    detail::fieldData(phi)[detail::meshData(mesh).patches[1].faces.front()] = 1;
    VectorField u(mesh,FieldLocation::Cell,"u");
    u.setBoundary(0,BoundaryCondition<Vec3>::inletOutlet({1,2,3}));
    u.setBoundary(1,fixedValue(Vec3{}));
    {
        // div(phi,c) == laplacian(1.0,c)；通量先绑定，扩散面处理据此判断入流/出流。
        auto equation = testEquation(c);
        equ::div(equation, phi);
        equ::laplacian(equation, 1.0, -1.0);
        require(equ::solve(equation).converged(), "mixed scalar solve failed");
    }
    require(std::abs(detail::fieldData(c)[0]-0.6)<1e-12,"inflow diffusion contribution missing");
    {
        auto equation = testEquation(u);
        equ::div(equation, phi);
        equ::laplacian(equation, 1.0, -1.0);
        require(equ::solve(equation).converged(), "mixed vector solve failed");
    }
    require(norm(detail::fieldData(u)[0]-Vec3{0.6,1.2,1.8})<1e-12,"vector mixed boundary mismatch");
    phi.assignScaled(-1,phi);
    {
        auto equation = testEquation(c);
        equ::div(equation, phi);
        equ::laplacian(equation, 1.0, -1.0);
        require(equ::solve(equation).converged(), "reversed flow solve failed");
    }
    require(std::abs(detail::fieldData(c)[0])<1e-12,"outflow retained inlet Dirichlet constraint");
}

void boundaryAndAlgebra(const Mesh& mesh) {
    RunTime time = RunTime::forMesh(mesh);
    const auto options = testEquationControl("boundary").spatial;
    require(!std::isfinite(EquationResidual{1.0, std::numeric_limits<double>::infinity()}.relative()),
            "nonfinite equation scale was reported as zero relative residual");
    ScalarField x(mesh, FieldLocation::Cell, "x", 2.0);
    x.setBoundary(0, fixedValue(0.0));
    x.setBoundary(1, fixedGradient(2.0));
    ScalarField derived(mesh, FieldLocation::Cell, "derived"), face(mesh, FieldLocation::Face);
    derived.useCalculatedBoundary();
    derived.evaluate(x, [](double value) { return value * value; });
    derived.assignScaled(3.0, derived);
    derived.addScaled(2.0, x);
    math::evaluate(math::interpolate(derived, options), face);
    for (Index f : detail::meshData(mesh).owned_faces) {
        if (!mesh.boundaryFace(f) || mesh.boundaryPatch(f) != 0) continue;
        require(detail::fieldData(face)[f] == 0.0, "nonlinear field algebra lost the boundary trace");
    }
    x.evaluate(x, [](double) { return 7.0; });
    require(x.boundary(0).value == 0.0, "unknown value update modified its constraint");

    ScalarField a(mesh, FieldLocation::Cell, "a", 4.0), u(mesh, FieldLocation::Cell, "u", 0.0);
    // a(x)*u == 12：显式体源与局部隐式线性项。
    auto equation = testEquation(u);
    equ::reaction(equation, a, 1.0);
    equ::source(equation, 12.0);
    // 欠松弛只写进副本；原方程保持未松弛，才能用它判断真实残差。
    equ::Equation<double> relaxed_equation = equation.copy();
    equ::relax(relaxed_equation, u, 0.3);
    auto r = equ::solve(relaxed_equation, u);
    require(r.converged(), "linear Sp equation failed");
    require(diagnostics::relativeResidual(equation, u) > 1e-10,
            "relaxed-system success mistaken for original-equation convergence");
    r = equ::solve(equation, u);
    require(r.converged() && diagnostics::relativeResidual(equation, u) <= 1e-10,
            "original scalar equation residual is incorrect");
    for (Index cell : detail::meshData(mesh).owned_cells)
        require(near(detail::fieldData(u)[cell], 3.0), "implicit source has wrong sign/volume");
    u.fill(1.0);
    require(diagnostics::relativeResidual(equation, u) > 1e-10,
            "post-update residual failed to detect a changed solution");
    u.fill(0.0);
    auto overflow = testEquation(u);
    equ::reaction(overflow, 1.0);
    equ::source(overflow, 1e200);
    require(!equ::solve(overflow).converged(),
            "overflowed residual norm was reported as linear convergence");
    VectorField vector_unknown(mesh, FieldLocation::Cell, "vectorSp");
    auto vector_equation = testEquation(vector_unknown);
    equ::reaction(vector_equation, a, 1.0);
    equ::source(vector_equation, Vec3{12,-8,20});
    require(equ::solve(vector_equation).converged() &&
            diagnostics::relativeResidual(vector_equation, vector_unknown) <= 1e-10,
            "vector implicit source or original-equation residual is incorrect");
    for (Index cell : detail::meshData(mesh).owned_cells)
        require(norm(detail::fieldData(vector_unknown)[cell]-Vec3{3,-2,5}) < 1e-12,
                "vector implicit source component or volume is incorrect");

    TensorField tensor(mesh, FieldLocation::Cell, "tensor");
    tensor.useCalculatedBoundary();
    tensor.evaluate([](Vec3 p) {
        Tensor3 value;
        value[0] = {p.x, 2*p.y, 3*p.z};
        value[1] = {2*p.x, -p.y, p.z};
        value[2] = {-p.x, p.y, 2*p.z};
        return value;
    });
    VectorField div(mesh, FieldLocation::Cell, "divTensor");
    div.useCalculatedBoundary();
    math::evaluate(math::div(tensor, options), div);
    for (Index cell : detail::meshData(mesh).owned_cells)
        require(norm(detail::fieldData(div)[cell] - Vec3{6,2,2}) < 1e-11,
                "tensor divergence component convention or face trace is wrong");
    VectorField face_div(mesh, FieldLocation::Face, "faceDivTensor");
    math::evaluate(math::interpolate(div, options), face_div);
    for (Index f : detail::meshData(mesh).owned_faces)
        require(norm(detail::fieldData(face_div)[f] - Vec3{6,2,2}) < 1e-11,
                "differential result retained a stale calculated boundary trace");
}

void preconditioners(const Mesh& mesh) {
    const bool distributed = ParallelContext::world().size > 1;
    for (auto solver : {LinearSolverType::ConjugateGradient, LinearSolverType::BiCGSTAB}) {
        const auto factor = solver == LinearSolverType::ConjugateGradient && !distributed
            ? PreconditionerType::IncompleteCholesky
            : solver == LinearSolverType::ConjugateGradient
                ? PreconditionerType::BlockJacobi : PreconditionerType::Jacobi;
        for (auto pc : {PreconditionerType::None, factor, PreconditionerType::Hypre}) {
            for (double scale : {1e-40, 1e-20, 1.0, 1e20, 1e40}) {
                RunTime time = RunTime::forMesh(mesh);
                ScalarField u(mesh, FieldLocation::Cell, "u", 0.0);
                u.setBoundary(0, fixedValue(0.0)); u.setBoundary(1, fixedValue(1.0));
                VectorField velocity(mesh, FieldLocation::Face, "velocity", {1,0,0});
                ScalarField phi(mesh, FieldLocation::Face, "phi");
                auto assembledControl = testEquationControl(
                    "manufactured", InterpolationMethod::Linear, GradientMethod::LeastSquares,
                    ConvectionMethod::Central, DiffusionMethod::Orthogonal);
                math::evaluate(math::flux(velocity, assembledControl.spatial), phi);
                auto& cfg = assembledControl.linear;
                cfg.solver = solver; cfg.preconditioner = pc;
                cfg.absolute_tolerance = scale * 1e-13;
                cfg.relative_tolerance = 1e-10;
                cfg.max_iterations = 2000;
                cfg.amg_coarse_size = 4;
                auto assembled = equ::createEquation(u, assembledControl);
                // 通量上下文先于扩散装配绑定：-div(scale grad u) [+ scale div(phi u)] == source。
                if (solver == LinearSolverType::BiCGSTAB) equ::div(assembled, phi, scale);
                equ::laplacian(assembled, scale, -1.0);
                equ::source(assembled, solver == LinearSolverType::BiCGSTAB ? scale : 0.0);
                const auto result = equ::solve(assembled);
                if (!result.converged()) std::cerr << "solver=" << int(solver) << " pc=" << int(pc)
                    << " scale=" << scale << " status=" << int(result.status)
                    << " iterations=" << result.iterations
                    << " initial=" << result.initial_residual
                    << " final=" << result.final_residual
                    << " relative=" << result.relative_residual << '\n';
                require(result.converged(), "scale/preconditioner convergence contract failed");
                require(result.final_residual <= std::max(cfg.absolute_tolerance,
                    cfg.relative_tolerance * result.initial_residual), "false linear convergence");
                require(diagnostics::relativeResidual(assembled, u) <= 1e-8,
                        "preconditioned solve did not satisfy original equation");
                for (Index cell : detail::meshData(mesh).owned_cells)
                    require(std::abs(detail::fieldData(u)[cell]-mesh.cellCentre(cell).x) < 1e-8,
                            "manufactured transport solution is incorrect");
            }
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    const auto parallel = ParallelContext::world();
    try {
        const Mesh global = makeHexBox({16,2,1}, {0,0,0}, {1,1,1});
        bool rejected = false;
        try { (void)decompose(global, parallel, 2); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "ghost layers below three were accepted");
        const Mesh mesh = decompose(global, parallel, argc > 1 ? std::stoi(argv[1]) : 3);
        if (parallel.size == 1) mixedBoundary();
        boundaryAndAlgebra(mesh);
        preconditioners(mesh);
        if (parallel.rank == 0) std::cout << "numerical_contract_test: field traces, Sp, true residuals, tensor divergence; 30 solver/preconditioner/scale combinations passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; MPI_Abort(MPI_COMM_WORLD, 1);
    }
    detail::finalizePetscSession();
    MPI_Finalize();
}
