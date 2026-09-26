#include "babelsim/equ.h"
#include "babelsim/geometry.h"
#include "babelsim/runtime.h"
#include "babelsim/parallel.h"
#include "babelsim/mpi_support.h"
#include "internal/field_access.h"
#include "internal/petsc_session.h"
#include "test_util.h"
#include <iostream>

using namespace babelsim;

void fluxContext(const Mesh& mesh) {
    const auto options = testEquationControl("flux").spatial;
    VectorField u(mesh, FieldLocation::Cell, "U", Vec3{1,0,0});
    u.boundary("plus_x") = BoundaryCondition<Vec3>::inletOutlet(Vec3{-2,0,0});
    ScalarField context(mesh, FieldLocation::Face, "phi", -1.0);
    u.setBoundaryFlux(context);
    const auto phi = math::flux(u, options);
    const auto trace = math::interpolate(u, options);
    ScalarField output(mesh, FieldLocation::Face, "output", 123.0);
    math::evaluate(math::FaceFlux{u, options}, output);
    for (Index f = 0; f < mesh.faceCount(); ++f) {
        if (!mesh.boundaryFace(f) || mesh.patchName(mesh.boundaryPatch(f)) != "plus_x") continue;
        require(near(detail::fieldData(phi)[f], dot(Vec3{-2,0,0}, mesh.faceAreaVector(f))),
                "flux lost the input's inflow context");
        require(near(detail::fieldData(trace)[f], Vec3{-2,0,0}), "flux changed inflow to outflow");
        require(near(detail::fieldData(output)[f], detail::fieldData(phi)[f]),
                "flux depends on output buffer contents");
    }
    // No prior context: owner velocity determines direction, also for corrected
    // interpolation's gradient. Poisoning the output must not choose a branch.
    VectorField initial(mesh, FieldLocation::Cell, "initial", Vec3{-1,0,0});
    initial.boundary("plus_x") = BoundaryCondition<Vec3>::inletOutlet(Vec3{-2,0,0});
    math::evaluate(math::FaceFlux{initial, options}, output);
    for (Index f = 0; f < mesh.faceCount(); ++f)
        if (mesh.boundaryFace(f) && mesh.patchName(mesh.boundaryPatch(f)) == "plus_x")
            require(near(detail::fieldData(output)[f], dot(Vec3{-2,0,0}, mesh.faceAreaVector(f))),
                    "first flux did not initialize direction from owner velocity");
}

void frozenDiffusion(const Mesh& mesh) {
    const auto V = geometry::cellVolumes(mesh);
    for (auto method : {DiffusionMethod::Orthogonal, DiffusionMethod::Corrected,
                        DiffusionMethod::LimitedCorrected}) {
        for (bool variable : {false, true}) {
            ScalarField p(mesh, FieldLocation::Cell, "p");
            for (Index patch = 0; patch < mesh.patchCount(); ++patch)
                p.setBoundary(patch, fixedValue(0.0));
            p.evaluate([](Vec3 x) { return 1 + x.x*x.y; });
            ScalarField k(mesh, FieldLocation::Cell, "k");
            k.useCalculatedBoundary();
            k.evaluate([](Vec3 x) { return 1 + 0.1*x.x; });
            auto control = testEquationControl("diffusion");
            control.spatial.diffusion = method;
            control.linear.relative_tolerance = 1e-12;
            auto a = equ::createEquation(p, control);
            if (variable) equ::laplacian(a, k, -1);
            else equ::laplacian(a, 1., -1);
            // Copy, scale and add must preserve frozen flux contributions too.
            auto combined = a.copy();
            equ::scale(combined, 0.4); equ::add(combined, a, 0.6);
            k.fill(20.0);
            require(equ::solve(combined).converged(), "diffusion solve failed");
            const auto check = [&] {
                const auto matrixFlux = equ::apply(combined, p) - combined.rhs();
                const auto recovered = V * math::div(equ::faceFlux(combined, p));
                require(math::normL2(recovered - matrixFlux) < 1e-10,
                        "recovered flux no longer belongs to the assembled matrix");
            };
            check();
            p.evaluate([](Vec3 x) { return 0.1 + 0.2*x.y*x.y; });
            check(); // Not only at a converged solution or the assembly iterate.
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    try {
        std::vector<Vec3> vertices;
        for (int k=0;k<=1;++k) for (int j=0;j<=3;++j) for (int i=0;i<=4;++i)
            vertices.push_back({i+0.7*j, double(j), double(k)});
        const auto global = makeHexFromVertices({4,3,1}, vertices);
        const auto mesh = decompose(global, ParallelContext::world(), 3);
        auto runtime = RunTime::forMesh(mesh);
        fluxContext(mesh);
        frozenDiffusion(mesh);
        if (primaryProcess()) std::cout << "physics_contract_test: inflow context and frozen diffusion flux passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; MPI_Abort(MPI_COMM_WORLD, 1);
    }
    detail::finalizePetscSession();
    MPI_Finalize();
}
