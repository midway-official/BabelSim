# BabelSim architecture and TaihoCFD migration map

## Audited starting point

`/home/midway/BabelSim` started with only a README and license. TaihoCFD is a
working C++17/Eigen/MPI two-dimensional incompressible finite-volume solver;
the following is its historical executable path used for the migration audit:

```text
case.cfg -> Taiho mesh reader -> x-decomposition -> SolverContext
         -> SIMPLE iteration
              momentum assembly (ddt + upwind + diffusion + grad(p))
              BiCGSTAB/ILUT velocity solves
              Rhie-Chow face velocity
              pressure-correction assembly and PCG/IC solve
              pressure, cell velocity, and face-flux correction
              continuity and field-change convergence tests
         -> owned-column output and MPI post-processing
```

The current TaihoCFD regression tests pass. The two-rank Poiseuille example
converges in 896 SIMPLE iterations; its outlet profile has approximately
`2.13e-4` RMS error against the fully-developed parabola. One- and two-rank
fields differ by roughly `1e-6` in relative maximum norm at the configured
iterative tolerances.

The bundled 64 by 64 cavity example is now driven by the native
`examples/meshes/cavity.mesh` file. Its converged Ghia comparison is recorded
in `docs/validation.md`; the old boundary-cell input convention is no longer
part of the BabelSim example path.

## What is retained

The following TaihoCFD algorithms were retained as numerical references and
migrated with equivalence tests rather than rewritten from memory:

- the segregated SIMPLE state transition and relaxation convention;
- first-order upwind convection and the orthogonal diffusion limit;
- Rhie-Chow momentum interpolation;
- pressure-outlet and closed-domain pressure-reference treatment;
- face-flux, pressure, and cell-velocity correction;
- continuity and relative field-change metrics;
- BiCGSTAB/ILUT and PCG/incomplete-Cholesky solver behaviour;
- compressed cell indexing, integer topology, and contiguous storage.

TaihoCFD's MPI owned/ghost semantics are retained in the distributed layer.
They are implemented as metadata and indexed halo maps over the same local
`Mesh`/`Field` storage, rather than as a second CFD data model.

## What is not migrated as framework structure

TaihoCFD's `Mesh` owns `u`, `v`, `p`, time history, face velocities, geometry,
and boundary state. This couples spatial topology to Navier-Stokes physics.
Its boundary patches contain boundary *cells*, while a general finite-volume
mesh needs patches of boundary faces. `Equation` is a fixed two-dimensional
five-point stencil, and operators directly address east/west/north/south Eigen
matrices. These structures are replaced rather than copied.

Legacy integer boundary codes remain an input-adapter concern only. Compatibility
overloads that merely reconstruct global MPI state, and the virtual
`FlowSolver` factory used for a single algorithm, are not part of the new core.

## Concrete TaihoCFD to BabelSim map

| TaihoCFD component | Decision | BabelSim component |
|---|---|---|
| `mesh/Mesh` geometry plus `u/v/p` arrays | split; keep indexing ideas, remove physics ownership | `Mesh` geometry/topology and independent `Field<T>` |
| boundary cells, integer `bctype/zoneid` | delete from the runtime core | native `Mesh` file patches and field boundary conditions |
| fixed E/W/N/S operator loops | retain verified formulas in the orthogonal limit; generalize by faces | free functions in `operators` with `Method` enums |
| 2D five-point `Equation` | replace | face-addressed scalar/vector LDU `Equation<T>` |
| Eigen solver wrappers | retain solver/preconditioner choices and residual checks | CFD-independent `PreparedLinearSolver` |
| `SimpleSolver`, momentum, Rhie-Chow, pressure correction | retain numerical state transition; split into named methods | `SimpleSolver` over generic fields/equations/operators |
| virtual `FlowSolver` plus one concrete implementation | delete | direct construction/composition |
| MPI x-decomposition and global reconstruction overloads | retain the owner/ghost idea, replace global reconstruction with owned output | `ParallelContext`, `decompose`, `HaloExchange`, global reductions, and rank-owned output |

No TaihoCFD source file is copied verbatim. TaihoCFD remains a numerical
reference only; its old boundary-cell input format is not a BabelSim runtime
dependency.

## Flat target modules

```text
Mesh        topology and geometry only
Field<T>    contiguous values on cells, faces, or vertices; field BC metadata
Operator    mathematical action on fields or equations
Method      small enums/configuration selecting a discretization
Equation<T> LDU coefficients and source accumulated from operators
Assembly    Equation -> sparse algebraic matrix/vector
Solver      Ax=b only; no CFD meaning
Physics     direct composition, initially incompressible SIMPLE
```

There are no manager, registry, provider, or deep factory layers. Hot paths use
free functions, small value types, and direct indexed arrays; runtime virtual
dispatch is not used by the mesh, fields, operators, assembly, or linear solver.

## Unified 3D mesh

The structured mesh stores logical dimensions `(nx, ny, nz)` but exposes an
explicit finite-volume topology:

- contiguous vertices and cell/face geometry;
- `face_owner`, `face_neighbour`, and four vertex indices per quadrilateral face;
- six face and neighbour indices per hexahedral cell;
- boundary patches containing face indices;
- cell centre and volume;
- face centre, area vector, area, orthogonal coefficient, non-orthogonal
  correction vector, and skewness vector.

The area vector of an internal face points from owner to neighbour; a boundary
face points out of its owner. Diffusion uses

```text
Sf = (Sf dot d / |d|^2) d + k_nonorth
```

so the first part is implicit and `k_nonorth dot grad(phi)` is an explicit
correction. Face skewness is the offset between the geometric face centre and
the owner-neighbour line intersection with the face plane.

A two-dimensional mesh is an extrusion with `nz=1`. Front and back patches use
symmetry conditions; volumes and in-plane areas retain a finite thickness.
No solver or operator has a two-dimensional branch.

## Native mesh files

Executable examples and mesh-file regression tests use the versioned BabelSim
text format, not a solver-specific legacy input format:

```text
BABELSIM_MESH 1
dimensions 62 70 1
geometry cartesian
bounds 0 0 0 10 1 1
patch xmin inlet inlet
patch xmax outlet outlet
patch ymin lower_wall wall
patch ymax upper_wall wall
patch zmin front symmetry
patch zmax back symmetry
end
```

`geometry vertices` may be used instead of `geometry cartesian`, followed by
the `(nx+1)(ny+1)(nz+1)` vertex triples in logical order. The file describes
geometry and patch roles only; physical field values remain Field boundary
conditions, so the same mesh can be reused by different PDEs.

MPI output is intentionally a separate concern: `writeOwnedCsv` writes one
CSV per rank with only owned cells. `tools/merge_parallel_csv.py` validates
global-id uniqueness and emits a legacy VTK point cloud plus a Tecplot ASCII
point zone. This keeps the runtime writer lightweight while providing a
standard visualization interchange format without serializing ghost cells.

## Field and boundary model

`Field<T>` supports `double` and `Vec3` without naming physical variables. Its
values are contiguous and its location is cell, face, or vertex. Cell fields
own one boundary-condition record per mesh patch. Supported mathematical
conditions are fixed value/Dirichlet, fixed normal gradient/Neumann,
zero-gradient, inlet-outlet, and symmetry/mirror. Wall, inlet, outlet, and
symmetry are patch roles; each physical field still chooses its own condition.

For vector symmetry, the boundary face value removes the normal component and
keeps tangential components. Scalar symmetry is zero normal gradient.
`InletOutlet` switches on the signed face flux in convection and boundary-face
evaluation; diffusion without a supplied flux uses its zero-gradient/outflow
branch. A known velocity or scalar inlet should use `FixedValue`.

## Operator, method, equation, and assembly interfaces

Operators are functions over fields or coefficient accumulators. Methods are
orthogonal configuration values, not subclasses:

```cpp
interpolate(cell, face, InterpolationMethod::Linear);
gradient(phi, grad_phi, GradientMethod::GreenGauss);
reconstruct(phi, grad_phi, face); // includes geometric skewness correction
flux(U, face_flux, InterpolationMethod::Linear);
divergence(U, div_u, InterpolationMethod::Linear);

addConvection(equation, face_flux, phi, ConvectionMethod::Upwind);
addDiffusion(equation, gamma, phi,
             GradientMethod::GreenGauss, DiffusionMethod::Corrected);
addDiffusion(equation, face_gamma, phi, // variable face coefficient
             GradientMethod::LeastSquares, DiffusionMethod::Corrected);
addTimeDerivative(equation, old, dt, rho, TimeMethod::Euler);

LinearSystem system = assemble(equation);
solve(system.A, system.b, x, solver_config);

PreparedLinearSolver prepared(solver_config);
prepared.compute(system.A);       // analyze pattern and factorize once
prepared.solve(rhs0, x0);
prepared.solve(rhs1, x1);         // reuse for another right-hand side

SparseAssembly cached(mesh);      // construct mesh-dependent CSR pattern once
cached.update(equation);          // update only coefficient values
prepared.factorize(cached.matrix());
```

`Equation<T>` has one diagonal per cell, lower/upper coupling per mesh face, and
a contiguous source. It works for scalar and segregated vector equations.
Sparse assembly knows topology and coefficient signs but no PDE or CFD
semantics. `SparseAssembly` precomputes the compressed nonzero structure and
the direct value-array positions of diagonal/upper/lower coefficients; SIMPLE
therefore performs no triplet allocation, sorting, or sparse insertion during
an iteration. The linear solver sees only `A`, `b`, and `x`.

## MPI data ownership and operator execution

`ParallelContext` is a small MPI-3 value object. A decomposed structured mesh
keeps all local cells needed by its stencil, marks a contiguous x-range as
owned, and stores two ghost layers on each internal partition boundary. The
`owned_cells`, `cell_owned_indices`, and `cell_global_ids` maps make the local
algebra and output compact without changing field indexing. Physical boundary
patches remain field boundary conditions; partition faces are explicit
`PatchKind::Processor` faces.

MPI is not hidden inside every scalar loop. Instead, each operator stage has an
explicit synchronization point: the producing cell/face field is exchanged
before a consumer reads a ghost value. This keeps generic operators reusable in
serial and distributed runs and avoids collectives in inner loops. Current
`HaloExchange` covers cell-centred scalar/vector/tensor fields, interface face
fields, and raw cell vectors used by distributed matrix-vector products. Global
dot products, norms, solver status, SIMPLE residuals, continuity, and field
change metrics use `MPI_Allreduce`; rank-owned CSV files provide deterministic
parallel output for post-processing. Vertex-field halo exchange is deliberately
not claimed yet because no current distributed operator consumes vertex values.

The intended call sequence is explicit and the same for every physical model:

```cpp
HaloExchange halo(mesh, parallel);
halo.exchange(phi);                 // cell scalar/vector/tensor
gradient(phi, grad_phi, method);    // local cells plus synchronized ghosts
halo.exchange(face_flux);            // when the next stage consumes face data
addConvection(equation, face_flux, phi, convection_method);
addDiffusion(equation, gamma, phi, gradient_method, diffusion_method);
```

For a transient term, the caller exchanges both `previous` and `older` before
`addTimeDerivative`; the operator itself remains a local algebraic operation.

The distributed Krylov solver assembles the owned block with the same cached
`SparseAssembly` used in serial. Its matvec adds remote face couplings after a
cell halo exchange, while incomplete-Cholesky/ILUT remains a local block
preconditioner. Thus momentum, pressure correction, gradients, fluxes,
diffusion, and continuity all operate on distributed fields; MPI is not a
SIMPLE-only wrapper.

SIMPLE-specific momentum interpolation and pressure correction remain named
algorithm components inside `SimpleSolver`. They consume the same Mesh, Field,
Equation, generic operators, assembly, and solver APIs; neither sparse assembly
nor the linear solver contains pressure, velocity, or CFD semantics.

## Migration status and remaining gates

1. **Complete:** core geometry, fields, patch conditions, LDU equations,
   assembly, solvers, skewed 3D geometry, and `nz=1` degeneration.
2. **Complete:** generic interpolation/reconstruction, scalar and vector
   gradients, divergence, flux, convection, corrected diffusion, and
   Euler/BDF2 terms. Affine linear and cross-quadratic manufactured tests cover
   non-orthogonality and skewness.
3. **Complete for single-node execution:** native BabelSim mesh-file reader and
   three-component SIMPLE/Rhie-Chow migration. The old boundary-cell
   representation is not retained in the core.
4. **Complete:** converged Ghia Re=100 and analytic Poiseuille gates. Reaching an
   iteration limit is never accepted as a result.
5. **Complete as a regression gate:** all three face orientations, corrected
   momentum diffusion, pressure/velocity/face-flux correction, a warped cavity,
   and a non-orthogonal manufactured Laplacian. A refined grid-convergence
   study remains pending.
6. **Partly complete:** one-cell-thick and true 3D cavities use the same solver;
   mass conservation and symmetry pass. The 3D case is still a coarse
   regression and needs an external quantitative reference before an accuracy
   claim.
7. **First profiling pass complete:** sparse ordering/preconditioner setup is
   the measured hotspot. Reuse across velocity components and iterations,
   precomputed sparse topology/index maps, and preallocated right-hand sides
   are implemented.
8. **Complete for the structured MPI baseline:** TaihoCFD-informed x-domain
   decomposition, two-layer cell and interface-face halo exchange, distributed
   Krylov matvec/global reductions, distributed SIMPLE metrics, and owned-cell
   rank output. One-, two-, and four-rank cavity tests plus one-/two-rank
   native Poiseuille equivalence are regression gates.
9. **Remaining:** general (non-structured) partitioning, vertex-field halos,
   scalable distributed preconditioners, refined 3D benchmark validation, and
   scaling/memory-bandwidth studies.
