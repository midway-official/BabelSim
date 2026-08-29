# Validation status

All values below are from converged runs. Small 3D and warped-grid cases are
regression tests for algorithms and geometry; only the Ghia and analytic
Poiseuille comparisons are quantitative accuracy claims.

## Automated regression suite

`make test` currently covers:

- Cartesian `nz=1` degeneration and a skewed 3D hexahedral geometry;
- scalar, vector, and symmetry/mirror boundary conditions;
- TaihoCFD boundary-cell to boundary-face import;
- Green-Gauss and least-squares scalar/vector gradients;
- scalar and vector non-orthogonal diffusion plus skewness reconstruction;
- a cross-quadratic manufactured solution that distinguishes corrected from
  orthogonal-only diffusion on an affine skew mesh;
- divergence, flux, Laplacian, upwind/central convection, Euler/BDF2 terms;
- LDU assembly, constant/face-variable diffusivity, cached sparse topology,
  CG/incomplete-Cholesky solution, and prepared-solver reuse;
- imported-channel, 2D cavity, true 3D cavity, and warped-mesh cavity SIMPLE.
- MPI ownership maps, two-cell-layer scalar/vector/tensor and interface-face
  halos, distributed sparse matvec/Krylov solves, global reductions, and
  owned-cell output.

Representative regression results:

| Case | Cells | SIMPLE iterations | Relative mass imbalance | Other gate |
|---|---:|---:|---:|---|
| imported small channel | 9 | 49 | `1.59e-16` | `max(abs(w)) < 1e-13` |
| 2D cavity | 12 x 12 x 1 | 137 | `1.61e-16` | centre `u=-0.149236` |
| 3D cavity | 6 x 6 x 6 | 57 | `1.59e-15` | centre `u=-0.114134`, `max(abs(w))=0.0192033` |
| warped cavity | 10 x 10 x 1 | 151 | about `6.0e-16` | `max(abs(k_nonorth))/abs(Sf)=0.412398` |

The 3D case additionally checks centre-plane velocity symmetry. The nonzero
spanwise velocity shows that it is a true wall-bounded 3D solution, not a
duplicated 2D field.

## Ghia lid-driven cavity

Run:

```bash
make validate-cavity
```

Configuration: steady incompressible SIMPLE, `Re=100`, 64 x 64 x 1 cells,
first-order upwind, orthogonal diffusion, velocity tolerance `1e-6`. The run
converges in 2355 iterations with relative mass imbalance `1.08e-14`.

Comparison with Ghia centreline samples:

| Metric | Value |
|---|---:|
| horizontal-velocity centre interpolation | `-0.19697720` |
| vertical-velocity centre interpolation | `0.05108802` |
| horizontal-velocity L-infinity error | `0.011111441` |
| horizontal-velocity L2 error | `0.005884307` |
| vertical-velocity L-infinity error | `0.007404026` |
| vertical-velocity L2 error | `0.003986802` |

The iteration count and horizontal-velocity metrics reproduce the converged
TaihoCFD reference at the configured discretization and tolerances.

## Poiseuille flow

Run:

```bash
make validate-poiseuille
```

The TaihoCFD mesh adapter solves 4340 interior cells. The run converges in 897
SIMPLE iterations (TaihoCFD reference: 896), with relative mass imbalance
`3.03e-15`. At the outlet, comparison with the fully developed parabola gives:

| Metric | Value |
|---|---:|
| maximum velocity | `1.49809557` |
| L-infinity error | `0.000371012` |
| L2 error | `0.000286418` |

The one-iteration difference is within the configured iterative tolerances and
the boundary representation changed from TaihoCFD ghost/boundary cells to
finite-volume boundary faces.

## MPI domain decomposition and equivalence

Run the structured MPI regression gates with:

```bash
make test-mpi
make test-mpi-poiseuille
```

`test-mpi` exercises the complete distributed operator path on one, two, and
four ranks for a 12 x 12 x 1 cavity. The 2-rank domain test additionally checks
global cell coverage, two ghost layers, scalar/vector/face halos, Green-Gauss
and least-squares operators at a partition interface, corrected diffusion,
owned sparse assembly, and a distributed affine diffusion solve. The cavity
results are:

| Ranks | SIMPLE iterations | Relative mass imbalance | centre `u` |
|---:|---:|---:|---:|
| 1 | 137 | `1.14e-14` | `-0.149236` |
| 2 | 137 | `1.49e-14` | `-0.149236` |
| 4 | 137 | `1.94e-14` | `-0.149236` |

The same target runs `parallel_cavity_3d_test` on two ranks for a 6 x 6 x 6
wall-bounded cavity. It converges in 57 iterations with relative mass
imbalance `3.09e-15`, centre `u=-0.114134`, and `max(abs(w))=0.0192033`; the
centre-plane symmetry gate passes.

`test-mpi-poiseuille` solves the imported 4340-cell channel on one and two
ranks and compares rank-owned CSV files by global cell id. Both runs converge
in 897 SIMPLE iterations; the maximum component-wise field difference is
`2.63e-7` with the comparison gate `atol=rtol=5e-6`. MPI reduction order and
distributed Krylov trajectories are therefore allowed to introduce bounded
floating-point differences; bitwise identity is neither required nor claimed.

Parallel output is one file per rank (`<prefix>_<rank>.csv`) containing only
owned cells and a global id. The files can be merged or compared without
duplicating ghost layers.

## Memory safety and profiling

The complete serial regression suite passes with AddressSanitizer and
UndefinedBehaviorSanitizer, including `UBSAN_OPTIONS=halt_on_error=1`. The MPI
sanitizer suite also passes with leak detection disabled; OpenMPI/UBSan reports
runtime-owned allocations from `MPI_Init`/`MPI_Finalize` when leak detection is
enabled, so those library leaks are not attributed to BabelSim.

GNU gprof identifies sparse ordering, ILUT/incomplete-Cholesky factorization,
and Krylov solves as the dominant work. Reusing one prepared momentum system
for all three velocity components, precomputing the compressed matrix topology
and coefficient positions, preallocating iteration workspaces, and avoiding
unused orthogonal-path gradient fields reduced the measured 64 x 64 cavity
wall time from 29.70 s to 19.34 s (34.9%) on the same host. The Poiseuille run
fell from 25.96 s to 15.46 s (40.4%).
Iteration counts and validation errors are unchanged, and the solver-only peak
resident set remains about 10-11 MB. The remaining optimization target is
numerical preconditioner work, not virtual dispatch or pointer-heavy mesh
traversal.

## Gates still pending

- refined 3D cavity comparison against an external published dataset;
- grid-convergence study for the warped/non-orthogonal flow case;
- transient temporal-order checks beyond operator-level Euler/BDF2 tests;
- general unstructured partitioning, vertex-field halos, and scalable global
  preconditioners;
- scaling and memory-bandwidth profiling on representative 3D grids.

The MPI baseline is therefore a validated structured-domain result, not a
claim that all future mesh locations or arbitrary partitioners are complete.
