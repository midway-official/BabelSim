# BabelSim

BabelSim is a compact C++17 finite-volume framework for CFD, PDE, and
multiphysics solvers. The framework separates mesh topology and geometry,
contiguous fields, mathematical operators, discretization methods, equation
assembly, algebraic solvers, and physics algorithms.

The current development baseline is a unified three-dimensional structured
mesh: a two-dimensional case is an extruded mesh with one cell in the
z-direction and uses the same geometry and operator code. Orthogonal meshes
are a special case of the non-orthogonal geometry representation.

Implemented:

- explicit cell/face/vertex topology with owner-neighbour and boundary patches;
- three-dimensional hexahedral geometry, face area vectors, non-orthogonal and
  skewness vectors;
- contiguous scalar/vector/tensor fields with field-specific fixed-value,
  fixed-gradient, zero-gradient, inlet-outlet, and symmetry/mirror conditions;
- Green-Gauss and least-squares gradients, interpolation, flux, divergence,
  skew-corrected reconstruction, convection, corrected diffusion, and
  Euler/BDF2 time terms;
- LDU finite-volume equations, constant or face-centred variable diffusion,
  precomputed sparse assembly structure, and CFD-independent Eigen linear
  solvers with reusable sparsity analysis and preconditioners;
- a versioned native mesh-file reader and a three-component SIMPLE solver
  with Rhie-Chow momentum interpolation, pressure correction, cell/face-flux
  correction, and continuity monitoring.

TaihoCFD remains the numerical reference rather than a copied framework. See
[`docs/architecture.md`](docs/architecture.md) for the audited mapping and
[`docs/validation.md`](docs/validation.md) for quantitative and regression
evidence.

## Build and test

Dependencies are a C++17 compiler, Eigen 3, OpenMPI (or another MPI-3
implementation), and GNU Make. The default Makefile compiler is `mpic++` so
serial tests and MPI executables use one ABI.

```bash
make -j
make test
make test-mpi
make test-mpi-poiseuille
make validate-cavity
make validate-poiseuille
# after the MPI target, merge its rank files for ParaView/Tecplot
make postprocess-mpi-poiseuille
```

`validate-poiseuille` uses the native `examples/meshes/poiseuille.mesh` file
by default; set `POISEUILLE_MESH=/path/to/mesh` to override it.

Both examples read a versioned BabelSim mesh file (`examples/meshes/*.mesh`);
the files contain dimensions, Cartesian or explicit vertex geometry, and six
boundary-patch records.  MPI field output is one owned-cell CSV per rank.  Use
`tools/merge_parallel_csv.py <directory> --prefix <prefix> --ranks <N>` to
write a merged legacy VTK point-cloud (`.vtk`) and Tecplot ASCII point zone
(`.dat`) that can be opened directly by ParaView or Tecplot.

For example, after `make test-mpi-poiseuille`, the merged files are written to
`build/postprocess/poiseuille.vtk` and `build/postprocess/poiseuille.dat`.

The validated baseline includes serial and MPI structured hexahedral finite volumes. Quantitative
2D cavity and Poiseuille gates pass; the 3D and non-orthogonal cases are
converged regression gates, not yet literature-accuracy claims. MPI currently
decomposes structured meshes in x, exchanges two cell layers and interface
face fields, reduces Krylov/SIMPLE metrics globally, and writes owned-cell
rank files. Refined 3D validation, turbulence models, and general unstructured
topology remain later stages.
