# Procedural numerical DSL migration

The target is a C++ embedded procedural DSL. Physics authors use ordinary
functions, local fields and explicit time/nonlinear/correction loops. Input field
and checked dictionary APIs remain available. `equ` assembles an integrated
linear system Ax=b immediately; `math` computes known fields. Solving must not
reassemble or silently iterate a physical algorithm.

Runtime output policy: solver programs own progress, residual, performance and
completion messages. Runtime APIs return data only; the launcher does not infer
physical convergence or automatically finalize/write a successful case. The application supplies any error-printing callback; the runtime itself has no
printing implementation.

## Required migration and acceptance

- [ ] Opaque scalar/vector matrix with immediate assembly, explicit histories,
      source/reaction signs, algebra, relaxation, reference and diagnostics.
- [ ] Direct field arithmetic, interpolation/differential operations, global
      reductions and well-defined assignment/alias semantics.
- [ ] Explicit time setting, history snapshots and writes; no physics loop()
      lifecycle controllers.
- [ ] Heat and transport expressed as procedural solver functions.
- [ ] Steady/transient SIMPLE expressed with explicit correction loops and
      ordinary functions; consistent pressure flux and momentum response.
- [ ] All RANS transport modules use the new DSL and explicit histories.
- [ ] SDK examples, documentation and architecture rules describe the new API.
- [ ] Numerical tests: assembly signs/volume scaling/frozen inputs, history,
      boundary/alias behavior, matrix operations and MPI ownership.
- [ ] Existing numerical and workflow regressions, SIMPLE and RANS acceptance.

Existing uncommitted numerical changes are preserved. This checklist records
remaining scope, not a claim that adding a compatibility facade completes it.

## Layer boundaries

- Physics: numerical algorithms, explicit iteration, stopping criteria, logs.
- Public DSL: typed fields, opaque integrated matrices, known-field arithmetic.
- Discretization: geometry, boundary integration, assembly and flux consistency.
- Backend: sparse algebra, linear solves, halo exchanges and global reductions.
- Case/IO: checked file configuration, field ownership, explicit result writes.

Physics must not include runtime implementation, MPI, Eigen, CSR or LDU. The
runtime must not identify physical algorithms or choose their stopping policy.

## Updated user requirements (2026-09-16)

- Explicitly loading schemes and enabling/advancing time are runtime services.
  BDF2 startup and endpoint arithmetic must not appear in solver programs.
- Steady SIMPLE has an explicit outer loop in its own main.cpp; transient SIMPLE
  has explicit time/outer loops in its own main.cpp. Optional pressure correction
  subloops are allowed. No shared SIMPLE implementation between the two.
- Each turbulence model independently owns transport, histories and closure code.
  A pure interaction interface/factory is allowed; a shared numerical base is not.
- Physics/application code owns every print operation and physical stopping test.

Progress evidence so far (not full acceptance): procedural matrix and eager math
unit tests passed; the standalone steady Poiseuille entry converged in 1155 outer
iterations (mass relative error 2.63e-15). The RANS verification initially found a
test fixture that discarded calculated boundary traces when assigning new eager
results; the fixture has been corrected to preserve its intended derived fields.
The RANS suite must be rerun. Tests under tests/support retain a test-only SIMPLE
reference, not a production algorithm and not evidence that production is tested.
