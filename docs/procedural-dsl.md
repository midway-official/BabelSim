# Procedural numerical DSL

The solver is an ordinary C++ function. It reads checked configuration and fields,
then explicitly assembles, solves, updates, checks convergence and writes results.
`src/physics/heat/main.cpp` is the minimal working example. Steady and transient
SIMPLE each live entirely in their own `main.cpp`; neither calls a shared SIMPLE
algorithm. SA, k-omega and k-epsilon each own their transport assembly and state.
Their common `RANS/api.h` is a pure interaction contract, not a numerical base class.

## Responsibilities

| Layer | Responsibility |
| --- | --- |
| Solver/model program | Physical equations, algorithm order, all nonlinear/correction loops, stopping criteria, logging |
| Public DSL | Fields, eager arithmetic, opaque integrated matrices, numerical histories |
| Discretization | Schemes, boundary integration, implicit/explicit contributions and consistent face fluxes |
| Backend | Sparse algebra, linear convergence, halo exchange, global reductions |
| Runtime and IO | MPI lifetime, explicitly loaded methods, explicitly enabled time service, file loading and explicit writes |

The runtime does not print, decide physical convergence, run a SIMPLE iteration,
or automatically write results based on the solver's return code. The application
may supply its own error-printing callback. `primaryProcess()` lets solver code
choose one-process logging without importing MPI.

## Loading and time

```cpp
auto& T = problem.scalarField("T");
const double k = problem.physics().nonnegative("conductivity");
loadMethods(problem);
const auto linear = readLinearControl(problem, T);
const int interval = readWriteInterval(problem);
auto time = enableTime(problem);
auto old = math::history(T);
auto A = equ::createEquation(T);
problem.validate();

while (time.value() < time.end()) {
    advance(time);
    math::saveOld(old, T, time.dt());
    equ::reset(A);
    equ::ddt(A, 1.0, old);
    equ::laplacian(A, k, -1.0);
    const auto result = equ::solve(A, T, linear);
    if (!result.converged()) return 2;
    if (time.step() % interval == 0 || time.finished()) write(problem, time);
}
```

`advance()` only updates time metadata. It handles floating-point endpoints and a
short final step; it never saves fields or writes output. `saveOld()` is the only
history-advancing operation. Call it once per physical step, outside SIMPLE or
other nonlinear loops. `equ::ddt()` selects the loaded time method, starts BDF2
with Euler, and uses the actual previous/current dt for variable-step BDF2.
A plain `math::copy(T)` plus `equ::ddt(A, capacity, oldField, dt)` remains available
for an explicitly selected Euler discretization.

Manual time control is also supported: `readTimeControl(problem)` returns
`{start,end,dt}`, `setTime(problem,t)` sets evaluation metadata, and
`write(problem,t,step)` explicitly writes a snapshot. Write calls update the time
series and latest-written snapshot; they do not declare the algorithm converged.
No write occurs implicitly when a solver returns zero.

## Fields and direct mathematics

```cpp
auto gradP = math::grad(p);
auto Uf = math::interpolate(U);
auto phi = math::flux(Uf);
auto stress = mu * (math::transpose(math::grad(U))
    - (2.0 / 3.0) * math::isotropic(math::trace(math::grad(U))));
U -= rAU * gradP;
double mass = math::integral(rho);
```

Public differential/interpolation functions return independent, evaluated fields;
`auto gradP = math::grad(p)` does not retain a deferred reference to p. Arithmetic
creates computed boundary traces. Assignment copies values while preserving the
destination's name, mesh and boundary constraints. Use `auto` for derived fields,
or explicitly mark a reusable derived destination with `useCalculatedBoundary()`.
Layouts must match; cell and face fields are not implicitly interchangeable.
Whole-field stencils synchronize through the backend. Reductions use owned
entities globally, excluding duplicate ghost values. `normL2` is an unweighted
Euclidean norm; `integral` applies cell volume or face area weights.

## Immediate equation assembly

Names distinguish construction, assembly, initialization and solution:

```cpp
auto pPrime = math::createHomogeneousField(p); // zero field, homogeneous BCs
auto pressureCorrectionEquation = equ::createEquation(pPrime); // empty system
pPrime.fill(0.0);                              // initialize the unknown values
equ::reset(pressureCorrectionEquation);        // remove matrix and RHS terms
equ::laplacian(pressureCorrectionEquation, rAU, -1.0);
equ::source(pressureCorrectionEquation, -math::div(phiH));
equ::reference(pressureCorrectionEquation, 0.0);
const auto pressureSolve = equ::solve(pressureCorrectionEquation, pPrime, linear);
```

Creating an equation does not infer any physics or assemble a PDE from the field.
Its terms define the equation. Resetting retains the unknown binding and storage,
and does not reset field values. Solving writes the solution to the named unknown
and returns linear-solve status; it does not reassemble or test physical convergence.
The older `matrix`, `clear` and `correction` spellings remain compatibility APIs.

In SIMPLE, the momentum equation is assembled once per outer iteration. Its
residual is measured immediately, before relaxation, row scaling or the velocity
predictor solve. The end-of-iteration stopping test combines this saved residual
with corrected continuity, field changes and turbulence diagnostics. Corrections
made during this iteration enter the next iteration's normal assembly. This is a
beginning-of-iteration nonlinear residual, not the linear solver's final residual
or a freshly reassembled end-of-iteration residual. Transient SIMPLE includes the
fixed physical-step history in this assembly and compares successive inner
iterates, not successive physical time levels.

Pressure correction fields preserve the homogeneous counterparts of the physical
boundary constraints:

```cpp
auto pPrime = math::createHomogeneousField(p);
auto pressureCorrectionEquation = equ::createEquation(pPrime);
```

After assembly, `equ::reference(pressureCorrectionEquation, 0.0)` adds a reference only when the assembled
matrix has a constant null mode. This is intended for pressure systems; it does
not diagnose arbitrary singularities or disconnected components.

## Turbulence interaction

```cpp
auto turbulence = rans::load(problem, U, phi);
const auto& muEff = turbulence.effectiveViscosity();
// Once per physical time step, before any inner corrections:
turbulence.saveOld(time.dt());
// At the explicitly chosen point in the SIMPLE iteration:
const auto result = turbulence.solveTransport();
```

The coupling object reads the configured model and molecular properties and
provides effective dynamic viscosity. Each concrete model independently implements
its transport equations and closure. Each result retains named per-equation linear status, normalized pre-relaxation
residual and bounded field change. Only dimensionless maxima are aggregated;
absolute residuals of different physical quantities are never added. The caller
owns logging and the physical convergence decision. Laminar configuration supplies molecular viscosity without
transport equations. Steady solvers do not save transient histories.

## Matrix operations

`equ::Matrix<T>` represents an integrated cell system **Ax=b**; its storage stays
opaque. Every assembly statement freezes its inputs. `solve()` does not assemble,
update coefficients, run outer iterations or select a physical stopping policy.

| Operation | Contribution |
| --- | --- |
| `ddt(A,capacity,history)` | LHS time derivative and RHS historical terms |
| `div(A,phi,scale)` | LHS `scale*div(phi*x)` |
| `laplacian(A,k,scale)` | LHS `scale*div(k*grad(x))`; normally use `-1` for diffusion |
| `reaction(A,c)` | LHS `c*x` |
| `source(A,s)` | RHS volumetric source, integrated by the framework |
| `addRhs(A,b)` | Already integrated algebraic RHS; no extra volume factor |
| `addDiagonal(A,d)` | Already integrated diagonal |

`clear` retains the matrix layout. `copy`, `scale`, `add`, `diagonal`, `rhs`,
`apply`, `residual` and `relativeResidual` operate on assembled coefficients.
Residual is `b-A*x`. `response(A)` returns `V/aP` for that exact matrix.
`relax(A,xOld,alpha)` increases the diagonal by `1/alpha` and compensates the RHS.
The SIMPLE programs explicitly scale the relaxed rows to preserve their established
response convention. A separated vector system is not a fully coupled block matrix.

`faceFlux(P,pPrime)` returns the signed diffusion contribution using the boundary
conditions, face coefficients and deferred-correction gradients frozen during
assembly. This keeps pressure flux correction consistent with the solved system.
Nonorthogonal correction iterations remain explicit in each SIMPLE program.

`readLinearControl` returns a value that can be changed and passed to a particular
`equ::solve(A,x,options)`. Backend caches are selected by these settings; the
argument is not ignored. Default-option overloads remain available.

## Current validation boundary

The migration is in progress. See `design/procedural-dsl-migration.md` for the
remaining acceptance scope. Existing expression APIs are retained for regression
compatibility; new solver examples and production physics use `equ` and eager
`math`. Restart serialization of multiple BDF2 history levels and full block
coupled matrices are not implemented by this migration yet.
