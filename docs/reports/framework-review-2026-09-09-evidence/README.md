# Review evidence

Baseline: `e112e597172b7fca5956504c4905e19add13d252`, 2026-09-09.
The production sources and checked-in tests were not edited during this review.
`source-sha256.txt` identifies the reviewed `src` and `include` files.

`serial.log`, `mpi.log`, `validation.log` preserve the completed repository
regressions and converged laminar validations. `baseline.txt` lists their commands.
`first-pass-smoke` preserves the earlier steady/Euler × none/SA/kOmega/kEpsilon
integration runs at 1 and 2 ranks. Those smoke runs are not turbulence validation.

Run `python3 docs/reports/framework-review-2026-09-09-evidence/reproduce.py`
from a checkout built with `make -j4`. This requires mpic++, MPI and Eigen3 in
the paths used by the Makefile. It writes only to a new temporary directory.
The script preserves the original probe code apart from relocating two case paths.
`-march=native` matches the reviewed library's Eigen configuration.

Observations to inspect (these programs intentionally print discrepancies and
generally exit zero; they are diagnostic probes, not passing acceptance tests):

| Probe | Observation |
|---|---|
| scale | CG/BiCGSTAB report converged, x=0, relativeResidual=1 for A=b=1e-20 |
| rans | clipping yields dTurb=0 while post-update k equation residual is -0.001 |
| sa-boundary | effective viscosity at zero-nuTilda wall is 0.00836425 instead of 0.001 |
| gg | default 2-layer halo, warped mesh/nonlinear field: serial/MPI max error 0.000246513 |
| gg3 | same Green–Gauss case with 3 layers: max error 0 |
| ls | same case with LeastSquares and 2 layers: max error 8.88178e-16 |
| inletoutlet | fixed-value inflow solution 0.6; inletOutlet inflow solution 0.333333 |
| transient | Euler/BDF2 discrete diffusion ODE errors below 9e-12 |
| post | accepts old results after mesh vertices are translated by +10 in x |

The archived `post-case` intentionally contains mismatched mesh/results. Its CSV
coordinates retain the original geometry; its current mesh has been translated.
It exists solely to reproduce missing geometry provenance validation.
The RANS full-stress finding is established analytically in the report.
