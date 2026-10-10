#!/usr/bin/env python3
"""Build a SIMPLE or PISO flow case using a generated triangular-prism mesh."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil


U_FIELD = """field U
{
 type vector
 location cell
 internal uniform (0 0 0)
 boundary
 {
  inlet { type fixedValue value (1 0 0) }
  outlet { type zeroGradient }
  lower_wall { type fixedValue value (0 0 0) }
  upper_wall { type fixedValue value (0 0 0) }
  front { type symmetry }
  back { type symmetry }
 }
}
"""

P_FIELD = """field p
{
 type scalar
 location cell
 internal uniform (0)
 boundary
 {
  inlet { type zeroGradient }
  outlet { type fixedValue value (0) }
  lower_wall { type zeroGradient }
  upper_wall { type zeroGradient }
  front { type symmetry }
  back { type symmetry }
 }
}
"""

COMMON_METHODS = """equation.momentum.interpolation linear
equation.momentum.gradient greenGauss
equation.momentum.convection upwind
equation.momentum.diffusion orthogonal
equation.pressureCorrection.interpolation linear
equation.pressureCorrection.gradient greenGauss
equation.pressureCorrection.convection upwind
equation.pressureCorrection.diffusion orthogonal
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-case", type=Path, required=True,
                        help="generated heat case whose mesh will be reused")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--solver", choices=("simple", "piso"), required=True)
    parser.add_argument("--max-iterations", type=int, default=100,
                        help="SIMPLE outer cap; PISO performs one step")
    args = parser.parse_args()
    if args.max_iterations < 1:
        parser.error("max-iterations must be positive")

    source_mesh = args.mesh_case.resolve() / "mesh.mesh"
    case = args.output.resolve()
    case.mkdir(parents=True, exist_ok=True)
    (case / "mesh").mkdir(exist_ok=True)
    mesh_path = case / "mesh/flow.mesh"
    rename = {
        "patch hot wall": "patch inlet inlet",
        "patch cold wall": "patch outlet outlet",
        "patch lower symmetry": "patch lower_wall wall",
        "patch upper symmetry": "patch upper_wall wall",
    }
    with source_mesh.open(encoding="ascii") as source, mesh_path.open("w", encoding="ascii") as target:
        for line in source:
            stripped = line.strip()
            replacement = rename.get(stripped)
            target.write(replacement + "\n" if replacement is not None else line)

    fields = case / "fields/initial"
    fields.mkdir(parents=True, exist_ok=True)
    (fields / "U.field").write_text(U_FIELD, encoding="utf-8")
    (fields / "p.field").write_text(P_FIELD, encoding="utf-8")
    physics = case / "physics"
    physics.mkdir(exist_ok=True)
    (physics / "flow.bs").write_text(
        "density 1\ndynamicViscosity 0.01\nturbulenceModel none\n", encoding="utf-8")
    numerics = case / "numerics"
    numerics.mkdir(exist_ok=True)
    if args.solver == "simple":
        methods = "time steady\n" + COMMON_METHODS
        solution = (
            f"maxIterations {args.max_iterations}\nvelocityRelaxation 0.5\npressureRelaxation 0.3\n"
            "continuityTolerance 1e-7\nvelocityTolerance 1e-6\n"
            "momentumTolerance 1e-6\npressureCorrectionTolerance 1e-6\n"
            "equation.momentum.kspType bcgs\nequation.momentum.pcType bjacobi\n"
            "equation.momentum.absoluteTolerance 1e-14\n"
            "equation.momentum.relativeTolerance 1e-8\nequation.momentum.maxIterations 1000\n"
            "equation.pressureCorrection.kspType cg\nequation.pressureCorrection.pcType hypre\n"
            "equation.pressureCorrection.absoluteTolerance 1e-14\n"
            "equation.pressureCorrection.relativeTolerance 1e-7\n"
            "equation.pressureCorrection.maxIterations 1200\n")
        control = "startTime 0\nendTime 0\ndeltaT 1\n"
        output = "directory results\ntimeName final\nwriteFields U p\n"
    else:
        methods = "time euler\n" + COMMON_METHODS
        solution = (
            "maxIterations 1\nnCorrectors 2\nnonOrthogonalCorrections 0\n"
            "velocityRelaxation 1\ncontinuityTolerance 1e-7\nvelocityTolerance 1e-7\n"
            "momentumTolerance 1e-6\npressureCorrectionTolerance 1e-7\n"
            "equation.momentum.kspType bcgs\nequation.momentum.pcType bjacobi\n"
            "equation.momentum.absoluteTolerance 1e-12\n"
            "equation.momentum.relativeTolerance 1e-8\nequation.momentum.maxIterations 1000\n"
            "equation.pressureCorrection.kspType cg\nequation.pressureCorrection.pcType hypre\n"
            "equation.pressureCorrection.absoluteTolerance 1e-13\n"
            "equation.pressureCorrection.relativeTolerance 1e-8\n"
            "equation.pressureCorrection.maxIterations 1200\n")
        control = "startTime 0\nendTime 0.01\ndeltaT 0.01\n"
        output = "directory results\ntimeName final\nwriteInterval 1\nwriteFields U p\n"
    (numerics / "methods.bs").write_text(methods, encoding="utf-8")
    (numerics / "solution.bs").write_text(solution, encoding="utf-8")
    (case / "control.bs").write_text(control, encoding="utf-8")
    (case / "output.bs").write_text(output, encoding="utf-8")
    (case / "case.bs").write_text(
        f"solver {args.solver}\nmesh mesh/flow.mesh\nfields fields/initial\n"
        "physics physics/flow.bs\nmethods numerics/methods.bs\n"
        "solution numerics/solution.bs\ncontrol control.bs\noutput output.bs\n"
        "ghostLayers 3\n", encoding="utf-8")
    info = json.loads((args.mesh_case.resolve() / "generation.json").read_text())
    info.update(solver=args.solver, case=str(case), simpleOuterCap=(
        args.max_iterations if args.solver == "simple" else None))
    (case / "generation.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
