#!/usr/bin/env python3
"""Create a small, self-contained BabelSim tutorial case without overwriting data."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]


def methods(names: list[str], time: str) -> str:
    return "time " + time + "\n" + "".join(
        f"equation.{name}.{key} {value}\n"
        for name in names for key, value in (
            ("interpolation", "linear"), ("gradient", "greenGauss"),
            ("convection", "upwind"), ("diffusion", "orthogonal")))


def linear(names: list[str]) -> str:
    return "".join(f"equation.{name}.{key} {value}\n"
                   for name in names for key, value in (
                       ("kspType", "bcgs"), ("pcType", "hypre"),
                       ("absoluteTolerance", "1e-14"),
                       ("relativeTolerance", "1e-10"), ("maxIterations", "1000")))


def field(name: str, kind: str, initial: str, boundaries: dict[str, str]) -> str:
    return (f"field {name}\n{{\n    type {kind}\n    location cell\n"
            f"    internal uniform ({initial})\n    boundary\n    {{\n" +
            "".join(f"        {patch} {{ {value} }}\n" for patch, value in boundaries.items()) +
            "    }\n}\n")


def create(template: str, target: Path) -> None:
    # lexists-like check also rejects an existing dangling symlink.
    if target.exists() or target.is_symlink():
        raise ValueError(f"output already exists: {target}; choose a new directory")
    target.mkdir(parents=True)
    try:
        for directory in ("mesh", "fields/initial", "physics", "numerics"):
            (target / directory).mkdir(parents=True, exist_ok=True)
        flow = template in {"simple", "piso"}
        names = ["temperature"] if template in {"heat", "custom-heat"} else (
            ["transport"] if template == "transport" else (
                ["momentum", "pressureCorrection"] if flow else ["temperature", "concentration"]))
        solver = {"custom-heat": "customHeat", "coupled": "coupledScalar"}.get(template, template)
        patches = (["inlet", "outlet", "lower_wall", "upper_wall", "front", "back"] if flow else
                   (["inlet", "outlet", "lower", "upper", "front", "back"] if template == "transport" else
                    ["hot", "cold", "lower", "upper", "front", "back"]))
        kinds = ["inlet", "outlet", "wall", "wall", "symmetry", "symmetry"] if flow else ["generic"] * 6
        command = [sys.executable, str(ROOT / "tools/generate_hexa_mesh.py"),
                   str(target / "mesh/domain.mesh"), "--cells", "32", "8" if flow else "4", "1",
                   "--bounds", "0", "0", "0", "4" if flow else "1", "1" if flow else "0.2", "0.1"]
        for face, patch, kind in zip(("minus_x", "plus_x", "minus_y", "plus_y", "minus_z", "plus_z"), patches, kinds):
            command += ["--patch", face, patch, kind]
        subprocess.run(command, check=True)
        (target / "case.bs").write_text(
            f"solver {solver}\nmesh mesh/domain.mesh\nfields fields/initial\n"
            "physics physics/model.bs\nmethods numerics/methods.bs\n"
            "solution numerics/solution.bs\ncontrol control.bs\noutput output.bs\n")
        (target / "control.bs").write_text("startTime 0\nendTime 0.05\ndeltaT 0.01\n")
        (target / "output.bs").write_text("directory results\ntimeName final\nwriteInterval 1\n")
        (target / "numerics/methods.bs").write_text(methods(names, "steady" if template == "simple" else "euler"))
        algorithm = ""
        if template == "simple":
            algorithm = ("maxIterations 1200\nvelocityRelaxation 0.5\npressureRelaxation 0.3\n"
                         "continuityTolerance 1e-7\nvelocityTolerance 1e-6\n"
                         "momentumTolerance 1e-6\npressureCorrectionTolerance 1e-6\n")
        elif template == "piso":
            algorithm = ("maxIterations 1\nnCorrectors 2\nmaxCorrectors 4\nvelocityRelaxation 1\n"
                         "continuityTolerance 1e-8\n")
        elif template == "coupled":
            algorithm = "couplingIterations 100\ncouplingTolerance 1e-10\n"
        (target / "numerics/solution.bs").write_text(algorithm + linear(names))
        if flow:
            (target / "physics/model.bs").write_text("density 1\ndynamicViscosity 0.01\nturbulenceModel none\n")
            velocity = dict(zip(patches, ("type fixedValue value (1 0 0)", "type zeroGradient",
                                         "type fixedValue value (0 0 0)", "type fixedValue value (0 0 0)",
                                         "type symmetry", "type symmetry")))
            pressure = {patch: "type symmetry" if patch in {"front", "back"} else "type zeroGradient"
                        for patch in patches}
            pressure["outlet"] = "type fixedValue value (0)"
            (target / "fields/initial/U.field").write_text(field("U", "vector", "0 0 0", velocity))
            (target / "fields/initial/p.field").write_text(field("p", "scalar", "0", pressure))
        elif template == "transport":
            (target / "physics/model.bs").write_text("storage 1\ndiffusivity 0.01\nsource 0\n")
            scalar = {patch: "type symmetry" for patch in patches}
            scalar.update(inlet="type fixedValue value (1)", outlet="type zeroGradient")
            velocity = {patch: "type symmetry" for patch in patches}
            velocity.update(inlet="type fixedValue value (1 0 0)", outlet="type zeroGradient")
            (target / "fields/initial/C.field").write_text(field("C", "scalar", "0", scalar))
            (target / "fields/initial/U.field").write_text(field("U", "vector", "1 0 0", velocity))
        else:
            scalar = {patch: "type symmetry" for patch in patches}
            scalar.update(hot="type fixedValue value (1)", cold="type fixedValue value (0)")
            (target / "fields/initial/T.field").write_text(field("T", "scalar", "0", scalar))
            if template == "coupled":
                (target / "physics/model.bs").write_text("diffusivity 0.1\ncoupling 0.2\n")
                scalar.update(hot="type fixedValue value (0)")
                (target / "fields/initial/C.field").write_text(field("C", "scalar", "0", scalar))
                (target / "output.bs").write_text("directory results\ntimeName final\nwriteInterval 1\nwriteFields T C total\n")
            else:
                (target / "physics/model.bs").write_text("density 1\nheatCapacity 1\nconductivity 0.1\nsource 0\n")
    except BaseException:
        # This directory was created by this invocation; no existing data is removed.
        shutil.rmtree(target)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", choices=("heat", "transport", "simple", "piso", "custom-heat", "coupled"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        create(args.template, args.output.absolute())
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"create_case: {error}\n")
    print(f"Created {args.template} case: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
