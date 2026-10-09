#!/usr/bin/env python3
"""Freeze the Re3900 production input for a bounded, reproducible speed study.

This does not edit the production run. An optional saved state is a new initial
condition with Euler startup, not a BDF2 history restart.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--turbulence-correctors", type=int)
    parser.add_argument("--dt", type=float)
    parser.add_argument("--set", action="append", default=[], help="solution key=value")
    parser.add_argument("--snapshot", type=Path, help="saved state to use as a new initial condition")
    args = parser.parse_args()
    if args.steps < 1 or args.output.exists():
        parser.error("positive step count and a new output directory are required")
    args.output.mkdir(parents=True)
    for name in ("case.bs", "control.bs", "output.bs"):
        shutil.copy2(args.source / name, args.output / name)
    for name in ("mesh", "fields", "physics", "numerics"):
        shutil.copytree(args.source / name, args.output / name)
    original_control = (args.source / "control.bs").read_text()
    dt = float(re.search(r"(?m)^deltaT\s+(\S+)", original_control)[1])
    if args.dt is not None:
        if not math.isfinite(args.dt) or args.dt <= 0:
            parser.error("dt must be positive and finite")
        dt = args.dt
    start = float(re.search(r"(?m)^startTime\s+(\S+)", original_control)[1])
    (args.output / "control.bs").write_text(
        f"startTime {start:.17g}\nendTime {start + dt * args.steps:.17g}\ndeltaT {dt:.17g}\n")
    (args.output / "output.bs").write_text(
        "directory results\ntimeName benchmark\n"
        f"writeInterval {args.steps}\nwriteFields U p k omega mut\n")
    if args.turbulence_correctors is not None:
        if args.turbulence_correctors < 1:
            parser.error("turbulence-correctors must be positive")
        solution = args.output / "numerics/solution.bs"
        solution.write_text(solution.read_text() +
                            f"\nturbulenceCorrectors {args.turbulence_correctors}\n")
    solution = args.output / "numerics/solution.bs"
    for item in args.set:
        key, value = item.split("=", 1)
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", key) or not value.strip():
            parser.error(f"invalid solution override: {item}")
        text = solution.read_text()
        pattern = rf"(?m)^{re.escape(key)}\s+[^\n]+$"
        text = re.sub(pattern, lambda _: f"{key} {value}", text) if re.search(pattern, text) else text + f"\n{key} {value}\n"
        solution.write_text(text)
    if args.snapshot is not None:
        ranks = sorted(args.snapshot.glob("rank-*"))
        if not ranks:
            parser.error("snapshot contains no rank directories")
        n = int(re.search(r"(?m)^global_cell_count\s+(\d+)",
                          (ranks[0] / "metadata.bs").read_text())[1])
        for definition in sorted((args.output / "fields/initial").glob("*.field")):
            name = definition.stem
            if not (ranks[0] / f"{name}.csv").exists():
                continue
            data = {}
            for rank in ranks:
                with (rank / f"{name}.csv").open() as stream:
                    reader = csv.DictReader(stream)
                    components = [k for k in reader.fieldnames if k.startswith("value")]
                    for row in reader:
                        gid = int(row["global_id"])
                        value = [float(row[c]) for c in components]
                        if gid in data or not all(map(math.isfinite, value)):
                            raise ValueError("duplicate ID or nonfinite initial state")
                        data[gid] = value
            if set(data) != set(range(n)):
                raise ValueError(f"incomplete snapshot field: {name}")
            filename = f"{name}.benchmark.dat"
            with (definition.parent / filename).open("w") as stream:
                for gid in range(n):
                    stream.write(f"{gid} " + " ".join(f"{v:.17g}" for v in data[gid]) + "\n")
            definition.write_text(re.sub(r"(?m)^\s*internal\s+[^\n]+",
                f"    internal file {filename}", definition.read_text()))
    manifest = {str(p.relative_to(args.output)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(args.output.rglob("*")) if p.is_file()}
    (args.output / "preparation.json").write_text(json.dumps(dict(
        source=str(args.source.resolve()), deltaT=dt, steps=args.steps,
        turbulenceCorrectors=args.turbulence_correctors, settings=args.set,
        snapshot=str(args.snapshot.resolve()) if args.snapshot else None,
        startup="Euler first step, then configured scheme", sha256=manifest), indent=2) + "\n")
    print(json.dumps(dict(case=str(args.output.resolve()), deltaT=dt, steps=args.steps)))


if __name__ == "__main__":
    main()
