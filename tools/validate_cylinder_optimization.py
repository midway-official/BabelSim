#!/usr/bin/env python3
"""Independent owned-field, wake-vorticity and wall-force comparison.

These gates establish numerical equivalence within a stated budget, not
experimental validation of 2D URANS. Wall forces use the existing first-cell
molecular shear estimate consistently on both sides.
"""

import argparse
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
from compare_parallel_results import read_result


def saved_time(root, cells):
    times = []
    for rank in sorted(root.glob("rank-*")):
        rows = [line.split() for line in (rank / "metadata.bs").read_text().splitlines()]
        count = next(int(row[1]) for row in rows if row[0] == "global_cell_count")
        times.append(next(float(row[1]) for row in rows if row[0] == "time"))
        if count != cells:
            raise ValueError("saved cell count differs from geometry")
    if not times or not all(math.isfinite(t) and abs(t - times[0]) <= 1e-12 for t in times):
        raise ValueError("missing or inconsistent saved time")
    return times[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--geometry", type=Path, default=Path(__file__).resolve().parents[1] / "cases/cylinder_re3900_urans")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-relative", type=float, default=1e-5)
    parser.add_argument("--max-force", type=float, default=1e-4)
    args = parser.parse_args()
    if not all(math.isfinite(t) and t >= 0 for t in (args.max_relative, args.max_force)):
        parser.error("error budgets must be nonnegative and finite")
    source = args.geometry / "postprocess.py"
    spec = importlib.util.spec_from_file_location("cylinder_post", source)
    post = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(post)
    geo = post.Geometry(args.geometry)
    affine_error = geo.check_affine()
    times = [saved_time(root, geo.n) for root in (args.reference, args.candidate)]
    if abs(times[0] - times[1]) > 1e-12:
        raise ValueError("reference and candidate have different physical times")
    parameters = json.loads((args.geometry / "case_parameters.json").read_text())
    area = np.load(args.geometry / "mesh/layout.npz")["areas"]
    reference, candidate = read_result(args.reference), read_result(args.candidate)
    if reference.keys() != candidate.keys():
        raise ValueError("field sets differ")
    arrays = []
    for fields in (reference, candidate):
        state = {}
        for name, rows in fields.items():
            if set(rows) != set(range(geo.n)):
                raise ValueError(f"incomplete global ID coverage: {name}")
            state[name] = np.array([rows[i] for i in range(geo.n)])
        for name in ("k", "omega", "mut"):
            if np.any(state[name] < 0) or (name == "omega" and np.any(state[name] == 0)):
                raise ValueError(f"invalid transported/closure field: {name}")
        gradient = geo.gradient(state["U"][:, :2])
        state["vorticity"] = (gradient[:, 0, 1] - gradient[:, 1, 0])[:, None]
        arrays.append(state)
    masks = {"global": np.ones(geo.n, dtype=bool),
             "wake": (geo.c[:, 0] > -.6) & (geo.c[:, 0] < 5) & (abs(geo.c[:, 1]) < 2)}
    metrics = {}
    valid = True
    for name in arrays[0]:
        a, b = arrays[0][name], arrays[1][name]
        row = {"maxAbsolute": float(np.max(abs(a - b)))}
        for region, mask in masks.items():
            for weighted in (False, True):
                weight = area[mask, None] if weighted else 1.0
                error = float(np.sqrt(np.sum(weight * (a[mask] - b[mask])**2)))
                scale = float(np.sqrt(np.sum(weight * a[mask]**2)))
                relative = error / scale if scale else (0.0 if error == 0 else float("inf"))
                row[region + ("AreaL2" if weighted else "L2")] = relative
                valid = valid and relative <= args.max_relative
        metrics[name] = row
    forces = [geo.wall_metrics(s["U"], s["p"], parameters)[0] for s in arrays]
    force_delta = abs(forces[0] - forces[1])
    valid = valid and bool(np.all(force_delta <= args.max_force))
    probes = []
    for position in ((1., .5), (2., .5), (4., 0.)):
        index = int(np.argmin(np.sum((geo.c - position)**2, axis=1)))
        probes.append(dict(globalId=index, position=geo.c[index].tolist(),
            reference={k: arrays[0][k][index].tolist() for k in arrays[0]},
            candidate={k: arrays[1][k][index].tolist() for k in arrays[1]}))
    result = dict(reference=str(args.reference.resolve()), candidate=str(args.candidate.resolve()),
        cells=geo.n, time=times[0], gradientAffineError=affine_error,
        maxRelative=args.max_relative, maxForce=args.max_force,
        fields=metrics, wallForceEstimate=dict(reference=forces[0].tolist(),
            candidate=forces[1].tolist(), absoluteDifference=force_delta.tolist()),
        probes=probes, valid=valid, note=__doc__)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(dict(valid=valid, output=str(args.output.resolve()),
        maxFieldRelative=max(r[k] for r in metrics.values() for k in r if k != "maxAbsolute"),
        maxForceDifference=float(np.max(force_delta)))))
    return 0 if valid else 1


if __name__ == "__main__":
    sys.exit(main())
