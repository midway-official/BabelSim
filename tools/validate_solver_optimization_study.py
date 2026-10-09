#!/usr/bin/env python3
"""Validate frozen binary A/B fields, coverage and saved time after measurement.

Settings overrides need their own accuracy budget and are reported as skipped.
This checks numerical equivalence; it does not validate a physical model.
"""

import argparse
import json
import math
from pathlib import Path
import re

from compare_parallel_results import read_result


def saved_state(run, require_end_time=False):
    root = Path(run["case"]) / "results" / f"benchmark-r{run['ranks']}-n{run['repeat']:02d}"
    ranks = sorted(p for p in root.glob("rank-*") if p.is_dir())
    if len(ranks) != run["ranks"]:
        raise ValueError(f"rank coverage differs: {root}")
    expected = run["performance"][0]["partition"]["globalCells"]
    times = []
    for rank in ranks:
        rows = [line.split() for line in (rank / "metadata.bs").read_text().splitlines()]
        count = next(int(row[1]) for row in rows if row[0] == "global_cell_count")
        times.append(next(float(row[1]) for row in rows if row[0] == "time"))
        if count != expected:
            raise ValueError(f"metadata global count differs: {rank}")
    if not all(math.isfinite(t) and abs(t - times[0]) <= 1e-12 for t in times):
        raise ValueError(f"inconsistent saved time: {root}")
    if require_end_time:
        match = re.search(r"(?m)^endTime\s+(\S+)", (Path(run["case"]) / "control.bs").read_text())
        expected_time = float(match[1]) if match else float("nan")
        if not math.isfinite(expected_time) or abs(times[0] - expected_time) > 1e-12 * max(1.0, abs(expected_time)):
            raise ValueError(f"saved PISO time does not reach the input end time: {root}")
    fields = read_result(root)
    for name, field in fields.items():
        if len(field) != expected or min(field) != 0 or max(field) != expected - 1:
            raise ValueError(f"incomplete global ID range: {root}/{name}")
    return root, times[0], fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("--reference", default="baseline")
    parser.add_argument("--atol", type=float, default=1e-10)
    parser.add_argument("--rtol", type=float, default=1e-9)
    args = parser.parse_args()
    if not all(math.isfinite(t) and t >= 0 for t in (args.atol, args.rtol)):
        parser.error("tolerances must be nonnegative and finite")
    runs = [r for r in json.loads((args.study / "records.json").read_text()) if not r["warmup"]]
    required_steps = json.loads((args.study / "study.json").read_text()).get("acceptedSteps", {})
    references = {(r["caseLabel"], r["ranks"]): r for r in reversed(runs)
                  if r["variant"] == args.reference and not r.get("settingsOverride")}
    comparisons, skipped = [], []
    for key, reference_run in references.items():
        matched = [r for r in runs if (r["caseLabel"], r["ranks"]) == key]
        reference_path, reference_time, reference = saved_state(reference_run, key[0] in required_steps)
        for run in matched:
            if run.get("settingsOverride"):
                skipped.append(dict(variant=run["variant"], case=key[0], ranks=key[1],
                                    repeat=run["repeat"], reason="settings accuracy checked separately"))
                continue
            if not run["validWork"] or not run.get("measurementValid", True):
                raise ValueError("unsuccessful or invalid measured work")
            if run["inputSha256"] != reference_run["inputSha256"] or run["workUnits"] != reference_run["workUnits"]:
                raise ValueError("binary A/B inputs or completed work differ")
            candidate_path, candidate_time, candidate = saved_state(run, key[0] in required_steps)
            if abs(candidate_time - reference_time) > 1e-12 or candidate.keys() != reference.keys():
                raise ValueError("saved time or field sets differ")
            fields = {}
            for name, expected in reference.items():
                actual = candidate[name]
                error2 = scale2 = maximum = maximum_ratio = 0.0
                for identifier, values in expected.items():
                    other = actual[identifier]
                    if len(values) != len(other):
                        raise ValueError("field component counts differ")
                    for a, b in zip(values, other):
                        delta = abs(a - b)
                        limit = args.atol + args.rtol * max(abs(a), abs(b))
                        if delta > limit:
                            raise ValueError(f"{candidate_path}/{name}[{identifier}] exceeds the component budget")
                        error2 += delta * delta
                        scale2 += a * a
                        maximum = max(maximum, delta)
                        maximum_ratio = max(maximum_ratio, delta / limit if limit else 0.0)
                fields[name] = dict(cells=len(actual), maxAbsolute=maximum,
                                    maxToleranceRatio=maximum_ratio,
                                    relativeL2=math.sqrt(error2 / scale2) if scale2 else None)
            comparisons.append(dict(variant=run["variant"], case=key[0], ranks=key[1],
                repeat=run["repeat"], time=candidate_time, reference=str(reference_path),
                candidate=str(candidate_path), fields=fields))
            print(f"PASS {run['variant']} {key[0]} ranks={key[1]} repeat={run['repeat']}", flush=True)
    if len(comparisons) + len(skipped) != len(runs) or not comparisons:
        raise ValueError("not all measured runs have a matching reference")
    result = dict(status="passed", atol=args.atol, rtol=args.rtol,
                  comparisons=comparisons, skipped=skipped, note=__doc__)
    (args.study / "field-validation.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"Validated {len(comparisons)} A/B runs; {len(skipped)} settings runs need separate validation")


if __name__ == "__main__":
    main()
