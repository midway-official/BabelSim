#!/usr/bin/env python3
"""Compare two BabelSim rank-directory results by global cell id."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_result(directory: Path) -> dict[str, dict[int, tuple[float, ...]]]:
    ranks = sorted(path for path in directory.glob("rank-*") if path.is_dir())
    if not ranks:
        raise ValueError(f"{directory} contains no rank directories")
    fields: dict[str, dict[int, tuple[float, ...]]] = {}
    for rank in ranks:
        for path in sorted(rank.glob("*.csv")):
            values = fields.setdefault(path.stem, {})
            with path.open(newline="") as stream:
                reader = csv.DictReader(stream)
                required = {"global_id", "x", "y", "z"}
                if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                    raise ValueError(f"{path} has an invalid BabelSim field header")
                components = [name for name in reader.fieldnames if name.startswith("value")]
                if not components:
                    raise ValueError(f"{path} has no field values")
                for row in reader:
                    identifier = int(row["global_id"])
                    if identifier in values:
                        raise ValueError(f"{path} duplicates global id {identifier}")
                    value = tuple(float(row[name]) for name in components)
                    if not all(math.isfinite(component) for component in value):
                        raise ValueError(f"{path} contains a non-finite value")
                    values[identifier] = value
    return fields


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--atol", type=float, default=5e-6)
    parser.add_argument("--rtol", type=float, default=5e-6)
    args = parser.parse_args()
    if args.atol < 0 or args.rtol < 0:
        raise ValueError("tolerances must be non-negative")
    reference = read_result(args.reference)
    candidate = read_result(args.candidate)
    if reference.keys() != candidate.keys():
        raise ValueError("field sets differ between reference and candidate")
    for name in sorted(reference):
        expected = reference[name]
        actual = candidate[name]
        if expected.keys() != actual.keys():
            raise ValueError(f"global-id coverage differs for field {name}")
        maximum = 0.0
        for identifier, reference_value in expected.items():
            candidate_value = actual[identifier]
            if len(reference_value) != len(candidate_value):
                raise ValueError(f"component count differs for field {name}")
            for lhs, rhs in zip(reference_value, candidate_value):
                maximum = max(maximum, abs(lhs - rhs))
                if abs(lhs - rhs) > args.atol + args.rtol * max(abs(lhs), abs(rhs)):
                    raise SystemExit(
                        f"{name}[{identifier}] differs: reference={lhs:.17g} candidate={rhs:.17g}"
                    )
        print(f"{name}: cells={len(expected)} max_abs_difference={maximum:.8g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
