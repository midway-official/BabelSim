#!/usr/bin/env python3
"""Compare owned-cell CSV files from two MPI runs by global cell id."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def load(directory: Path, prefix: str, ranks: int) -> dict[int, np.ndarray]:
    result: dict[int, np.ndarray] = {}
    for rank in range(ranks):
        path = directory / f"{prefix}_{rank}.csv"
        data = np.genfromtxt(path, delimiter=",", names=True)
        if data.dtype.names is None or "global_id" not in data.dtype.names:
            raise ValueError(f"{path} has no global_id column")
        for row in np.atleast_1d(data):
            identifier = int(row["global_id"])
            if identifier in result:
                raise ValueError(f"duplicate global id {identifier}")
            result[identifier] = np.array(
                [row[name] for name in ("x", "y", "z", "u", "v", "w", "p")],
                dtype=float,
            )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--reference-ranks", type=int, required=True)
    parser.add_argument("--candidate-ranks", type=int, required=True)
    parser.add_argument("--prefix", default="poiseuille")
    parser.add_argument("--atol", type=float, default=5e-6)
    parser.add_argument("--rtol", type=float, default=5e-6)
    args = parser.parse_args()
    reference = load(args.reference, args.prefix + str(args.reference_ranks), args.reference_ranks)
    candidate = load(args.candidate, args.prefix + str(args.candidate_ranks), args.candidate_ranks)
    if set(reference) != set(candidate):
        raise ValueError("MPI runs do not cover the same global cell ids")
    errors = np.array(
        [np.abs(candidate[index] - reference[index]) for index in sorted(reference)])
    scales = np.array(
        [np.maximum(np.abs(reference[index]), 1.0) for index in sorted(reference)])
    limit = args.atol + args.rtol * scales
    maximum = float(np.max(errors))
    print(f"cells={errors.shape[0]} max_abs_difference={maximum:.8g}")
    if np.any(errors > limit):
        raise SystemExit("MPI field difference exceeds configured tolerance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
