#!/usr/bin/env python3
"""Compare the outlet of a BabelSim cell CSV with a Poiseuille parabola."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--y-min", type=float, required=True)
    parser.add_argument("--y-max", type=float, required=True)
    parser.add_argument("--mean-velocity", type=float, default=1.0)
    parser.add_argument("--l2-limit", type=float, default=5e-3)
    args = parser.parse_args()
    if not args.y_max > args.y_min:
        raise ValueError("y-max must exceed y-min")

    data = np.genfromtxt(args.csv, delimiter=",", names=True)
    required = {"x", "y"}
    if data.dtype.names is None or not required.issubset(data.dtype.names):
        raise ValueError(f"{args.csv} must contain columns {sorted(required)}")
    if "u" in data.dtype.names:
        velocity = data["u"]
    elif "value0" in data.dtype.names:
        velocity = data["value0"]
    else:
        raise ValueError(f"{args.csv} must contain u or value0")
    if not all(np.isfinite(data[name]).all() for name in required) or not np.isfinite(velocity).all():
        raise ValueError("Poiseuille CSV contains non-finite values")
    outlet = np.isclose(data["x"], np.max(data["x"]), rtol=0.0, atol=1e-13)
    if not np.any(outlet):
        raise ValueError("Poiseuille CSV contains no outlet cells")
    order = np.argsort(data["y"][outlet])
    y = data["y"][outlet][order]
    u = velocity[outlet][order]
    eta = (y - args.y_min) / (args.y_max - args.y_min)
    exact = 6.0 * args.mean_velocity * eta * (1.0 - eta)
    error = u - exact
    l2 = float(np.sqrt(np.mean(error * error)))
    linf = float(np.max(np.abs(error)))
    print(
        f"outlet_cells={u.size} u_max={np.max(u):.8f} "
        f"Linf={linf:.8g} L2={l2:.8g}"
    )
    if l2 > args.l2_limit:
        raise SystemExit("Poiseuille outlet error exceeds the configured limit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
