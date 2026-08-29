#!/usr/bin/env python3
"""Validate a BabelSim Re=100 cavity CSV against Ghia centreline data."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


GHIA_Y = np.array(
    [1.0, .9766, .9688, .9609, .9531, .8516, .7344, .6172,
     .5, .4531, .2813, .1719, .1016, .0703, .0625, .0547, 0.0]
)
GHIA_U = np.array(
    [1.0, .84123, .78871, .73722, .68717, .23151, .00332, -.13641,
     -.20581, -.21090, -.15662, -.10150, -.06434, -.04775, -.04192,
     -.03717, 0.0]
)
GHIA_X = np.array(
    [1.0, .9688, .9609, .9531, .9453, .9063, .8594, .8047,
     .5, .2344, .2266, .1563, .0938, .0781, .0703, .0625, 0.0]
)
GHIA_V = np.array(
    [0.0, -.05906, -.07391, -.08864, -.10313, -.16914, -.22445,
     -.24533, .05454, .17527, .17507, .16077, .12317, .10890, .10091,
     .09233, 0.0]
)


def _grid(csv: Path) -> tuple[np.ndarray, ...]:
    data = np.genfromtxt(csv, delimiter=",", names=True)
    required = {"x", "y", "u", "v"}
    if data.dtype.names is None or not required.issubset(data.dtype.names):
        raise ValueError(f"{csv} must contain columns {sorted(required)}")
    if not all(np.isfinite(data[name]).all() for name in required):
        raise ValueError("cavity CSV contains non-finite values")
    # Polyhedral centroids can differ by a few ulps between otherwise aligned
    # rows, so topology inference must not rely on bitwise-equal coordinates.
    rounded_x = np.round(data["x"], 12)
    rounded_y = np.round(data["y"], 12)
    xs = np.unique(rounded_x)
    ys = np.unique(rounded_y)
    if xs.size * ys.size != data.size:
        raise ValueError("cavity CSV is not one structured z layer")
    order = np.lexsort((rounded_x, rounded_y))
    shape = (ys.size, xs.size)
    return (
        data["x"][order].reshape(shape),
        data["y"][order].reshape(shape),
        data["u"][order].reshape(shape),
        data["v"][order].reshape(shape),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--u-l2-limit", type=float, default=0.01)
    parser.add_argument("--v-l2-limit", type=float, default=0.01)
    args = parser.parse_args()
    if GHIA_Y.shape != GHIA_U.shape or GHIA_X.shape != GHIA_V.shape:
        raise RuntimeError("Ghia coordinate and reference arrays have different lengths")
    x, y, u, v = _grid(args.csv)

    u_line = np.array([np.interp(0.5, x[row], u[row]) for row in range(y.shape[0])])
    v_line = np.array([np.interp(0.5, y[:, column], v[:, column])
                       for column in range(x.shape[1])])
    calculated_u = np.interp(GHIA_Y, np.r_[0.0, y[:, 0], 1.0],
                             np.r_[0.0, u_line, 1.0])
    calculated_v = np.interp(GHIA_X, np.r_[0.0, x[0], 1.0],
                             np.r_[0.0, v_line, 0.0])
    u_error = calculated_u - GHIA_U
    v_error = calculated_v - GHIA_V
    u_l2 = float(np.sqrt(np.mean(u_error * u_error)))
    v_l2 = float(np.sqrt(np.mean(v_error * v_error)))
    print(
        f"grid={x.shape[1]}x{y.shape[0]} "
        f"u_center={calculated_u[8]:.8f} v_center={calculated_v[8]:.8f} "
        f"u_Linf={np.max(np.abs(u_error)):.8g} u_L2={u_l2:.8g} "
        f"v_Linf={np.max(np.abs(v_error)):.8g} v_L2={v_l2:.8g}"
    )
    if u_l2 > args.u_l2_limit or v_l2 > args.v_l2_limit:
        raise SystemExit("Ghia centreline error exceeds the configured limit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
