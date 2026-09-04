#!/usr/bin/env python3
"""绘制 BabelSim 方腔中心线速度与 Ghia 基准数据的对照图。"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from validate_cavity import GHIA_U, GHIA_V, GHIA_X, GHIA_Y, centreline


def _result(value: str) -> tuple[int, Path]:
    try:
        reynolds, path = value.split("=", 1)
        result = int(reynolds), Path(path)
    except ValueError as error:
        raise argparse.ArgumentTypeError("结果必须写成 Re=U.csv") from error
    if result[0] not in GHIA_U:
        raise argparse.ArgumentTypeError(f"没有 Re={result[0]} 的 Ghia 数据")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", action="append", type=_result, required=True,
                        help="计算结果，格式为 Re=U.csv；可重复指定")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    results = dict(args.result)
    reynolds_numbers = sorted(results)
    figure, axes = plt.subplots(
        len(reynolds_numbers), 2,
        figsize=(10, 3.35 * len(reynolds_numbers)), squeeze=False)

    for row, reynolds in enumerate(reynolds_numbers):
        y, u_line, x, v_line = centreline(results[reynolds])

        u_axis, v_axis = axes[row]
        u_axis.plot(u_line, y,
                    color="#1764ab", linewidth=2, label="BabelSim")
        u_axis.scatter(GHIA_U[reynolds], GHIA_Y, color="#d62728", marker="o",
                       facecolors="none", linewidths=1.4, label="Ghia et al.", zorder=3)
        u_axis.axvline(0.0, color="0.75", linewidth=0.8)
        u_axis.set(xlabel=r"$u(x=0.5,y)/U_{lid}$", ylabel=r"$y/L$",
                   title=fr"$Re={reynolds}$: vertical centreline")

        v_axis.plot(x, v_line,
                    color="#1764ab", linewidth=2, label="BabelSim")
        v_axis.scatter(GHIA_X, GHIA_V[reynolds], color="#d62728", marker="o",
                       facecolors="none", linewidths=1.4, label="Ghia et al.", zorder=3)
        if reynolds == 400:
            anomaly = int(np.flatnonzero(np.isclose(GHIA_X, .9063))[0])
            v_axis.scatter([GHIA_X[anomaly]], [GHIA_V[reynolds][anomaly]],
                           color="#ff7f0e", marker="x", s=55,
                           label="published anomalous point", zorder=4)
        v_axis.axhline(0.0, color="0.75", linewidth=0.8)
        v_axis.set(xlabel=r"$x/L$", ylabel=r"$v(x,y=0.5)/U_{lid}$",
                   title=fr"$Re={reynolds}$: horizontal centreline")

        for axis in (u_axis, v_axis):
            axis.grid(True, color="0.9", linewidth=0.7)
            axis.legend(frameon=False, fontsize=8)

    figure.suptitle("Lid-driven cavity: BabelSim and Ghia et al. (1982)", fontsize=14)
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
