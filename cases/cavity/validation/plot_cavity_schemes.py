#!/usr/bin/env python3
"""绘制同一 Reynolds 数下多组方腔中心线结果。"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from validate_cavity import GHIA_U, GHIA_V, GHIA_X, GHIA_Y, centreline


def _result(value: str) -> tuple[str, Path]:
    try:
        label, path = value.split("=", 1)
    except ValueError as error:
        raise argparse.ArgumentTypeError("结果必须写成 标签=结果路径") from error
    if not label:
        raise argparse.ArgumentTypeError("格式标签不能为空")
    return label, Path(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--re", type=int, choices=sorted(GHIA_U), required=True)
    parser.add_argument("--result", action="append", type=_result, required=True)
    parser.add_argument("--title", help="覆盖默认图标题")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    figure, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for label, path in args.result:
        y, u, x, v = centreline(path)
        axes[0].plot(u, y, linewidth=1.8, label=label)
        axes[1].plot(x, v, linewidth=1.8, label=label)
    axes[0].scatter(GHIA_U[args.re], GHIA_Y, marker="o", facecolors="none",
                    color="black", label="Ghia et al.", zorder=4)
    axes[1].scatter(GHIA_X, GHIA_V[args.re], marker="o", facecolors="none",
                    color="black", label="Ghia et al.", zorder=4)
    if args.re == 400:
        anomaly = int(np.flatnonzero(np.isclose(GHIA_X, .9063))[0])
        axes[1].scatter([GHIA_X[anomaly]], [GHIA_V[args.re][anomaly]],
                        color="#ff7f0e", marker="x", s=55,
                        label="published anomalous point", zorder=5)
    axes[0].set(xlabel=r"$u/U_{lid}$", ylabel=r"$y/L$",
                title="Vertical centreline")
    axes[1].set(xlabel=r"$x/L$", ylabel=r"$v/U_{lid}$",
                title="Horizontal centreline")
    for axis in axes:
        axis.grid(True, color="0.9", linewidth=0.7)
        axis.legend(frameon=False)
    figure.suptitle(args.title or f"Lid-driven cavity profiles, Re={args.re}")
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
