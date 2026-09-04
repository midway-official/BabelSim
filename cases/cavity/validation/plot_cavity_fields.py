#!/usr/bin/env python3
"""绘制多个 Reynolds 数方腔结果的速度流线与涡量。"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from validate_cavity import GHIA_U, _grid


def _result(value: str) -> tuple[int, Path]:
    try:
        reynolds, path = value.split("=", 1)
        result = int(reynolds), Path(path)
    except ValueError as error:
        raise argparse.ArgumentTypeError("结果必须写成 Re=结果路径") from error
    if result[0] not in GHIA_U:
        raise argparse.ArgumentTypeError(f"没有 Re={result[0]} 的基准定义")
    return result


def _uniform_field(
    x: np.ndarray, y: np.ndarray, value: np.ndarray, points: int = 161,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """把可分离的非均匀结构网格双线性重采样到均匀绘图网格。"""
    uniform_x = np.linspace(float(x[0, 0]), float(x[0, -1]), points)
    uniform_y = np.linspace(float(y[0, 0]), float(y[-1, 0]), points)
    along_x = np.array([np.interp(uniform_x, x[row], value[row])
                        for row in range(value.shape[0])])
    uniform = np.array([np.interp(uniform_y, y[:, 0], along_x[:, column])
                        for column in range(points)]).T
    return uniform_x, uniform_y, uniform


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", action="append", type=_result, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    results = dict(args.result)
    figure, axes = plt.subplots(len(results), 2, figsize=(10, 3.6 * len(results)),
                                squeeze=False)
    for row, reynolds in enumerate(sorted(results)):
        x, y, u, v = _grid(results[reynolds])
        uniform_x, uniform_y, uniform_u = _uniform_field(x, y, u)
        _, _, uniform_v = _uniform_field(x, y, v)
        speed = np.hypot(uniform_u, uniform_v)
        vorticity = (
            np.gradient(uniform_v, uniform_x, axis=1) -
            np.gradient(uniform_u, uniform_y, axis=0))

        speed_axis, vorticity_axis = axes[row]
        speed_plot = speed_axis.contourf(
            uniform_x, uniform_y, speed, levels=36, cmap="viridis")
        speed_axis.streamplot(
            uniform_x, uniform_y, uniform_u, uniform_v,
            color="white", density=1.25, linewidth=0.55, arrowsize=0.65)
        figure.colorbar(speed_plot, ax=speed_axis, label=r"$|U|/U_{lid}$")
        speed_axis.set_title(f"Re={reynolds}: speed and streamlines")

        limit = max(float(np.percentile(np.abs(vorticity), 99.0)), 1e-12)
        vorticity_plot = vorticity_axis.contourf(
            uniform_x, uniform_y, np.clip(vorticity, -limit, limit),
            levels=np.linspace(-limit, limit, 41), cmap="coolwarm", extend="both")
        figure.colorbar(vorticity_plot, ax=vorticity_axis,
                        label=r"$(\partial v/\partial x-\partial u/\partial y)L/U_{lid}$")
        vorticity_axis.set_title(f"Re={reynolds}: vorticity (99% clipped)")

        for axis in (speed_axis, vorticity_axis):
            axis.set(xlabel=r"$x/L$", ylabel=r"$y/L$", aspect="equal")
    figure.suptitle("BabelSim lid-driven cavity flow fields")
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
