#!/usr/bin/env python3
"""计算并绘制方腔中心线的 Ghia 误差与逐级网格差异。"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from validate_cavity import GHIA_U, centreline, ghia_errors


def _result(value: str) -> tuple[int, int, Path]:
    try:
        label, path = value.split("=", 1)
        reynolds, cells = (int(item) for item in label.split(":", 1))
    except ValueError as error:
        raise argparse.ArgumentTypeError("结果必须写成 Re:N=结果路径") from error
    if reynolds not in GHIA_U or cells <= 0:
        raise argparse.ArgumentTypeError("Re 没有 Ghia 数据，或 N 不是正整数")
    return reynolds, cells, Path(path)


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(values * values)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", action="append", type=_result, required=True,
                        help="格式 Re:N=U.csv（并行时也可给结果目录）")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cases: dict[int, dict[int, Path]] = defaultdict(dict)
    for reynolds, cells, path in args.result:
        if cells in cases[reynolds]:
            raise ValueError(f"Re={reynolds}, N={cells} 被重复指定")
        cases[reynolds][cells] = path

    figure, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    dense = np.linspace(0.0, 1.0, 1001)
    for reynolds in sorted(cases):
        resolutions = sorted(cases[reynolds])
        if len(resolutions) < 2:
            raise ValueError(f"Re={reynolds} 至少需要两个网格")
        ghia_u: list[float] = []
        ghia_v: list[float] = []
        lines: dict[int, tuple[np.ndarray, ...]] = {}
        for cells in resolutions:
            path = cases[reynolds][cells]
            _, _, u_error, v_error, v_mask = ghia_errors(path, reynolds)
            ghia_u.append(_rms(u_error))
            ghia_v.append(_rms(v_error[v_mask]))
            lines[cells] = centreline(path)
            print(
                f"Re={reynolds} N={cells} "
                f"Ghia_u_L2={ghia_u[-1]:.9g} Ghia_v_L2={ghia_v[-1]:.9g}")

        axes[0, 0].loglog(resolutions, ghia_u, "o-", label=f"Re={reynolds}")
        axes[0, 1].loglog(resolutions, ghia_v, "o-", label=f"Re={reynolds}")
        fine_resolutions: list[int] = []
        delta_u: list[float] = []
        delta_v: list[float] = []
        for coarse, fine in zip(resolutions, resolutions[1:]):
            coarse_y, coarse_u, coarse_x, coarse_v = lines[coarse]
            fine_y, fine_u, fine_x, fine_v = lines[fine]
            delta_u.append(_rms(
                np.interp(dense, coarse_y, coarse_u) -
                np.interp(dense, fine_y, fine_u)))
            delta_v.append(_rms(
                np.interp(dense, coarse_x, coarse_v) -
                np.interp(dense, fine_x, fine_v)))
            fine_resolutions.append(fine)
            print(
                f"Re={reynolds} N={coarse}->{fine} "
                f"grid_u_L2={delta_u[-1]:.9g} grid_v_L2={delta_v[-1]:.9g}")
        axes[1, 0].loglog(fine_resolutions, delta_u, "s--", label=f"Re={reynolds}")
        axes[1, 1].loglog(fine_resolutions, delta_v, "s--", label=f"Re={reynolds}")

    titles = (
        "Ghia RMS: vertical-centreline u",
        "Ghia RMS: horizontal-centreline v",
        "Successive-grid RMS: u",
        "Successive-grid RMS: v",
    )
    for axis, title in zip(axes.flat, titles):
        axis.set(title=title, xlabel="cells per direction, N", ylabel="RMS error")
        axis.grid(True, which="both", color="0.9", linewidth=0.7)
        axis.legend(frameon=False)
    figure.suptitle("Lid-driven cavity grid convergence")
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
