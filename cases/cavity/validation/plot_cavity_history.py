#!/usr/bin/env python3
"""从 SIMPLE 日志绘制外迭代收敛历史。"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt


LINE = re.compile(
    r"SIMPLE\s+(\d+)\s+mass=([0-9.eE+-]+)\s+dU=([0-9.eE+-]+)\s+"
    r"linP=([0-9.eE+-]+)")


def _log(value: str) -> tuple[int, Path]:
    try:
        reynolds, path = value.split("=", 1)
        return int(reynolds), Path(path)
    except ValueError as error:
        raise argparse.ArgumentTypeError("日志必须写成 Re=文件") from error


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", action="append", type=_log, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    figure, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    for reynolds, path in sorted(args.log):
        values = [LINE.search(line) for line in path.read_text().splitlines()]
        rows = [match.groups() for match in values if match]
        if not rows:
            raise ValueError(f"{path} 中没有 SIMPLE 收敛记录")
        iteration, mass, velocity, pressure = zip(*(
            (int(item[0]), float(item[1]), float(item[2]), float(item[3]))
            for item in rows))
        for axis, data in zip(axes, (velocity, mass, pressure)):
            axis.semilogy(iteration, data, label=f"Re={reynolds}")

    titles = ("Outer-iteration velocity change", "Relative mass imbalance",
              "Pressure linear relative residual")
    labels = (r"$\Delta U_{rel}$", "mass imbalance", "linear residual")
    for axis, title, label in zip(axes, titles, labels):
        axis.set(title=title, xlabel="SIMPLE iteration", ylabel=label)
        axis.grid(True, which="both", color="0.9", linewidth=0.7)
        axis.legend(frameon=False)
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
