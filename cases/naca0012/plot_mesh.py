#!/usr/bin/env python3
"""Render a full-domain and airfoil-neighbourhood view of the generated mesh."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


ROOT = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "mesh" / "mesh_overview.png")
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location("naca_mesh", ROOT / "generate_mesh.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load generate_mesh.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    points, cells, boundary, _metadata = module.build_2d_mesh()

    edges: set[tuple[int, int]] = set()
    for cell in cells:
        for a, b in zip(cell, cell[1:] + cell[:1]):
            edges.add(tuple(sorted((a, b))))
    segments = [(points[a], points[b]) for a, b in edges]
    wall = [(points[a], points[b]) for (a, b), name in boundary.items() if name == "airfoil"]

    fig, (overview, detail) = plt.subplots(
        1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [1.25, 1]})
    for ax, xlim, ylim in (
        (overview, (module.X_MIN, module.X_MAX), (module.Y_MIN, module.Y_MAX)),
        (detail, (-1.5, 1.5), (-1.5, 1.5)),
    ):
        ax.add_collection(LineCollection(
            segments, colors="#a9b4bf", linewidths=0.22, alpha=0.8, rasterized=True))
        ax.add_collection(LineCollection(wall, colors="#c4473a", linewidths=1.15, zorder=4))
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x / c")
        ax.set_ylabel("y / c")

    overview.set_title("Rectangular domain and Cartesian background")
    overview.annotate(
        "U∞ = (-1, 0)  ←", xy=(1.5, 4.2), xytext=(3.0, 4.2),
        arrowprops={"arrowstyle": "->", "lw": 1.8, "color": "#225ea8"},
        color="#225ea8", ha="center", va="bottom")
    detail.set_title("Airfoil, boundary layers and transition")
    detail.plot([], [], color="#c4473a", label="NACA0012 wall")
    detail.legend(loc="upper right", frameon=True)
    fig.suptitle(f"NACA0012 at 15° — {len(cells):,} cells (hex + prism)", fontsize=15)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=170, bbox_inches="tight")
    fig2, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, limits, title in zip(axes, [(.19,.26,-.10,.015),(-.25,-.05,.06,.16),(-.80,-.67,.15,.23)],
                                  ["Leading edge", "Hexahedral wall-normal layers", "Rounded trailing cap and transition"]):
        ax.add_collection(LineCollection(segments, colors="#506a80", linewidths=.35))
        ax.add_collection(LineCollection(wall, colors="#c4473a", linewidths=1.3))
        ax.set_xlim(*limits[:2]); ax.set_ylim(*limits[2:]); ax.set_aspect("equal")
        ax.set_title(title); ax.set_xlabel("x / c"); ax.set_ylabel("y / c")
    fig2.tight_layout()
    fig2.savefig(args.output.with_name("mesh_wall_details.png"), dpi=180, bbox_inches="tight")
    print(args.output)


if __name__ == "__main__":
    main()
