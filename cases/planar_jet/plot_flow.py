#!/usr/bin/env python3
"""Plot final planar-jet vorticity/velocity and the saved transverse-velocity probe."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE = Path(__file__).resolve().parent
RESULTS = CASE / "results"
NX, NY = 500, 200
X0, X1 = 0.0, 20.0
Y0, Y1 = -5.0, 5.0


def field(time_dir: Path, name: str, width: int) -> np.ndarray:
    rank_dirs = sorted(p for p in time_dir.glob("rank-*") if p.is_dir())
    arrays = []
    for rank_dir in rank_dirs:
        path = rank_dir / f"{name}.csv"
        if path.exists():
            arrays.append(np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2))
    if not arrays:
        raise FileNotFoundError(f"no {name}.csv under {time_dir}")
    data = np.concatenate(arrays, axis=0)
    data = data[np.argsort(data[:, 0].astype(np.int64))]
    if data.shape[0] != NX * NY or data.shape[1] != width:
        raise ValueError(f"unexpected {name} shape {data.shape} in {time_dir}")
    ids = data[:, 0].astype(np.int64)
    if not np.array_equal(ids, np.arange(NX * NY)):
        raise ValueError(f"global cell IDs are incomplete in {time_dir}")
    return data


def main() -> None:
    time_dirs = []
    for path in RESULTS.iterdir():
        if not path.is_dir():
            continue
        try:
            float(path.name)
        except ValueError:
            continue
        time_dirs.append(path)
    time_dirs.sort(key=lambda p: float(p.name))
    if not time_dirs:
        raise FileNotFoundError(f"no saved solver snapshots under {RESULTS}")
    final_dir = time_dirs[-1]
    udata = field(final_dir, "U", 7)
    centres = udata[:, 1:4]
    velocity = udata[:, 4:7].reshape(NY, NX, 3)
    x = centres[:, 0].reshape(NY, NX)
    y = centres[:, 1].reshape(NY, NX)
    ux, uy = velocity[:, :, 0], velocity[:, :, 1]
    dy = (Y1 - Y0) / NY
    dx = (X1 - X0) / NX
    omega = np.gradient(uy, dx, axis=1, edge_order=2) - np.gradient(ux, dy, axis=0, edge_order=2)
    limit = max(float(np.quantile(np.abs(omega), 0.995)), 1.0e-12)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), constrained_layout=True)
    levels = np.linspace(-limit, limit, 61)
    image = axes[0].contourf(x, y, omega, levels=levels, cmap="RdBu_r", extend="both")
    axes[0].streamplot(x[0], y[:, 0], ux, uy, color="black", density=1.1, linewidth=0.45, arrowsize=0.55)
    axes[0].set(
        xlim=(X0, X1), ylim=(Y0, Y1), xlabel="x/D", ylabel="y/D",
        title=rf"$\omega_z=\partial U_y/\partial x-\partial U_x/\partial y$ at t={float(final_dir.name):g}",
    )
    fig.colorbar(image, ax=axes[0], label=r"$\omega_z D/U_j$")

    for target in [1.0, 2.0, 4.0, 8.0, 12.0, 16.0]:
        i = int(np.argmin(np.abs(x[0] - target)))
        axes[1].plot(ux[:, i], y[:, 0], label=f"x/D={x[0, i]:.1f}")
    axes[1].set(
        xlabel="$U_x/U_j$", ylabel="y/D", xlim=(0.0, 1.1), ylim=(Y0, Y1),
        title="Instantaneous velocity profiles",
    )
    axes[1].grid(alpha=0.25)
    axes[1].legend(ncol=2, fontsize=8)
    figure_path = CASE / "validation/planar_jet_final.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_path, dpi=170)
    plt.close(fig)

    probe_target = (5.0, 0.5)
    probe = []
    for time_dir in time_dirs:
        data = field(time_dir, "U", 7)
        xy = data[:, 1:3]
        index = int(np.argmin((xy[:, 0] - probe_target[0]) ** 2 + (xy[:, 1] - probe_target[1]) ** 2))
        probe.append((float(time_dir.name), float(data[index, 5])))
    probe = np.asarray(probe)
    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    ax.plot(probe[:, 0], probe[:, 1], marker="o", markersize=2.5, linewidth=1)
    ax.set(xlabel="t Uj/D", ylabel="$U_y/U_j$", title="Transverse-velocity probe at (5D, 0.5D)")
    ax.grid(alpha=0.25)
    probe_path = CASE / "validation/planar_jet_probe.png"
    fig.savefig(probe_path, dpi=170)
    plt.close(fig)

    diagnostics = {
        "case": "2-D laminar plane jet, no turbulence model, transient PISO/BDF2",
        "final_time_D_over_Uj": float(final_dir.name),
        "snapshot_count": len(time_dirs),
        "first_snapshot_time_D_over_Uj": float(time_dirs[0].name),
        "last_snapshot_time_D_over_Uj": float(final_dir.name),
        "mesh_cells": NX * NY,
        "Ux_range_Uj": [float(np.min(ux)), float(np.max(ux))],
        "Uy_range_Uj": [float(np.min(uy)), float(np.max(uy))],
        "maximum_abs_omega_z_D_over_Uj": float(np.max(np.abs(omega))),
        "omega_plot_symmetric_limit_99_5_percentile": limit,
        "probe": {
            "location_D": list(probe_target),
            "component": "Uy/Uj",
            "samples": int(len(probe)),
            "sample_times_and_values": probe.tolist(),
            "sampling_interval_D_over_Uj": float(np.median(np.diff(probe[:, 0]))) if len(probe) > 1 else None,
        },
        "artifacts": [str(figure_path.relative_to(CASE)), str(probe_path.relative_to(CASE))],
    }
    output = CASE / "validation/diagnostics.json"
    output.write_text(json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
