#!/usr/bin/env python3
"""Create a fixed-color-scale GIF from saved planar-jet vorticity snapshots."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

CASE = Path(__file__).resolve().parent
RESULTS = CASE / "results"
OUTPUT = CASE / "validation/planar_jet_vorticity.gif"
NX, NY = 500, 200
X0, X1 = 0.0, 20.0
Y0, Y1 = -5.0, 5.0
DX = (X1 - X0) / NX
DY = (Y1 - Y0) / NY
FPS = 6


def numeric_time_directories() -> list[Path]:
    directories = []
    for path in RESULTS.iterdir():
        if not path.is_dir():
            continue
        try:
            float(path.name)
        except ValueError:
            continue
        directories.append(path)
    directories.sort(key=lambda path: float(path.name))
    if len(directories) != 60:
        raise RuntimeError(f"expected 60 saved time directories, found {len(directories)}")
    return directories


def read_velocity(time_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    pieces = []
    for rank_dir in sorted(time_dir.glob("rank-*")):
        csv_path = rank_dir / "U.csv"
        if csv_path.is_file():
            pieces.append(
                np.loadtxt(
                    csv_path,
                    delimiter=",",
                    skiprows=1,
                    usecols=(0, 4, 5),
                    ndmin=2,
                )
            )
    if not pieces:
        raise FileNotFoundError(f"no rank-sharded U.csv files under {time_dir}")
    data = np.concatenate(pieces, axis=0)
    data = data[np.argsort(data[:, 0].astype(np.int64))]
    expected = NX * NY
    if data.shape != (expected, 3):
        raise ValueError(f"expected {expected} global cells in {time_dir}, found {data.shape[0]}")
    ids = data[:, 0].astype(np.int64)
    if not np.array_equal(ids, np.arange(expected)):
        raise ValueError(f"global cell IDs are incomplete in {time_dir}")
    ux = data[:, 1].reshape(NY, NX)
    uy = data[:, 2].reshape(NY, NX)
    return ux, uy


def main() -> None:
    time_dirs = numeric_time_directories()
    times = np.asarray([float(path.name) for path in time_dirs])
    fields = []
    for path in time_dirs:
        ux, uy = read_velocity(path)
        omega = np.gradient(uy, DX, axis=1, edge_order=2)
        omega -= np.gradient(ux, DY, axis=0, edge_order=2)
        fields.append(omega.astype(np.float32))
    all_vorticity = np.concatenate([frame.ravel() for frame in fields])
    limit = max(float(np.quantile(np.abs(all_vorticity), 0.995)), 1.0e-12)

    fig, ax = plt.subplots(figsize=(12, 5.8), dpi=100, constrained_layout=True)
    image = ax.imshow(
        fields[0],
        extent=(X0, X1, Y0, Y1),
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
    )
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label(r"$\omega_z D/U_j$")
    ax.set(
        xlim=(X0, X1),
        ylim=(Y0, Y1),
        xlabel="x/D",
        ylabel="y/D",
        title=r"Planar-jet vorticity, $\omega_z=\partial U_y/\partial x-\partial U_x/\partial y$",
    )
    stamp = ax.text(
        0.98,
        0.96,
        "",
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )

    def update(index: int):
        image.set_data(fields[index])
        stamp.set_text(f"t = {times[index]:.1f} D/Uj")
        return image, stamp

    animation = FuncAnimation(
        fig,
        update,
        frames=len(fields),
        interval=1000.0 / FPS,
        blit=True,
        repeat=True,
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    animation.save(OUTPUT, writer=PillowWriter(fps=FPS), dpi=100)
    plt.close(fig)
    report = {
        "animation": str(OUTPUT.relative_to(CASE)),
        "frames": len(fields),
        "fps": FPS,
        "duration_seconds": len(fields) / FPS,
        "time_range_D_over_Uj": [float(times[0]), float(times[-1])],
        "fixed_symmetric_vorticity_range_D_over_Uj": [-limit, limit],
        "dimensions_pixels": [1200, 580],
        "file_size_bytes": OUTPUT.stat().st_size,
    }
    (CASE / "validation/animation_summary.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
