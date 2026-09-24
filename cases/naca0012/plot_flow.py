#!/usr/bin/env python3
"""Plot cell-resolved speed and streamlines from a saved NACA0012 result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.tri import LinearTriInterpolator

ROOT = Path(__file__).resolve().parent


def read_field(snapshot: Path, name: str) -> np.ndarray:
    files = sorted(snapshot.glob(f"rank-*/{name}.csv"))
    if not files:
        raise FileNotFoundError(f"no {name}.csv files under {snapshot}")
    data = np.vstack([np.loadtxt(path, delimiter=",", skiprows=1) for path in files])
    data = data[np.argsort(data[:, 0])]
    ids = data[:, 0].astype(np.int64)
    if not np.array_equal(ids, np.arange(len(ids))):
        raise ValueError(f"{name}: expected one field value per global cell")
    return data


def polygon_centres(points: list[tuple[float, float]], cells: list[tuple[int, ...]]) -> np.ndarray:
    centres = np.empty((len(cells), 2), dtype=np.float64)
    for i, cell in enumerate(cells):
        polygon = np.asarray([points[index] for index in cell])
        following = np.roll(polygon, -1, axis=0)
        cross = polygon[:, 0] * following[:, 1] - following[:, 0] * polygon[:, 1]
        twice_area = cross.sum()
        centres[i, 0] = ((polygon[:, 0] + following[:, 0]) * cross).sum() / (3 * twice_area)
        centres[i, 1] = ((polygon[:, 1] + following[:, 1]) * cross).sum() / (3 * twice_area)
    return centres


def read_mesh(path: Path) -> tuple[list[tuple[float, float]], list[tuple[int, ...]], list[tuple[int, int]]]:
    """Read the exact 2-D cell polygons and airfoil edges from the solver mesh."""
    lines = iter(path.read_text(encoding="utf-8").splitlines())
    if next(lines).strip() != "BABELSIM_MESH 3":
        raise ValueError(f"unsupported mesh format in {path}")
    if next(lines).strip() != "vertices":
        raise ValueError("mesh is missing its vertex section")
    vertex_count = int(next(lines))
    vertices = np.asarray([[float(x) for x in next(lines).split()] for _ in range(vertex_count)])
    if next(lines).strip() != "faces":
        raise ValueError("mesh is missing its face section")
    face_count = int(next(lines))
    cell_vertices: list[set[int]] = [set() for _ in range(int(json.loads((path.parent / "mesh_quality.json").read_text())["cell_count"]))]
    wall_faces: list[tuple[int, ...]] = []
    for _ in range(face_count):
        fields = next(lines).split()
        if fields[0] != "face":
            raise ValueError("invalid mesh face record")
        count, owner, neighbour, patch = map(int, fields[1:5])
        ring = tuple(map(int, fields[5:5 + count]))
        bottom = tuple(v for v in ring if abs(vertices[v, 2]) < 1e-15)
        if owner >= 0:
            cell_vertices[owner].update(bottom)
        if neighbour >= 0:
            cell_vertices[neighbour].update(bottom)
        if neighbour < 0:
            wall_faces.append((patch, *bottom))
    if next(lines).strip() != "patches":
        raise ValueError("mesh is missing its patch section")
    patch_count = int(next(lines))
    patches = {}
    for _ in range(patch_count):
        fields = next(lines).split()
        if len(fields) != 3 or fields[0] != "patch":
            raise ValueError("invalid mesh patch record")
        patches[int(len(patches))] = fields[1]
    if next(lines).strip() != "end":
        raise ValueError("mesh is missing its end marker")

    points = [tuple(row[:2]) for row in vertices if abs(row[2]) < 1e-15]
    cells = []
    for cell_id, cell in enumerate(cell_vertices):
        if len(cell) < 3:
            raise ValueError(f"cell {cell_id} has fewer than three planar vertices")
        centre = vertices[list(cell), :2].mean(axis=0)
        ordered = sorted(cell, key=lambda v: np.arctan2(vertices[v, 1] - centre[1], vertices[v, 0] - centre[0]))
        cells.append(tuple(ordered))
    airfoil_id = next((patch_id for patch_id, name in patches.items() if name == "airfoil"), None)
    if airfoil_id is None:
        raise ValueError("mesh has no airfoil patch")
    wall_edges = [tuple(face[1:]) for face in wall_faces if face[0] == airfoil_id]
    if not wall_edges:
        raise ValueError("mesh has no airfoil wall faces")
    return points, cells, wall_edges


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time", default="50", help="saved time directory (default: 50)")
    parser.add_argument("--mesh", type=Path, default=ROOT / "mesh" / "naca0012.mesh")
    parser.add_argument("--results-root", type=Path, default=ROOT / "results" / "dt0.5_100steps")
    parser.add_argument("--output", type=Path, default=ROOT / "validation" / "flow_t50_dt0.5.png")
    args = parser.parse_args()

    snapshot = args.results_root / args.time
    velocity = read_field(snapshot, "U")
    points, cells, wall_edges = read_mesh(args.mesh)
    if len(cells) != len(velocity):
        raise ValueError(f"mesh has {len(cells)} cells but result has {len(velocity)}")
    centres = polygon_centres(points, cells)
    mismatch = np.max(np.abs(centres - velocity[:, 1:3]))
    if mismatch > 2e-8:
        raise ValueError(f"result and plotting mesh do not match (centre error {mismatch:g})")
    time_label = args.time
    if args.time == "final":
        metadata = snapshot / "rank-0000" / "metadata.bs"
        for line in metadata.read_text(encoding="utf-8").splitlines():
            if line.startswith("time "):
                time_label = line.split(maxsplit=1)[1]
                break

    components = velocity[:, 4:6]
    speed = np.linalg.norm(components, axis=1)
    polygons = [np.asarray([points[index] for index in cell]) for cell in cells]
    wall_segments = np.asarray([[points[a], points[b]] for a, b in wall_edges])

    # Build a valid interpolation triangulation from the actual convex cell
    # polygons. Each cell centre is connected to its edges, so no Delaunay
    # triangle can bridge the airfoil cutout or cross a mesh boundary.
    vertex_points = np.asarray(points)
    vertex_values = np.zeros((len(points), 2), dtype=np.float64)
    vertex_counts = np.zeros(len(points), dtype=np.float64)
    triangle_list: list[tuple[int, int, int]] = []
    cell_offset = len(points)
    for cell_id, cell in enumerate(cells):
        cell_vertices = np.asarray(cell, dtype=np.int64)
        np.add.at(vertex_values, cell_vertices, components[cell_id])
        np.add.at(vertex_counts, cell_vertices, 1.0)
        center_id = cell_offset + cell_id
        triangle_list.extend(
            (center_id, cell[i], cell[(i + 1) % len(cell)])
            for i in range(len(cell))
        )
    vertex_values /= vertex_counts[:, None]
    wall_vertices = {vertex for edge in wall_edges for vertex in edge}
    if wall_vertices:
        vertex_values[list(wall_vertices)] = 0.0
    triangulation_points = np.vstack((vertex_points, centres))
    triangulation_values = np.vstack((vertex_values, components))
    triangulation = mtri.Triangulation(
        triangulation_points[:, 0], triangulation_points[:, 1],
        triangles=np.asarray(triangle_list, dtype=np.int32),
    )
    u_interp = LinearTriInterpolator(triangulation, triangulation_values[:, 0])
    v_interp = LinearTriInterpolator(triangulation, triangulation_values[:, 1])

    fig = plt.figure(figsize=(16, 7))
    grid = fig.add_gridspec(
        1, 3, width_ratios=(1.05, 1.05, 0.035),
        left=0.065, right=0.96, bottom=0.10, top=0.88, wspace=0.16,
    )
    domain_ax = fig.add_subplot(grid[0, 0])
    zoom_ax = fig.add_subplot(grid[0, 1])
    colorbar_ax = fig.add_subplot(grid[0, 2])
    vmax = max(1.8, float(np.ceil(speed.max() * 10) / 10))
    artist = None
    for ax in (domain_ax, zoom_ax):
        artist = PolyCollection(
            polygons,
            array=speed,
            cmap="turbo",
            edgecolors="none",
            linewidths=0,
            antialiaseds=False,
            rasterized=True,
            clim=(0, vmax),
        )
        ax.add_collection(artist)
        ax.add_collection(LineCollection(wall_segments, colors="#151c24", linewidths=1.0, zorder=5))
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x / c")
        ax.set_ylabel("y / c")

    domain_ax.set_xlim(min(p[0] for p in points), max(p[0] for p in points))
    domain_ax.set_ylim(min(p[1] for p in points), max(p[1] for p in points))
    domain_ax.set_title("Full rectangular domain")
    domain_ax.annotate(
        "U∞ = (−1, 0)", xy=(-3.2, 4.15), xytext=(-1.1, 4.15),
        arrowprops={"arrowstyle": "->", "color": "#162c40", "lw": 1.6},
        ha="center", va="bottom", color="#162c40",
    )

    zoom_ax.set_xlim(-1.5, 1.0)
    zoom_ax.set_ylim(-1.0, 1.0)
    zoom_ax.set_title("Airfoil neighbourhood with streamlines")
    gx = np.linspace(-1.5, 1.0, 700)
    gy = np.linspace(-1.0, 1.0, 560)
    grid_x, grid_y = np.meshgrid(gx, gy)
    u_grid = np.ma.masked_invalid(np.asarray(u_interp(grid_x, grid_y).filled(np.nan)))
    v_grid = np.ma.masked_invalid(np.asarray(v_interp(grid_x, grid_y).filled(np.nan)))
    zoom_ax.streamplot(
        gx, gy, u_grid, v_grid, density=1.55, color="#182a3b",
        linewidth=0.55, arrowsize=0.7, minlength=0.08, zorder=4,
    )
    fig.colorbar(artist, cax=colorbar_ax, label="Speed |U| / U∞")
    fig.suptitle(
        f"NACA0012 transient flow — t = {time_label} s, "
        f"{len(cells):,} cells",
        fontsize=15,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=190, bbox_inches="tight")
    print(f"saved {args.output} (max cell-centre mismatch {mismatch:.3g})")


if __name__ == "__main__":
    main()
