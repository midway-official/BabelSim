#!/usr/bin/env python3
"""Generate a one-step heat case on a mildly warped triangular-prism mesh.

The face-based v3 mesh is deliberately written in a streaming fashion so the
generator's memory use stays small for multi-million-cell cases.  Prisms are
formed by extruding a checkerboard-triangulated square grid; a smooth interior
warp makes the nonzero-warp cases geometrically unstructured without poor
elements.  Meshes and case inputs are written outside the repository by default.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil


def grid_triangles(n: int):
    """Yield counter-clockwise triangle vertex triples on an n by n grid."""
    width = n + 1
    for j in range(n):
        for i in range(n):
            a = i + width * j
            b = a + 1
            d = a + width
            c = d + 1
            if (i + j) % 2 == 0:
                yield (a, b, c)
                yield (a, c, d)
            else:
                yield (a, b, d)
                yield (b, c, d)


def warped_point(i: int, j: int, n: int, warp: float) -> tuple[float, float]:
    x = i / n
    y = j / n
    if warp == 0.0 or i in (0, n) or j in (0, n):
        return x, y
    x += warp / n * math.sin(math.pi * x) * math.sin(2.0 * math.pi * y)
    y += warp / n * math.sin(2.0 * math.pi * x) * math.sin(math.pi * y)
    return x, y


def write_mesh(path: Path, n: int, warp: float, solver: str) -> dict[str, int | float | str]:
    width = n + 1
    plane = width * width
    triangle_count = 2 * n * n
    # Only one 2-D layer of edges is resident, independent of extrusion depth.
    edges: dict[tuple[int, int], list[int]] = {}
    for cell, triangle in enumerate(grid_triangles(n)):
        for index in range(3):
            a = triangle[index]
            b = triangle[(index + 1) % 3]
            key = (min(a, b), max(a, b))
            record = edges.get(key)
            if record is None:
                edges[key] = [a, b, cell, -1]
            elif record[3] == -1:
                record[3] = cell
            else:
                raise RuntimeError(f"nonmanifold triangulation edge {key}")

    # Face count: extruded 2-D edges plus horizontal interfaces and two caps.
    horizontal_faces = triangle_count * (n + 1)
    face_count = len(edges) * n + horizontal_faces
    vertex_count = plane * (n + 1)
    with path.open("w", encoding="ascii", buffering=4 * 1024 * 1024) as out:
        out.write(f"BABELSIM_MESH 3\nvertices\n{vertex_count}\n")
        for k in range(n + 1):
            z = k / n
            for j in range(n + 1):
                for i in range(n + 1):
                    x, y = warped_point(i, j, n, warp)
                    out.write(f"{x:.12g} {y:.12g} {z:.12g}\n")

        out.write(f"faces\n{face_count}\n")
        for k in range(n):
            offset = k * plane
            top_offset = offset + plane
            for a, b, owner, neighbour in edges.values():
                ring = (a + offset, b + offset, b + top_offset, a + top_offset)
                if neighbour >= 0:
                    patch = -1
                    global_neighbour = k * triangle_count + neighbour
                else:
                    ai, aj = a % width, a // width
                    bi, bj = b % width, b // width
                    if ai == bi == 0:
                        patch = 0  # hot / x-min
                    elif ai == bi == n:
                        patch = 1  # cold / x-max
                    elif aj == bj == 0:
                        patch = 2  # lower / y-min
                    elif aj == bj == n:
                        patch = 3  # upper / y-max
                    else:
                        raise RuntimeError(f"unassigned boundary edge {(a, b)}")
                    global_neighbour = -1
                global_owner = k * triangle_count + owner
                out.write("face 4 {} {} {} {} {} {} {}\n".format(
                    global_owner, global_neighbour, patch, *ring))

        # Horizontal triangle interfaces.  A CCW ring points out of the lower
        # prism toward +z; the same oriented face connects it to the upper prism.
        for triangle_id, triangle in enumerate(grid_triangles(n)):
            reversed_bottom = (triangle[0], triangle[2], triangle[1])
            ring = tuple(v for v in reversed_bottom)
            out.write("face 3 {} -1 4 {} {} {}\n".format(
                triangle_id, *ring))  # front / z-min
        for k in range(1, n):
            offset = k * plane
            for triangle_id, triangle in enumerate(grid_triangles(n)):
                owner = (k - 1) * triangle_count + triangle_id
                neighbour = k * triangle_count + triangle_id
                ring = tuple(v + offset for v in triangle)
                out.write("face 3 {} {} -1 {} {} {}\n".format(
                    owner, neighbour, *ring))
        offset = n * plane
        for triangle_id, triangle in enumerate(grid_triangles(n)):
            owner = (n - 1) * triangle_count + triangle_id
            ring = tuple(v + offset for v in triangle)
            out.write("face 3 {} -1 5 {} {} {}\n".format(
                owner, *ring))  # back / z-max

        if solver == "heat":
            patches = (
                ("hot", "wall"), ("cold", "wall"), ("lower", "symmetry"),
                ("upper", "symmetry"), ("front", "symmetry"), ("back", "symmetry"),
            )
        else:
            patches = (
                ("inlet", "inlet"), ("outlet", "outlet"),
                ("lower_wall", "wall"), ("upper_wall", "wall"),
                ("front", "symmetry"), ("back", "symmetry"),
            )
        out.write(f"patches\n{len(patches)}\n")
        for name, kind in patches:
            out.write(f"patch {name} {kind}\n")
        out.write("end\n")
    return {
        "gridN": n,
        "cells": triangle_count * n,
        "vertices": vertex_count,
        "faces": face_count,
        "warpFraction": warp,
    }


def write_case(case: Path, mesh: Path) -> None:
    case.mkdir(parents=True, exist_ok=True)
    for name in ("fields", "physics", "numerics"):
        source = Path(__file__).resolve().parents[1] / "cases" / "heat" / name
        target = case / name
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(source, target)
    mesh_target = case / "mesh.mesh"
    if mesh.resolve() != mesh_target.resolve():
        shutil.copy2(mesh, mesh_target)
    (case / "case.bs").write_text(
        "solver heat\nmesh mesh.mesh\nfields fields/initial\nphysics physics/thermal.bs\n"
        "methods numerics/methods.bs\nsolution numerics/solution.bs\ncontrol control.bs\n"
        "output output.bs\nghostLayers 3\n", encoding="utf-8")
    (case / "control.bs").write_text(
        "startTime 0\nendTime 1\ndeltaT 1\n", encoding="utf-8")
    (case / "output.bs").write_text(
        "directory results\ntimeName final\nwriteInterval 1\nwriteFields T\n",
        encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="study data directory")
    parser.add_argument("--grid-n", type=int, required=True,
                        help="same cell count in each direction; cells = 2*n^3")
    parser.add_argument("--warp", type=float, default=0.15,
                        help="smooth displacement as a fraction of planar grid spacing")
    parser.add_argument("--name", default=None)
    args = parser.parse_args()
    if args.grid_n < 2 or not 0.0 <= args.warp <= 0.25:
        parser.error("grid-n must be >= 2 and warp must be between 0 and 0.25")
    root = args.output.resolve()
    name = args.name or f"prism-n{args.grid_n}-w{args.warp:g}"
    case = root / name
    case.mkdir(parents=True, exist_ok=True)
    mesh = case / "mesh.mesh"
    metadata = write_mesh(mesh, args.grid_n, args.warp, "heat")
    write_case(case, mesh)
    metadata["case"] = str(case)
    metadata["meshBytes"] = mesh.stat().st_size
    (case / "generation.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
