#!/usr/bin/env python3
"""Generate an explicit face-based BABELSIM_MESH 3 mesh file."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from polyhedral_mesh import write_hex_v3


FACES = (
    ("minus_x", (0, 4, 7, 3)),
    ("plus_x", (1, 2, 6, 5)),
    ("minus_y", (0, 1, 5, 4)),
    ("plus_y", (3, 7, 6, 2)),
    ("minus_z", (0, 3, 2, 1)),
    ("plus_z", (4, 5, 6, 7)),
)


def coordinate(index: int, count: int, clustering: float) -> float:
    if clustering == 0.0:
        return index / count
    scale = math.tanh(clustering)
    return 0.5 * (1.0 + math.tanh(clustering * (2.0 * index / count - 1.0)) / scale)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cells", type=int, nargs=3, metavar=("A", "B", "C"), required=True)
    parser.add_argument("--bounds", type=float, nargs=6, metavar=("X0", "Y0", "Z0", "X1", "Y1", "Z1"),
                        required=True)
    parser.add_argument("--tanh-clustering", type=float, default=0.0)
    parser.add_argument("--patch", nargs=3, action="append", metavar=("FACE", "NAME", "KIND"), default=[])
    args = parser.parse_args()
    a_count, b_count, c_count = args.cells
    if min(args.cells) <= 0 or args.tanh_clustering < 0.0:
        parser.error("cell counts must be positive and tanh clustering must be non-negative")
    x0, y0, z0, x1, y1, z1 = args.bounds
    if not (x1 > x0 and y1 > y0 and z1 > z0):
        parser.error("bounds must increase on every axis")
    patch_values = {face: (face, "generic") for face, _ in FACES}
    for face, name, kind in args.patch:
        if face not in patch_values:
            parser.error(f"unknown face {face}")
        if kind not in {"generic", "wall", "inlet", "outlet", "symmetry"}:
            parser.error(f"unknown patch kind {kind}")
        patch_values[face] = (name, kind)

    def vertex(a: int, b: int, c: int) -> int:
        return a + (a_count + 1) * (b + (b_count + 1) * c)

    vertices: list[tuple[float, float, float]] = []
    for c in range(c_count + 1):
        z = z0 + (z1 - z0) * coordinate(c, c_count, args.tanh_clustering)
        for b in range(b_count + 1):
            y = y0 + (y1 - y0) * coordinate(b, b_count, args.tanh_clustering)
            for a in range(a_count + 1):
                x = x0 + (x1 - x0) * coordinate(a, a_count, args.tanh_clustering)
                vertices.append((x, y, z))
    cells: list[tuple[int, ...]] = []
    faces: dict[str, list[tuple[int, int, int, int]]] = {name: [] for name, _ in FACES}
    for c in range(c_count):
        for b in range(b_count):
            for a in range(a_count):
                cell = (
                    vertex(a, b, c), vertex(a + 1, b, c), vertex(a + 1, b + 1, c), vertex(a, b + 1, c),
                    vertex(a, b, c + 1), vertex(a + 1, b, c + 1),
                    vertex(a + 1, b + 1, c + 1), vertex(a, b + 1, c + 1),
                )
                cells.append(cell)
                if a == 0:
                    faces["minus_x"].append(tuple(cell[index] for index in FACES[0][1]))
                if a + 1 == a_count:
                    faces["plus_x"].append(tuple(cell[index] for index in FACES[1][1]))
                if b == 0:
                    faces["minus_y"].append(tuple(cell[index] for index in FACES[2][1]))
                if b + 1 == b_count:
                    faces["plus_y"].append(tuple(cell[index] for index in FACES[3][1]))
                if c == 0:
                    faces["minus_z"].append(tuple(cell[index] for index in FACES[4][1]))
                if c + 1 == c_count:
                    faces["plus_z"].append(tuple(cell[index] for index in FACES[5][1]))
    patches = [(patch_values[face][0], patch_values[face][1], faces[face])
               for face, _ in FACES]
    write_hex_v3(args.output, vertices, cells, patches)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
