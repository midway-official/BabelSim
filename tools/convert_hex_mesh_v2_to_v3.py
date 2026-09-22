#!/usr/bin/env python3
"""One-time migration of BabelSim's old fixed-hex mesh files to face topology.

The old format is intentionally supported only by this migration utility.  The
solver reader accepts version 3, whose records contain variable-length faces
and explicit owner/neighbour/patch labels.
"""

from __future__ import annotations

import argparse
from pathlib import Path


HEX_FACES = (
    (0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4),
    (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7),
)


def read_v2(path: Path):
    tokens = path.read_text().split()
    cursor = 0

    def take(expected: str | None = None) -> str:
        nonlocal cursor
        if cursor >= len(tokens):
            raise ValueError(f"{path}: unexpected end of file")
        value = tokens[cursor]
        cursor += 1
        if expected is not None and value != expected:
            raise ValueError(f"{path}: expected {expected}, got {value}")
        return value

    if take() != "BABELSIM_MESH" or take() != "2":
        raise ValueError(f"{path}: not a BABELSIM_MESH version 2 file")
    take("vertices")
    vertex_count = int(take())
    vertices = [(take(), take(), take()) for _ in range(vertex_count)]
    take("cells")
    cell_count = int(take())
    cells = [tuple(int(take()) for _ in range(8)) for _ in range(cell_count)]
    take("patches")
    patch_count = int(take())
    patches = []
    boundary = {}
    for patch in range(patch_count):
        take("patch")
        name = take()
        kind = take()
        patches.append((name, kind))
        face_count = int(take())
        for _ in range(face_count):
            ring = tuple(int(take()) for _ in range(4))
            key = tuple(sorted(ring))
            if key in boundary:
                raise ValueError(f"{path}: duplicate boundary face {key}")
            boundary[key] = patch
    if take() != "end" or cursor != len(tokens):
        raise ValueError(f"{path}: malformed trailing data")
    return vertices, cells, patches, boundary


def convert(path: Path) -> None:
    vertices, cells, patches, boundary = read_v2(path)
    faces = {}
    for owner, cell in enumerate(cells):
        for local_face in HEX_FACES:
            ring = tuple(cell[index] for index in local_face)
            key = tuple(sorted(ring))
            if key not in faces:
                faces[key] = [list(ring), owner, -1, -1]
            elif faces[key][2] != -1:
                raise ValueError(f"{path}: non-manifold face {key}")
            else:
                faces[key][2] = owner
    records = []
    for key, (ring, owner, neighbour, patch) in faces.items():
        if neighbour == -1:
            if key not in boundary:
                raise ValueError(f"{path}: boundary face has no patch {key}")
            patch = boundary[key]
        elif key in boundary:
            raise ValueError(f"{path}: internal face has boundary patch {key}")
        records.append((ring, owner, neighbour, patch))

    output = ["BABELSIM_MESH 3", "vertices", str(len(vertices))]
    output.extend("{} {} {}".format(*vertex) for vertex in vertices)
    output.extend(("faces", str(len(records))))
    for ring, owner, neighbour, patch in records:
        output.append("face {} {} {} {} {}".format(
            len(ring), owner, neighbour, patch, " ".join(map(str, ring))))
    output.extend(("patches", str(len(patches))))
    for name, kind in patches:
        output.extend(("patch", name, kind))
    output.append("end")
    path.write_text("\n".join(output) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", type=Path)
    args = parser.parse_args()
    for path in args.paths:
        convert(path)


if __name__ == "__main__":
    main()
