"""Writers for BabelSim's face-based BABELSIM_MESH 3 format.

The hex helper is intentionally a producer-side convenience for migrating
structured case generators.  The solver and Mesh API consume only explicit
face records; no fixed-hex reader or constructor is involved at runtime.
"""

from __future__ import annotations

from pathlib import Path


HEX_FACES = (
    (0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4),
    (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7),
)


def face_records(cells, patch_faces):
    """Return ``(ring, owner, neighbour, patch)`` records from hex producer data."""
    faces = {}
    for owner, cell in enumerate(cells):
        for local_face in HEX_FACES:
            ring = tuple(cell[index] for index in local_face)
            key = tuple(sorted(ring))
            if key not in faces:
                faces[key] = [list(ring), owner, -1, -1]
            elif faces[key][2] != -1:
                raise ValueError(f"non-manifold face {key}")
            else:
                faces[key][2] = owner
    boundary = {}
    for patch, rings in enumerate(patch_faces):
        for ring in rings:
            key = tuple(sorted(ring))
            if key in boundary:
                raise ValueError(f"duplicate boundary face {key}")
            boundary[key] = patch
    records = []
    for key, (ring, owner, neighbour, patch) in faces.items():
        if neighbour == -1:
            if key not in boundary:
                raise ValueError(f"boundary face has no patch {key}")
            patch = boundary[key]
        elif key in boundary:
            raise ValueError(f"internal face has boundary patch {key}")
        records.append((ring, owner, neighbour, patch))
    return records


def write_v3(path: Path, vertices, records, patches) -> None:
    output = ["BABELSIM_MESH 3", "vertices", str(len(vertices))]
    output.extend("{} {} {}".format(*vertex) for vertex in vertices)
    output.extend(("faces", str(len(records))))
    for ring, owner, neighbour, patch in records:
        output.append("face {} {} {} {} {}".format(
            len(ring), owner, neighbour, patch, " ".join(map(str, ring))))
    output.extend(("patches", str(len(patches))))
    for name, kind in patches:
        output.append(f"patch {name} {kind}")
    output.append("end")
    path.write_text("\n".join(output) + "\n", encoding="utf-8")


def write_hex_v3(path: Path, vertices, cells, patches) -> None:
    """Write a v3 file from legacy generator-side hex arrays."""
    records = face_records(cells, [rings for _, _, rings in patches])
    write_v3(path, vertices, records, [(name, kind) for name, kind, _ in patches])
