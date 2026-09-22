#!/usr/bin/env python3
"""Generate a z-extruded NACA0012 O-grid in BabelSim mesh format."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from polyhedral_mesh import write_hex_v3


def naca0012_thickness(x: float) -> float:
    # Finite trailing-edge gap avoids a degenerate wrap face.
    return 5.0 * 0.12 * (0.2969 * math.sqrt(max(x, 0.0)) - 0.1260 * x
                         - 0.3516 * x * x + 0.2843 * x**3 - 0.0950 * x**4)


def surface_points(n_half: int, angle_deg: float = 15.0,
                   pivot_x: float = 0.25) -> list[tuple[float, float]]:
    # Uniform spacing avoids the very short cells at the leading/trailing
    # edges created by cosine clustering.  The radial offset supplies the
    # boundary-layer resolution while keeping the O-grid well conditioned.
    x = [i / n_half for i in range(n_half + 1)]
    lower = [(x[i], -naca0012_thickness(x[i])) for i in range(n_half, -1, -1)]
    upper = [(x[i], naca0012_thickness(x[i])) for i in range(1, n_half + 1)]
    angle = math.radians(angle_deg)
    cosine, sine = math.cos(angle), math.sin(angle)
    return [((point[0] - pivot_x) * cosine - point[1] * sine + pivot_x,
             (point[0] - pivot_x) * sine + point[1] * cosine)
            for point in lower + upper]


def _line(a: tuple[float, float], b: tuple[float, float], count: int,
          fractions: list[float] | None = None) -> list[tuple[float, float]]:
    if fractions is None:
        fractions = [i / count for i in range(count + 1)]
    return [(a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]))
            for f in fractions]


def _fractions(count: int, growth: float) -> list[float]:
    """Fractions clustered at the first point (the wall or trailing edge)."""
    weights = [growth**i for i in range(count)]
    total = sum(weights)
    running = 0.0
    result = [0.0]
    for weight in weights:
        running += weight
        result.append(running / total)
    return result


def _coons(bottom, top, left, right, ni: int, nj: int):
    """Generate a logically rectangular 2-D block from four boundaries."""
    result = []
    for j in range(nj + 1):
        v = j / nj
        for i in range(ni + 1):
            u = i / ni
            # Bilinear Coons patch.  All four boundary arrays are supplied
            # explicitly so neighbouring blocks share the same coordinates.
            x = ((1.0 - v) * bottom[i][0] + v * top[i][0]
                 + (1.0 - u) * left[j][0] + u * right[j][0]
                 - ((1.0 - u) * (1.0 - v) * left[0][0]
                    + u * (1.0 - v) * right[0][0]
                    + (1.0 - u) * v * left[-1][0]
                    + u * v * right[-1][0]))
            y = ((1.0 - v) * bottom[i][1] + v * top[i][1]
                 + (1.0 - u) * left[j][1] + u * right[j][1]
                 - ((1.0 - u) * (1.0 - v) * left[0][1]
                    + u * (1.0 - v) * right[0][1]
                    + (1.0 - u) * v * left[-1][1]
                    + u * v * right[-1][1]))
            result.append((x, y))
    return result


def geometry(n_half: int, n_radial: int, growth: float, thickness: float,
             n_z: int, angle_deg: float, pivot_x: float,
             x_min: float, x_max: float, y_min: float, y_max: float,
             n_left: int = 70, n_wake: int = 320,
             n_wake_normal: int = 20):
    """Build a rectangular-domain C-grid with a separately refined wake.

    The upper/lower airfoil blocks provide the wall-normal boundary-layer
    spacing.  Two wake blocks and a thin middle wake block continue from the
    finite trailing edge to the outlet, so no single O-grid closure edge is
    stretched across the 24-chord-height rectangle.
    """
    surface = surface_points(n_half, angle_deg, pivot_x)
    half = n_half
    nose = surface[half]
    upper = surface[half:2 * half + 1]
    lower = list(reversed(surface[:half + 1]))
    upper_te, lower_te = upper[-1], lower[-1]
    wall_fractions = _fractions(n_radial, growth)
    wake_fractions = _fractions(n_wake, 1.012)

    # The wake is centred on the trailing-edge height and expands to a
    # controllable strip before the outlet.  The outer boundary remains the
    # exact rectangle x=[x_min,x_max], y=[y_min,y_max].
    wake_lower, wake_upper = -2.0, 2.0
    xnose, xte = nose[0], 0.5 * (upper_te[0] + lower_te[0])
    nj = n_radial
    nl = n_left
    nw = n_wake
    nm = n_wake_normal
    vertices: list[tuple[float, float, float]] = []
    vertex_ids: dict[tuple[float, float, int], int] = {}
    cells = []
    patches = {name: [] for name in
               ("airfoil", "inlet", "outlet", "farfield", "front", "back")}

    def vertex(point, plane):
        key = (round(point[0], 12), round(point[1], 12), plane)
        if key not in vertex_ids:
            vertex_ids[key] = len(vertices)
            vertices.append((point[0], point[1], thickness * plane / n_z))
        return vertex_ids[key]

    def add_block(bottom, top, left, right, ni, nj_block, boundary):
        points = _coons(bottom, top, left, right, ni, nj_block)
        ids = []
        for j in range(nj_block + 1):
            row = []
            for i in range(ni + 1):
                point = points[j * (ni + 1) + i]
                row.append(tuple(vertex(point, plane) for plane in range(n_z + 1)))
            ids.append(row)
        for j in range(nj_block):
            for i in range(ni):
                a, b = ids[j][i], ids[j][i + 1]
                c, d = ids[j + 1][i + 1], ids[j + 1][i]
                for plane in range(n_z):
                    cells.append((a[plane], b[plane], c[plane], d[plane],
                                  a[plane + 1], b[plane + 1],
                                  c[plane + 1], d[plane + 1]))

        def add_faces(side, name):
            if name is None:
                return
            if side == "bottom":
                edges = [(ids[0][i], ids[0][i + 1]) for i in range(ni)]
            elif side == "top":
                edges = [(ids[-1][i], ids[-1][i + 1]) for i in range(ni)]
            elif side == "left":
                edges = [(ids[j][0], ids[j + 1][0]) for j in range(nj_block)]
            else:
                edges = [(ids[j][-1], ids[j + 1][-1]) for j in range(nj_block)]
            for a, b in edges:
                for plane in range(n_z):
                    patches[name].append((a[plane], b[plane],
                                          b[plane + 1], a[plane + 1]))

        for side, name in boundary.items():
            add_faces(side, name)

    # Shared boundaries are passed by reference to the same line builders;
    # this makes all block interfaces conformal after vertex deduplication.
    nose_top = _line(nose, (xnose, y_max), nj,
                     _fractions(nj, 1.06))
    nose_bottom = _line((xnose, y_min), nose, nj,
                        _fractions(nj, 1.06))
    ute_top = _line(upper_te, (upper_te[0], y_max), nj,
                    _fractions(nj, 1.06))
    lte_bottom = _line((lower_te[0], y_min), lower_te, nj,
                       _fractions(nj, 1.06))

    # Left upstream blocks.
    add_block(_line((x_min, nose[1]), nose, nl),
              _line((x_min, y_max), (xnose, y_max), nl),
              _line((x_min, nose[1]), (x_min, y_max), nj), nose_top,
              nl, nj, {"left": "inlet", "top": "farfield"})
    add_block(_line((x_min, y_min), (xnose, y_min), nl),
              _line((x_min, nose[1]), nose, nl),
              _line((x_min, y_min), (x_min, nose[1]), nj), nose_bottom,
              nl, nj, {"left": "inlet", "bottom": "farfield"})

    # Body-fitted upper and lower blocks.  The wall-normal fractions are
    # geometric only in this direction; streamwise wall spacing is uniform.
    add_block(upper, _line((xnose, y_max), (upper_te[0], y_max), half),
              nose_top, ute_top, half, nj,
              {"bottom": "airfoil", "top": "farfield"})
    add_block(_line((xnose, y_min), (lower_te[0], y_min), half), lower,
              nose_bottom, lte_bottom, half, nj,
              {"top": "airfoil", "bottom": "farfield"})

    # Streamwise wake lines use a mild geometric expansion from the trailing
    # edge, providing substantially more cells in the near wake than in the
    # far downstream rectangle.
    wake_top = _line(upper_te, (x_max, wake_upper), nw, wake_fractions)
    wake_bottom = _line(lower_te, (x_max, wake_lower), nw, wake_fractions)
    upper_top = _line((upper_te[0], y_max), (x_max, y_max), nw, wake_fractions)
    lower_bottom = _line((lower_te[0], y_min), (x_max, y_min), nw, wake_fractions)
    upper_wake_left = ute_top
    lower_wake_left = lte_bottom
    upper_wake_right = _line((x_max, wake_upper), (x_max, y_max), nj,
                             _fractions(nj, 1.06))
    lower_wake_right = _line((x_max, y_min), (x_max, wake_lower), nj,
                             _fractions(nj, 1.06))
    add_block(wake_top, upper_top, upper_wake_left, upper_wake_right,
              nw, nj, {"top": "farfield"})
    add_block(lower_bottom, wake_bottom, lower_wake_left, lower_wake_right,
              nw, nj, {"bottom": "farfield"})

    # Thin middle wake block closes the finite trailing-edge gap and expands
    # to the refined wake strip.  Its left edge is part of the airfoil wall.
    middle_right = _line((x_max, wake_lower), (x_max, wake_upper), nm)
    add_block(wake_bottom, wake_top, _line(lower_te, upper_te, nm), middle_right,
              nw, nm, {"left": "airfoil"})

    # Complete the external rectangular boundary: the upper/lower wake outlet
    # segments and the middle strip together form x=x_max exactly.
    def append_line_patch(points, name):
        for a, b in zip(points, points[1:]):
            for plane in range(n_z):
                patches[name].append((vertex(a, plane), vertex(b, plane),
                                      vertex(b, plane + 1), vertex(a, plane + 1)))

    append_line_patch(upper_wake_right, "outlet")
    append_line_patch(middle_right, "outlet")
    append_line_patch(lower_wake_right, "outlet")
    # Front/back are the two symmetry planes of the extrusion.  They are
    # collected from every block face by identifying the z=0 and z=thickness
    # quads after all cells have been created.
    face_slots = ((0, 1, 2, 3), (4, 5, 6, 7))
    seen = {}
    for cell in cells:
        for slot, name in zip(face_slots, ("front", "back")):
            face = tuple(cell[k] for k in slot)
            key = tuple(sorted(face))
            entry = seen.setdefault((name, key), [0, face])
            entry[0] += 1
    # Only exterior z faces occur once; interfaces occur twice and are removed.
    for name in ("front", "back"):
        patches[name].extend(entry[1] for (which, _key), entry in seen.items()
                             if which == name and entry[0] == 1)
    outer = [(x_max, y_min), (x_min, y_min), (x_min, y_max), (x_max, y_max)]
    return surface, outer, wall_fractions, vertices, cells, patches


def signed_hex_volume(cell, vertices) -> float:
    faces = ((0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4),
             (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7))
    result = 0.0
    for a, b, c, d in faces:
        va, vb, vc, vd = (vertices[cell[k]] for k in (a, b, c, d))
        def triple(x, y, z):
            return (x[0] * (y[1] * z[2] - y[2] * z[1])
                    + x[1] * (y[2] * z[0] - y[0] * z[2])
                    + x[2] * (y[0] * z[1] - y[1] * z[0]))
        result += (triple(va, vb, vc) + triple(va, vc, vd)) / 6.0
    return result


def vsub(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def vnorm(a):
    return math.sqrt(sum(value * value for value in a))


def vdot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def quality(surface, outer, vertices, cells, n_z, patches):
    volumes = [signed_hex_volume(cell, vertices) for cell in cells]
    planar_edges = ((0, 1), (1, 2), (2, 3), (3, 0))
    all_edges = planar_edges + ((4, 5), (5, 6), (6, 7), (7, 4),
                                (0, 4), (1, 5), (2, 6), (3, 7))
    centres = [tuple(sum(vertices[v][q] for v in cell) / 8.0 for q in range(3))
               for cell in cells]
    planar_lengths = [vnorm(vsub(vertices[cell[a]], vertices[cell[b]]))
                      for cell in cells for a, b in planar_edges]
    all_lengths = [vnorm(vsub(vertices[cell[a]], vertices[cell[b]]))
                   for cell in cells for a, b in all_edges]
    planar_aspect = []
    for cell in cells:
        lengths = [vnorm(vsub(vertices[cell[a]], vertices[cell[b]]))
                   for a, b in planar_edges]
        planar_aspect.append(max(lengths) / max(min(lengths), 1.0e-300))

    face_slots = ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
                  (3, 2, 6, 7), (0, 4, 7, 3), (1, 2, 6, 5))
    faces = {}
    for cell_id, cell in enumerate(cells):
        for slot in face_slots:
            face = tuple(cell[k] for k in slot)
            faces.setdefault(tuple(sorted(face)), []).append((cell_id, face))
    nonorth, skew = [], []
    for entries in faces.values():
        if len(entries) != 2:
            continue
        owner, face = entries[0]
        neighbour = entries[1][0]
        points = [vertices[v] for v in face]
        a, b, c, d = points
        ab, ac, ad = vsub(b, a), vsub(c, a), vsub(d, a)
        area = (ab[1] * ac[2] - ab[2] * ac[1] + ac[1] * ad[2] - ac[2] * ad[1],
                ab[2] * ac[0] - ab[0] * ac[2] + ac[2] * ad[0] - ac[0] * ad[2],
                ab[0] * ac[1] - ab[1] * ac[0] + ac[0] * ad[1] - ac[1] * ad[0])
        direction = vsub(centres[neighbour], centres[owner])
        cosine = abs(vdot(area, direction)) / max(vnorm(area) * vnorm(direction), 1.0e-300)
        nonorth.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
        fc = tuple(sum(point[q] for point in points) / 4.0 for q in range(3))
        midpoint = tuple(0.5 * (centres[owner][q] + centres[neighbour][q]) for q in range(3))
        skew.append(vnorm(vsub(fc, midpoint)) / max(vnorm(direction), 1.0e-300))

    wall_keys = {tuple(sorted(face)) for face in patches.get("airfoil", [])}
    first_layer = []
    for key in wall_keys:
        entries = faces.get(key, [])
        if len(entries) != 1:
            continue
        cell_id, face = entries[0]
        fc = tuple(sum(vertices[v][q] for v in face) / 4.0 for q in range(3))
        first_layer.append(2.0 * vnorm(vsub(centres[cell_id], fc)))

    def percentile(values, fraction):
        values = sorted(values)
        return values[min(len(values) - 1, int(fraction * (len(values) - 1)))] if values else 0.0

    return {
        "cell_count": len(cells), "vertex_count": len(vertices),
        "wall_surface_points": len(surface), "spanwise_cells": n_z,
        "positive_volume_cells": sum(value > 0.0 for value in volumes),
        "negative_or_zero_volume_cells": sum(value <= 0.0 for value in volumes),
        "cell_volume_min": min(volumes), "cell_volume_max": max(volumes),
        "edge_length_min": min(planar_lengths), "edge_length_max": max(planar_lengths),
        "edge_length_3d_min": min(all_lengths), "edge_length_3d_max": max(all_lengths),
        "planar_aspect_ratio_max": max(planar_aspect),
        "planar_aspect_ratio_p95": percentile(planar_aspect, 0.95),
        "nonorthogonality_deg_max": max(nonorth) if nonorth else 0.0,
        "nonorthogonality_deg_p95": percentile(nonorth, 0.95),
        "face_skewness_max": max(skew) if skew else 0.0,
        "face_skewness_p95": percentile(skew, 0.95),
        "first_layer_spacing_min": min(first_layer) if first_layer else 0.0,
        "first_layer_spacing_max": max(first_layer) if first_layer else 0.0,
        "outer_x_range": [min(point[0] for point in outer), max(point[0] for point in outer)],
        "outer_y_range": [min(point[1] for point in outer), max(point[1] for point in outer)],
        "extrusion_thickness": max(point[2] for point in vertices) - min(point[2] for point in vertices),
    }


def write_mesh(path: Path, vertices, cells, patches):
    kind = {"airfoil": "wall", "inlet": "inlet", "outlet": "outlet", "farfield": "generic",
            "front": "symmetry", "back": "symmetry"}
    write_hex_v3(path, vertices, cells,
                 [(name, kind[name], faces) for name, faces in patches.items()])


def write_surface_map(path: Path, surface, n_z: int):
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(["segment", "owner_global_id", "x", "y", "nx", "ny", "length"])
        ns = len(surface)
        for i in range(ns):
            ip = (i + 1) % ns
            x0, y0, x1, y1 = *surface[i], *surface[ip]
            dx, dy = x1 - x0, y1 - y0
            length = math.hypot(dx, dy)
            writer.writerow([i, i * n_z, 0.5 * (x0 + x1), 0.5 * (y0 + y1),
                             -dy / length, dx / length, length])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "mesh" / "naca0012.mesh")
    parser.add_argument("--n-half", type=int, default=200)
    parser.add_argument("--n-radial", type=int, default=80,
                        help="wall-normal cells in the body and outer blocks")
    parser.add_argument("--n-left", type=int, default=140,
                        help="streamwise cells between inlet and leading edge")
    parser.add_argument("--n-wake", type=int, default=320,
                        help="streamwise cells from trailing edge to outlet")
    parser.add_argument("--n-wake-normal", type=int, default=20,
                        help="cross-stream cells in the middle wake strip")
    parser.add_argument("--growth", type=float, default=1.06)
    parser.add_argument("--thickness", type=float, default=1.0)
    parser.add_argument("--angle-deg", type=float, default=15.0,
                        help="airfoil angle of attack relative to horizontal flow")
    parser.add_argument("--pivot-x", type=float, default=0.25,
                        help="chordwise rotation pivot")
    parser.add_argument("--x-min", type=float, default=-3.0)
    parser.add_argument("--x-max", type=float, default=15.0)
    parser.add_argument("--y-min", type=float, default=-6.0)
    parser.add_argument("--y-max", type=float, default=6.0)
    parser.add_argument("--n-z", type=int, default=1,
                        help="spanwise extrusion cells; symmetry keeps the solution two-dimensional")
    args = parser.parse_args()
    if (args.n_half < 4 or args.n_radial < 2 or args.n_z < 1
            or args.growth <= 1.0 or args.x_max <= args.x_min
            or args.y_max <= args.y_min or args.n_left < 2 or args.n_wake < 4
            or args.n_wake_normal < 2):
        parser.error("invalid grid dimensions or growth")
    surface, _outer, _fractions, vertices, cells, patches = geometry(
        args.n_half, args.n_radial, args.growth, args.thickness, args.n_z,
        args.angle_deg, args.pivot_x, args.x_min, args.x_max, args.y_min,
        args.y_max, args.n_left, args.n_wake, args.n_wake_normal)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_mesh(args.output, vertices, cells, patches)
    data = quality(surface, _outer, vertices, cells, args.n_z, patches)
    data.update({
        "angle_of_attack_deg": args.angle_deg,
        "rotation_pivot_x": args.pivot_x,
        "farfield_box": {"x_min": args.x_min, "x_max": args.x_max,
                          "y_min": args.y_min, "y_max": args.y_max},
        "farfield_geometry": "exact_rectangle_c_grid",
        "block_counts": {"left": args.n_left, "wall_normal": args.n_radial,
                         "wake": args.n_wake,
                         "wake_normal": args.n_wake_normal},
    })
    args.output.with_name("mesh_quality.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    write_surface_map(args.output.with_name("surface_faces.csv"), surface, args.n_z)
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
