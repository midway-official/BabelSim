import math
from pathlib import Path
from collections import defaultdict
SPAN_THICKNESS=0.02
def cross(a,b): return a[0]*b[1]-a[1]*b[0]
def signed_area(p): return 0.5*sum(cross(a,b) for a,b in zip(p,p[1:]+p[:1]))
def newell_area(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    nx = ny = nz = 0.0
    for a, b in zip(points, points[1:] + points[:1]):
        nx += (a[1] - b[1]) * (a[2] + b[2])
        ny += (a[2] - b[2]) * (a[0] + b[0])
        nz += (a[0] - b[0]) * (a[1] + b[1])
    return (0.5 * nx, 0.5 * ny, 0.5 * nz)


def write_volume_mesh(path: Path, vertices2d: list[tuple[float, float]],
                      cells2d: list[tuple[int, ...]],
                      boundary_edges: dict[tuple[int, int], str],
                      region_cell_ranges: dict[str, list[int]]) \
        -> tuple[list[tuple[float, float, float]], list[tuple[list[int], int, int, int]], dict[str, object]]:
    vertices = [(x, y, z) for z in (0.0, SPAN_THICKNESS)
                for x, y in vertices2d]
    n2 = len(vertices2d)
    patch_names = ["airfoil", "inlet", "outlet", "farfield", "front", "back"]
    patch_ids = {name: i for i, name in enumerate(patch_names)}
    faces: dict[tuple[int, ...], list[int]] = {}
    cell_centres: list[tuple[float, float, float]] = []

    def add_face(ring: list[int], owner: int, patch: int = -1) -> None:
        key = tuple(sorted(ring))
        if key not in faces:
            faces[key] = [ring, owner, -1, patch]
        else:
            record = faces[key]
            if record[2] != -1:
                raise ValueError("non-manifold volume face")
            record[2] = owner
            if record[3] != -1:
                raise ValueError("boundary face became internal")

    for cell_id, cell in enumerate(cells2d):
        n = len(cell)
        bottom = list(cell)
        top = [v + n2 for v in cell]
        points = [vertices[v] for v in bottom + top]
        polygon = [vertices2d[v] for v in cell]
        area = signed_area(polygon)
        centroid = [sum((a[k]+b[k])*cross(a,b) for a,b in zip(polygon,polygon[1:]+polygon[:1]))/(6*area) for k in range(2)]
        cell_centres.append((centroid[0],centroid[1],SPAN_THICKNESS/2))
        add_face(bottom, cell_id, patch_ids["front"])
        add_face(top, cell_id, patch_ids["back"])
        for i in range(n):
            a, b = cell[i], cell[(i + 1) % n]
            edge = tuple(sorted((a, b)))
            patch = patch_ids[boundary_edges[edge]] if edge in boundary_edges else -1
            add_face([a, b, b + n2, a + n2], cell_id, patch)

    records = [(record[0], record[1], record[2], record[3]) for record in faces.values()]
    if any(neighbour == -1 and patch == -1 for _, _, neighbour, patch in records):
        raise ValueError("unassigned boundary face in volume mesh")
    if any(neighbour != -1 and patch != -1 for _, _, neighbour, patch in records):
        raise ValueError("internal face has a boundary patch")

    output = ["BABELSIM_MESH 3", "vertices", str(len(vertices))]
    output.extend(f"{x:.15g} {y:.15g} {z:.15g}" for x, y, z in vertices)
    output.extend(("faces", str(len(records))))
    for ring, owner, neighbour, patch in records:
        output.append("face {} {} {} {} {}".format(
            len(ring), owner, neighbour, patch, " ".join(map(str, ring))))
    kinds = {"airfoil": "wall", "inlet": "inlet", "outlet": "outlet",
             "farfield": "generic", "front": "symmetry", "back": "symmetry"}
    output.extend(("patches", str(len(patch_names))))
    output.extend(f"patch {name} {kinds[name]}" for name in patch_names)
    output.append("end")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(output) + "\n", encoding="utf-8")

    face_counts = defaultdict(int)
    nonorthogonality = []
    skewness = []
    planar_areas = []
    planar_aspect = []
    planar_edges = []
    minimum_corner_angles = []
    concave_cells = 0
    concave_cell_flags = []
    for cell in cells2d:
        pts = [vertices2d[v] for v in cell]
        area = signed_area(pts)
        planar_areas.append(area)
        lengths = [math.dist(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
        planar_edges.extend(lengths)
        planar_aspect.append(max(lengths) / max(min(lengths), 1e-300))
        angles = []
        is_concave = False
        for i, current in enumerate(pts):
            previous, following = pts[i - 1], pts[(i + 1) % len(pts)]
            a = (previous[0] - current[0], previous[1] - current[1])
            b = (following[0] - current[0], following[1] - current[1])
            turn = cross((current[0] - previous[0], current[1] - previous[1]),
                         (following[0] - current[0], following[1] - current[1]))
            if turn < -1e-12:
                is_concave = True
            cosine = sum(a[k] * b[k] for k in range(2)) / max(
                math.hypot(*a) * math.hypot(*b), 1e-300)
            angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
        minimum_corner_angles.append(min(angles))
        concave_cells += int(is_concave)
        concave_cell_flags.append(is_concave)
    for ring, owner, neighbour, _patch in records:
        if neighbour < 0:
            continue
        pts = [vertices[v] for v in ring]
        area_vec = newell_area(pts)
        delta = tuple(cell_centres[neighbour][k] - cell_centres[owner][k]
                      for k in range(3))
        an = math.sqrt(sum(v*v for v in area_vec))
        dn = math.sqrt(sum(v*v for v in delta))
        cosine = abs(sum(area_vec[k] * delta[k] for k in range(3))) / max(an * dn, 1e-300)
        nonorthogonality.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
        fc = tuple(sum(p[k] for p in pts) / len(pts) for k in range(3))
        fraction = sum(area_vec[k]*(fc[k]-cell_centres[owner][k]) for k in range(3)) / sum(area_vec[k]*delta[k] for k in range(3))
        midpoint = tuple(cell_centres[owner][k]+fraction*delta[k] for k in range(3))
        skewness.append(math.sqrt(sum((fc[k] - midpoint[k])**2 for k in range(3)))
                        / max(dn, 1e-300))
    for _, _, _, patch in records:
        face_counts[patch] += 1

    cell_region = ["unclassified"] * len(cells2d)
    for name, (start, end) in region_cell_ranges.items():
        for cell_id in range(start, end):
            cell_region[cell_id] = name
    cell_region_quality = {}
    for name, (start, end) in region_cell_ranges.items():
        region_areas = planar_areas[start:end]
        region_aspect = planar_aspect[start:end]
        region_angles = minimum_corner_angles[start:end]
        region_concave = sum(concave_cell_flags[start:end])
        cell_region_quality[name] = {
            "cells": end - start,
            "min_volume": min(region_areas) * SPAN_THICKNESS,
            "max_volume": max(region_areas) * SPAN_THICKNESS,
            "planar_aspect_ratio_max": max(region_aspect),
            "planar_aspect_ratio_p95": percentile(region_aspect, 0.95),
            "minimum_corner_angle_deg": min(region_angles),
            "cells_with_concave_corners": region_concave,
        }
    internal_face_groups: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: {"nonorthogonality_deg": [], "skewness": []})
    for (_ring, owner, neighbour, _patch), nonorth, skew in zip(
            (record for record in records if record[2] >= 0),
            nonorthogonality, skewness):
        regions = sorted((cell_region[owner], cell_region[neighbour]))
        internal_face_groups[" <-> ".join(regions)]["nonorthogonality_deg"].append(nonorth)
        internal_face_groups[" <-> ".join(regions)]["skewness"].append(skew)
    internal_face_quality = {}
    for name, metrics in internal_face_groups.items():
        no = metrics["nonorthogonality_deg"]
        sk = metrics["skewness"]
        internal_face_quality[name] = {
            "faces": len(no),
            "nonorthogonality_deg_max": max(no),
            "nonorthogonality_deg_p95": percentile(no, 0.95),
            "nonorthogonality_over_20_deg": sum(v > 20.0 for v in no),
            "nonorthogonality_over_45_deg": sum(v > 45.0 for v in no),
            "nonorthogonality_over_70_deg": sum(v > 70.0 for v in no),
            "skewness_max": max(sk),
            "skewness_p95": percentile(sk, 0.95),
            "skewness_over_0_5": sum(v > 0.5 for v in sk),
        }
    quality = {
        "cell_count": len(cells2d),
        "vertex_count": len(vertices),
        "face_count": len(records),
        "cell_type": "extruded quadrilateral hexahedra and triangular prisms",
        "non_hexahedral_cells": sum(len(c) != 4 for c in cells2d),
        "extrusion_thickness": SPAN_THICKNESS,
        "positive_volume_cells": sum(area > 0.0 for area in planar_areas),
        "negative_or_zero_volume_cells": sum(area <= 0.0 for area in planar_areas),
        "cell_volume_min": min(planar_areas) * SPAN_THICKNESS,
        "cell_volume_max": max(planar_areas) * SPAN_THICKNESS,
        "planar_edge_length_min": min(planar_edges),
        "planar_edge_length_max": max(planar_edges),
        "planar_aspect_ratio_max": max(planar_aspect),
        "planar_aspect_ratio_p95": percentile(planar_aspect, 0.95),
        "minimum_cell_corner_angle_deg": min(minimum_corner_angles),
        "cells_with_minimum_corner_angle_below_5_deg": sum(v < 5.0 for v in minimum_corner_angles),
        "cells_with_minimum_corner_angle_below_10_deg": sum(v < 10.0 for v in minimum_corner_angles),
        "concave_cells": concave_cells,
        "cell_quality_by_region": cell_region_quality,
        "face_count_by_patch": {patch_names[k]: face_counts[k] for k in range(len(patch_names))},
        "internal_face_count": len(nonorthogonality),
        "internal_nonorthogonality_deg_max": max(nonorthogonality),
        "internal_nonorthogonality_deg_p95": percentile(nonorthogonality, 0.95),
        "internal_faces_nonorthogonality_over_20_deg": sum(v > 20.0 for v in nonorthogonality),
        "internal_faces_nonorthogonality_over_45_deg": sum(v > 45.0 for v in nonorthogonality),
        "internal_faces_nonorthogonality_over_70_deg": sum(v > 70.0 for v in nonorthogonality),
        "internal_face_skewness_max": max(skewness),
        "internal_face_skewness_p95": percentile(skewness, 0.95),
        "internal_faces_skewness_over_0_5": sum(v > 0.5 for v in skewness),
        "internal_face_quality_by_region_pair": internal_face_quality,
    }
    return vertices, records, quality


def percentile(values: list[float], p: float) -> float:
    sorted_values = sorted(values)
    return sorted_values[min(len(sorted_values) - 1, int(p * (len(sorted_values) - 1)))]
