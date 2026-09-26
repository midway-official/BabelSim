#!/usr/bin/env python3
"""Build the orthogonal 2-D planar jet hex mesh and its seeded initial velocity."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

CASE = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE.parents[1] / "tools"))
from polyhedral_mesh import write_hex_v3  # noqa: E402

NX, NY, NZ = 500, 200, 1
X0, X1 = 0.0, 20.0
Y0, Y1 = -5.0, 5.0
Z0, Z1 = 0.0, 0.05
D = 1.0
DX = (X1 - X0) / NX
DY = (Y1 - Y0) / NY
DZ = (Z1 - Z0) / NZ
RAMP_START, RAMP_END, RAMP_WIDTH = 0.25, 0.75, 0.05
RAMP_BANDS = int(round((RAMP_END - RAMP_START) / RAMP_WIDTH))
JET_WIDTH = 1.0
EDGE_THICKNESS = 0.06


def jet_speed(abs_y: float) -> float:
    return 0.5 * (1.0 + math.tanh((0.5 - abs_y) / EDGE_THICKNESS))


def main() -> None:
    xs = [X0 + i * DX for i in range(NX + 1)]
    ys = [Y0 + j * DY for j in range(NY + 1)]
    zs = [Z0 + k * DZ for k in range(NZ + 1)]

    def vertex(i: int, j: int, k: int) -> int:
        return i + (NX + 1) * (j + (NY + 1) * k)

    vertices = [(x, y, z) for z in zs for y in ys for x in xs]
    cells: list[tuple[int, ...]] = []
    inlet = {"jet_core": [], "coflow_inlet": []}
    inlet.update({f"inlet_ramp_{r:02d}": [] for r in range(RAMP_BANDS)})
    outlet, farfield, front_back = [], [], []
    initial_path = CASE / "fields/initial/U.initial.dat"
    initial_path.parent.mkdir(parents=True, exist_ok=True)
    inlet_mass_flow = 0.0

    with initial_path.open("w", encoding="utf-8") as initial:
        initial.write("# globalCellId Ux Uy Uz\n")
        for k in range(NZ):
            for j in range(NY):
                yc = 0.5 * (ys[j] + ys[j + 1])
                abs_y = abs(yc)
                ux = jet_speed(abs_y)
                inlet_mass_flow += ux * DY * DZ
                for i in range(NX):
                    xc = 0.5 * (xs[i] + xs[i + 1])
                    cell = (
                        vertex(i, j, k), vertex(i + 1, j, k),
                        vertex(i + 1, j + 1, k), vertex(i, j + 1, k),
                        vertex(i, j, k + 1), vertex(i + 1, j, k + 1),
                        vertex(i + 1, j + 1, k + 1), vertex(i, j + 1, k + 1),
                    )
                    cells.append(cell)
                    if i == 0:
                        face = (cell[0], cell[4], cell[7], cell[3])
                        if abs_y < RAMP_START:
                            inlet["jet_core"].append(face)
                        elif abs_y < RAMP_END:
                            band = int(math.floor((abs_y - RAMP_START - 1.0e-12) / RAMP_WIDTH))
                            band = max(0, min(RAMP_BANDS - 1, band))
                            inlet[f"inlet_ramp_{band:02d}"].append(face)
                        else:
                            inlet["coflow_inlet"].append(face)
                    if i == NX - 1:
                        outlet.append((cell[1], cell[2], cell[6], cell[5]))
                    if j == 0:
                        farfield.append((cell[0], cell[1], cell[5], cell[4]))
                    if j == NY - 1:
                        farfield.append((cell[3], cell[7], cell[6], cell[2]))
                    if k == 0:
                        front_back.append((cell[0], cell[3], cell[2], cell[1]))
                    if k == NZ - 1:
                        front_back.append((cell[4], cell[5], cell[6], cell[7]))

                    perturbation = 0.01 * math.sin(2.0 * math.pi * xc / 4.0)
                    perturbation *= math.exp(-((abs_y - 0.5) / 0.08) ** 2)
                    cell_id = len(cells) - 1
                    initial.write(f"{cell_id} {ux:.12e} {perturbation:.12e} 0\n")

    patches = [
        ("jet_core", "inlet", inlet["jet_core"]),
        *((f"inlet_ramp_{r:02d}", "inlet", inlet[f"inlet_ramp_{r:02d}"])
          for r in range(RAMP_BANDS)),
        ("coflow_inlet", "inlet", inlet["coflow_inlet"]),
        ("outlet", "outlet", outlet),
        ("farfield", "symmetry", farfield),
        ("frontBack", "symmetry", front_back),
    ]
    cell_count = NX * NY * NZ
    if cell_count != 100_000 or len(vertices) != (NX + 1) * (NY + 1) * (NZ + 1):
        raise RuntimeError("unexpected structured mesh dimensions")
    if not all(value > 0.0 for value in (DX, DY, DZ)):
        raise RuntimeError("all structured mesh spacings must be positive")
    if any(not faces for _, _, faces in patches):
        raise RuntimeError("every mesh patch must contain at least one face")

    mesh_dir = CASE / "mesh"
    mesh_dir.mkdir(parents=True, exist_ok=True)
    write_hex_v3(mesh_dir / "planar_jet.mesh", vertices, cells, patches)
    report = {
        "mesh_format": "BABELSIM_MESH 3",
        "cells": cell_count,
        "hexahedra": cell_count,
        "non_hexahedral_cells": 0,
        "vertices": len(vertices),
        "dimensions": [NX, NY, NZ],
        "bounds_D": {"x": [X0 / D, X1 / D], "y": [Y0 / D, Y1 / D], "z": [Z0 / D, Z1 / D]},
        "spacing_D": {"dx": DX / D, "dy": DY / D, "dz": DZ / D},
        "cell_volume_D3": DX * DY * DZ / D**3,
        "minimum_cell_volume_D3": DX * DY * DZ / D**3,
        "maximum_cell_volume_D3": DX * DY * DZ / D**3,
        "maximum_nonorthogonality_deg": 0.0,
        "maximum_skewness": 0.0,
        "orthogonal_hex_by_structured_construction": True,
        "patch_face_counts": {name: len(faces) for name, _, faces in patches},
        "jet_width_D": JET_WIDTH / D,
        "jet_cells_across": JET_WIDTH / DY,
        "discrete_inlet_mass_flow_rho_Uj_D2": inlet_mass_flow * D / 0.05,
        "initial_perturbation": {
            "type": "transverse velocity near both shear layers",
            "amplitude_Uj": 0.01,
            "wavelength_D": 4.0,
            "field": "Uy",
            "file": str(initial_path.relative_to(CASE)),
        },
        "note": "Orthogonality follows from axis-aligned uniform Cartesian coordinates; the runtime mesh reader is also exercised by the solver run.",
    }
    (mesh_dir / "mesh_quality.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
