#!/usr/bin/env python3
"""从内置方腔案例生成指定 Re、网格和迎风格式的可复现实验案例。"""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path


def coordinates(cells: int, clustering: float) -> list[float]:
    if clustering == 0.0:
        return [index / cells for index in range(cells + 1)]
    scale = math.tanh(clustering)
    return [
        0.5 * (1.0 + math.tanh(
            clustering * (2.0 * index / cells - 1.0)) / scale)
        for index in range(cells + 1)
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cells", type=int, required=True)
    parser.add_argument("--re", type=float, required=True)
    parser.add_argument("--cluster", type=float, default=1.5,
                        help="双曲正切壁面加密强度，0 表示均匀网格")
    parser.add_argument("--convection", choices=("upwind", "linearUpwind"),
                        default="linearUpwind")
    parser.add_argument("--gradient", choices=("greenGauss", "leastSquares"),
                        default="greenGauss")
    parser.add_argument("--velocity-relaxation", type=float, default=0.3)
    parser.add_argument("--pressure-relaxation", type=float, default=0.3)
    parser.add_argument("--max-iterations", type=int, default=30000)
    args = parser.parse_args()
    if args.cells < 4 or args.re <= 0.0 or args.cluster < 0.0:
        parser.error("cells 必须至少为 4，Re 必须为正数，cluster 不能为负")
    if not 0.0 < args.velocity_relaxation <= 1.0 or \
       not 0.0 < args.pressure_relaxation <= 1.0:
        parser.error("松弛因子必须位于 (0, 1]")
    if args.max_iterations <= 0:
        parser.error("max-iterations 必须为正数")

    repository = Path(__file__).resolve().parents[3]
    if args.output.exists():
        parser.error(f"输出目录已存在：{args.output}")
    shutil.copytree(
        repository / "cases/cavity", args.output,
        ignore=shutil.ignore_patterns("results", "post", "validation"))

    (args.output / "physics/simple.bs").write_text(
        f"density 1.0\ndynamicViscosity {1.0 / args.re:.17g}\n",
        encoding="utf-8")
    (args.output / "numerics/methods.bs").write_text(
        "interpolation linear\n"
        f"gradient {args.gradient}\n"
        f"convection {args.convection}\n"
        "diffusion orthogonal\ntime steady\n",
        encoding="utf-8")
    (args.output / "numerics/solution.bs").write_text(
        f"maxIterations {args.max_iterations}\n"
        f"velocityRelaxation {args.velocity_relaxation:.17g}\n"
        f"pressureRelaxation {args.pressure_relaxation:.17g}\n"
        "continuityTolerance 1e-10\nvelocityTolerance 2e-7\n"
        "vectorSolver bicgstab ilut 1e-13 1e-9 2000\n"
        "scalarSolver cg incompleteCholesky 1e-13 1e-9 2000\n",
        encoding="utf-8")

    xy = coordinates(args.cells, args.cluster)
    mesh = args.output / "mesh/cavity.mesh"
    with mesh.open("w", encoding="utf-8") as output:
        count = args.cells
        def vertex(i: int, j: int, k: int) -> int:
            return i + (count + 1) * (j + (count + 1) * k)
        points = [(x, y, z) for z in (0.0, 1.0) for y in xy for x in xy]
        cells = []
        boundaries = [[] for _ in range(6)]
        for k in range(1):
            for j in range(count):
                for i in range(count):
                    cell = (vertex(i, j, k), vertex(i + 1, j, k),
                            vertex(i + 1, j + 1, k), vertex(i, j + 1, k),
                            vertex(i, j, k + 1), vertex(i + 1, j, k + 1),
                            vertex(i + 1, j + 1, k + 1), vertex(i, j + 1, k + 1))
                    cells.append(cell)
                    if i == 0: boundaries[0].append((cell[0], cell[4], cell[7], cell[3]))
                    if i + 1 == count: boundaries[1].append((cell[1], cell[2], cell[6], cell[5]))
                    if j == 0: boundaries[2].append((cell[0], cell[1], cell[5], cell[4]))
                    if j + 1 == count: boundaries[3].append((cell[3], cell[7], cell[6], cell[2]))
                    boundaries[4].append((cell[0], cell[3], cell[2], cell[1]))
                    boundaries[5].append((cell[4], cell[5], cell[6], cell[7]))
        output.write(f"BABELSIM_MESH 2\nvertices {len(points)}\n")
        for point in points:
            output.write(" ".join(f"{value:.17g}" for value in point) + "\n")
        output.write(f"cells {len(cells)}\n")
        for cell in cells:
            output.write(" ".join(map(str, cell)) + "\n")
        output.write("patches 6\n")
        for name, kind, faces in zip(
                ("cavity_left", "cavity_right", "cavity_bottom", "lid", "front", "back"),
                ("wall", "wall", "wall", "wall", "symmetry", "symmetry"), boundaries):
            output.write(f"patch {name} {kind} {len(faces)}\n")
            for face in faces:
                output.write(" ".join(map(str, face)) + "\n")
        output.write("end\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
