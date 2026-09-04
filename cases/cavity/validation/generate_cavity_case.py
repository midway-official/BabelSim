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
        output.write(
            "BABELSIM_MESH 1\n"
            f"dimensions {args.cells} {args.cells} 1\n")
        if args.cluster == 0.0:
            output.write("geometry cartesian\nbounds 0 0 0 1 1 1\n")
        else:
            output.write("geometry vertices\n")
            for z in (0.0, 1.0):
                for y in xy:
                    for x in xy:
                        output.write(f"{x:.17g} {y:.17g} {z:.17g}\n")
        output.write(
            "patch xmin cavity_left wall\n"
            "patch xmax cavity_right wall\n"
            "patch ymin cavity_bottom wall\n"
            "patch ymax lid wall\n"
            "patch zmin front symmetry\n"
            "patch zmax back symmetry\nend\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
