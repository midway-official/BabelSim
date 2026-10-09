#!/usr/bin/env python3
"""Recreate the accepted 2026-09-28 solver-only experiments in a fresh directory."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import generate_mpi_scaling_case as heat


def variant(base, target, old, new, extra=''):
    shutil.copytree(base, target)
    solution = target / 'numerics/solution.bs'
    text = solution.read_text()
    if text.count(old) != 1:
        raise ValueError(f'expected exactly one configuration entry: {old}')
    solution.write_text(text.replace(old, new) + extra)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--include-naca', action='store_true', help='also copy current NACA inputs for a five-step pressure-PC comparison')
    args = parser.parse_args()
    root = args.output.resolve()
    if root.exists():
        parser.error('output must be a new directory')
    root.mkdir(parents=True)
    for label, n in [('101', 37), ('500', 63)]:
        base = root / f'heat{label}-hypre'
        base.mkdir()
        mesh = base / 'mesh.mesh'
        info = heat.write_mesh(mesh, n, 0.15, 'heat')
        heat.write_case(base, mesh)
        info.update(case=str(base), meshBytes=mesh.stat().st_size)
        (base / 'generation.json').write_text(json.dumps(info, indent=2) + '\n')
        variant(base, root / f'heat{label}-bjacobi', 'pcType hypre', 'pcType bjacobi')
        if label == '101':
            variant(base, root / 'heat101-gamg1', 'pcType hypre', 'pcType gamg',
                    'equation.temperature.amgSmoothingSteps 1\n')
    base = root / 'piso101-hypre'
    subprocess.run([sys.executable, str(Path(__file__).with_name('generate_mpi_flow_case.py')),
                    '--mesh-case', str(root / 'heat101-hypre'), '--output', str(base),
                    '--solver', 'piso'], check=True)
    (base / 'control.bs').write_text('startTime 0\nendTime 0.05\ndeltaT 0.01\n')
    solution = base / 'numerics/solution.bs'
    solution.write_text(solution.read_text().replace('nCorrectors 2', 'nCorrectors 8')
                        + 'couplingTolerance 1e-5\nmaxCorrectors 50\n')
    output = base / 'output.bs'
    output.write_text(output.read_text().replace('writeInterval 1', 'writeInterval 5'))
    variant(base, root / 'piso101-gamg1',
            'equation.pressureCorrection.pcType hypre',
            'equation.pressureCorrection.pcType gamg',
            'equation.pressureCorrection.amgSmoothingSteps 1\n')
    # Preserve the exact dictionary order of the measured exploratory case.
    candidate = root / 'piso101-gamg1/numerics/solution.bs'
    coupling = 'couplingTolerance 1e-5\nmaxCorrectors 50\n'
    candidate.write_text(candidate.read_text().replace(coupling, '') + coupling)
    if args.include_naca:
        source = Path(__file__).resolve().parents[1] / 'cases/naca0012'
        base = root / 'naca-hypre'
        base.mkdir()
        for name in ('case.bs', 'control.bs', 'output.bs'):
            shutil.copy2(source / name, base / name)
        for name in ('mesh', 'fields/initial', 'physics', 'numerics'):
            shutil.copytree(source / name, base / name)
        (base / 'control.bs').write_text('startTime 0\nendTime 0.05\ndeltaT 0.01\n')
        (base / 'output.bs').write_text('directory results\ntimeName final\nwriteInterval 5\n')
        variant(base, root / 'naca-gamg1',
                'equation.pressureCorrection.pcType hypre',
                'equation.pressureCorrection.pcType gamg',
                'equation.pressureCorrection.amgSmoothingSteps 1\n')
    print(f'Created {9 if args.include_naca else 7} cases in {root}')


if __name__ == '__main__':
    main()
