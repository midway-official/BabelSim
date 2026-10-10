#!/usr/bin/env python3
"""Prepare identical, fixed-work production cases for every scaling variant."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'tools'))
from benchmark_backend import case_manifest, sha256


def clone(template, target, mesh_source=None):
    assert not target.exists(), target
    manifest = case_manifest(template)
    target.mkdir(parents=True)
    for item in manifest['files']:
        relative = Path(item['path'])
        destination = target / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if relative.suffix == '.mesh' and mesh_source is not None:
            rename = {'patch hot wall': 'patch inlet inlet',
                      'patch cold wall': 'patch outlet outlet',
                      'patch lower symmetry': 'patch lower_wall wall',
                      'patch upper symmetry': 'patch upper_wall wall'}
            with mesh_source.open() as source, destination.open('w') as out:
                for line in source:
                    replacement = rename.get(line.strip())
                    out.write(replacement+'\n' if replacement else line)
        elif relative.suffix in ('.mesh', '.npz'):
            destination.hardlink_to(template / relative)
        else:
            shutil.copy2(template / relative, destination)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=10)
    args = parser.parse_args()
    assert args.steps >= 2
    root = args.root.resolve()
    specs = [('prism', root/'formal-inputs/piso6', root/'mesh-inputs/prism80/mesh.mesh',
              dict(cells=1024000, faces=2585600, gridN=80, minCellFaces=5, maxCellFaces=5)),
             ('voronoi', root/'formal-inputs/voronoi-piso6', None,
              dict(cells=512000, faces=2064032, minCellFaces=5, maxCellFaces=15))]
    result = {}
    for name, template, mesh, geometry in specs:
        target = root/f'strong-inputs/{name}-piso{args.steps}'
        manifest = clone(template, target, mesh)
        (target/'control.bs').write_text(f'startTime 0\nendTime {args.steps*.01:.12g}\ndeltaT 0.01\n')
        text = (target/'output.bs').read_text()
        text = '\n'.join(f'writeInterval {args.steps}' if line.startswith('writeInterval ') else line
                         for line in text.splitlines())+'\n'
        (target/'output.bs').write_text(text)
        info = dict(geometry, template=str(template), templateInputSha256=manifest['sha256'],
                    acceptedSteps=args.steps, measuredIntervals=args.steps-1, deltaT=.01,
                    unchangedSolutionControls=True, meshSha256=sha256(target/'mesh/flow.mesh'))
        (target/'generation.json').write_text(json.dumps(info, indent=2)+'\n')
        result[name] = dict(case=str(target), **info)
        print(json.dumps(result[name]), flush=True)
    Path(__file__).with_name(f'strong-case-inputs-{args.steps}.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
