#!/usr/bin/env python3
"""Create a production flow case on an independently checked partition graph.

The solver dictionaries come from an existing case. Only its mesh is replaced;
all equation, coupling and convergence controls are retained.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

from benchmark_backend import case_manifest
from partition_model_graphs import make_graph, native_fixture


BOX_PATCHES = [('inlet', 'inlet'), ('outlet', 'outlet'),
               ('lower_wall', 'wall'), ('upper_wall', 'wall'),
               ('front', 'symmetry'), ('back', 'symmetry')]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--template', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--family', choices=('voronoi_prisms', 'graded_voronoi',
                                            'tet_cube', 'hex_cube'), required=True)
    parser.add_argument('--scale', choices=('fixture', 'small', 'engineering'),
                        default='engineering')
    args = parser.parse_args()
    template, output = args.template.resolve(), args.output.resolve()
    if output.exists():
        parser.error('output already exists; use a new directory')
    manifest = case_manifest(template)
    meshes = [entry for entry in manifest['files'] if entry['path'].endswith('.mesh')]
    if len(meshes) != 1:
        parser.error('template must contain exactly one native .mesh input')
    output.mkdir(parents=True)
    for entry in manifest['files']:
        if entry == meshes[0] or entry['path'] == 'generation.json':
            continue
        relative = Path(entry['path'])
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(template / relative, target)
    graph = make_graph(args.family, args.scale)
    graph.validate()
    mesh = output / meshes[0]['path']
    mesh.parent.mkdir(parents=True, exist_ok=True)
    native_fixture(graph, mesh, box_boundary_patches=BOX_PATCHES)
    graph.save(output / 'dual-graph.npz')
    info = dict(graph.metadata, template=str(template),
                templateInputSha256=manifest['sha256'], mesh=str(mesh),
                meshSha256=hashlib.sha256(mesh.read_bytes()).hexdigest(),
                csrAndFaceCountsIndependentlyChecked=True,
                boundaryPlanes=[name for name, role in BOX_PATCHES])
    (output / 'generation.json').write_text(json.dumps(info, indent=2) + '\n')
    print(json.dumps(info, indent=2), flush=True)


if __name__ == '__main__':
    main()
