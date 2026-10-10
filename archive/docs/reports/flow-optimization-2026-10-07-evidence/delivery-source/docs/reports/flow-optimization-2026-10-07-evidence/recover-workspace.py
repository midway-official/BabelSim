#!/usr/bin/env python3
"""Restore checksummed engineering inputs after the temporary workspace was lost.

Run before any solver benchmark. Persistent copies remain outside the source
checkout; temporary paths use hard links so the original study can resume.
"""
from pathlib import Path
import hashlib
import json
import shutil
import sys

EVIDENCE = Path(__file__).resolve().parent
REPOSITORY = EVIDENCE.parents[2]
sys.path.insert(0, str(REPOSITORY/'tools'))
from benchmark_backend import case_manifest, sha256
from generate_mpi_scaling_case import write_mesh
from partition_model_graphs import make_graph, native_fixture
from generate_partition_flow_case import BOX_PATCHES

PERSISTENT = Path('/home/midway/babelsim-flow-opt-20261007')
TEMPORARY = Path('/tmp/babelsim-flow-opt-20261007')


def restore_case(study_name, label, mesh, graph=None):
    study = EVIDENCE/study_name
    specification = json.loads((study/'study.json').read_text())['cases'][label]
    relative = Path(specification['path']).relative_to(TEMPORARY)
    persistent = PERSISTENT/relative
    temporary = TEMPORARY/relative
    for entry in specification['manifest']['files']:
        rel = Path(entry['path'])
        source = mesh if rel.suffix == '.mesh' else graph if rel.suffix == '.npz' else study/'inputs'/label/rel
        assert source and source.is_file(), (rel, source)
        assert source.stat().st_size == entry['bytes'] and sha256(source) == entry['sha256'], rel
        target = persistent/rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            target.hardlink_to(source) if rel.suffix in ('.mesh', '.npz') else shutil.copy2(source, target)
        link = temporary/rel
        link.parent.mkdir(parents=True, exist_ok=True)
        if not link.exists():
            link.hardlink_to(target)
    assert case_manifest(temporary) == specification['manifest']
    print('RESTORED', label, str(persistent), specification['manifest']['sha256'], flush=True)


def main():
    meshes = PERSISTENT/'recovered-meshes'
    meshes.mkdir(parents=True, exist_ok=True)
    prisms = {}
    for n in (63, 80):
        path = meshes/f'prism{n}.mesh'
        if not path.exists():
            print('GENERATING PRISM', n, flush=True)
            write_mesh(path, n, .15, 'piso')
        prisms[n] = path
    voronoi, graph_path = meshes/'voronoi.mesh', meshes/'dual-graph.npz'
    if not voronoi.exists() or not graph_path.exists():
        print('GENERATING VORONOI', flush=True)
        graph = make_graph('voronoi_prisms', 'engineering')
        graph.validate()
        if not voronoi.exists():
            native_fixture(graph, voronoi, box_boundary_patches=BOX_PATCHES)
        if not graph_path.exists():
            graph.save(graph_path)
    restore_case('final-strong-scaling', 'simple', prisms[63])
    restore_case('final-strong-scaling', 'piso', prisms[63])
    restore_case('amg-strong-scaling', 'prism', prisms[80])
    restore_case('amg-strong-scaling', 'voronoi', voronoi, graph_path)
    binary_name = 'babelsim-memory-partition-normalized'
    reference = PERSISTENT/'reference-source/build-petsc/babelsim-solve'
    expected = '7a29412a574f912f6d67c7b3514b7aa15c60b854c222a455392ed977cb61a6b2'
    assert sha256(reference) == expected
    persistent = PERSISTENT/'bin'/binary_name
    temporary = TEMPORARY/'bin'/binary_name
    if not persistent.exists():
        shutil.copy2(reference, persistent)
    if not temporary.exists():
        temporary.hardlink_to(persistent)
    info = json.loads((EVIDENCE/'recovery-2026-10-08.json').read_text())
    info.update(referenceBinarySha256=expected, referenceRecovery='byte-identical rebuild confirmed',
                restoredInputManifests=True,
                meshes={p.name:sha256(p) for p in (*prisms.values(), voronoi, graph_path)})
    (EVIDENCE/'recovery-2026-10-08.json').write_text(json.dumps(info, indent=2)+'\n')
    print('ALL INPUT AND BINARY CHECKSUMS MATCH', flush=True)


if __name__ == '__main__':
    main()
