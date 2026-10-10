#!/usr/bin/env python3
"""Create continuous near-wall k/omega seeds from the existing extruded mesh.

This does not change the mesh or generate a developed boundary-layer solution.
U and p are initialized by the PISO potential projection at runtime.
"""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
mesh_path = ROOT / 'mesh/naca0012.mesh'
with mesh_path.open() as stream:
    if next(stream).strip() != 'BABELSIM_MESH 3':
        raise ValueError('expected the existing explicit-face mesh')
    assert next(stream).strip() == 'vertices'
    vertices = np.array([[float(x) for x in next(stream).split()]
                         for _ in range(int(next(stream)))])
    assert next(stream).strip() == 'faces'
    centres = {}
    wall_edges = []
    for _ in range(int(next(stream))):
        row = next(stream).split()
        owner, patch = int(row[2]), int(row[4])
        points = vertices[[int(x) for x in row[5:]], :2]
        if patch == 0:  # airfoil, checked below
            wall_edges.append(np.unique(points, axis=0))
        if patch == 4:  # front plane: exact extruded-cell planar centroid
            following = np.roll(points, -1, axis=0)
            cross = points[:, 0]*following[:, 1] - following[:, 0]*points[:, 1]
            centres[owner] = ((points+following)*cross[:, None]).sum(axis=0)/(3*cross.sum())
    assert next(stream).strip() == 'patches'
    patches = [next(stream).split()[1] for _ in range(int(next(stream)))]
    assert patches == ['airfoil', 'inlet', 'outlet', 'farfield', 'front', 'back']
ids = np.array(sorted(centres))
assert np.array_equal(ids, np.arange(len(ids)))
xy = np.array([centres[i] for i in ids])
distance = np.full(len(ids), np.inf)
for edge in wall_edges:
    assert edge.shape == (2, 2)
    start, end = edge
    vector = end-start
    fraction = np.clip(((xy-start)*vector).sum(axis=1)/(vector@vector), 0, 1)
    distance = np.minimum(distance, np.linalg.norm(xy-start-fraction[:, None]*vector, axis=1))

# k vanishes continuously at the wall. A shifted viscous omega asymptote
# matches the finite wall value at d=0 and approaches the inlet value away
# from the wall. These are numerical initial seeds, not measured turbulence.
k_min, k_inlet, omega_inlet = 1e-12, 1.5e-6, 1.5
wall_offset, k_transition = 1.5e-5, 0.001
k = k_min+(k_inlet-k_min)*(-np.expm1(-distance/k_transition))**2
omega = np.maximum(omega_inlet, (6e-6/0.075)/(distance+wall_offset)**2)
folder = ROOT/'fields/initial'
for name, values in [('k', k), ('omega', omega)]:
    np.savetxt(folder/(name+'.dat'), np.column_stack((ids, values)), fmt=['%d', '%.16g'])
report = dict(mesh_sha256=hashlib.sha256(mesh_path.read_bytes()).hexdigest(),
              cells=len(ids), minimum_wall_distance=float(distance.min()),
              k_transition=k_transition, omega_wall_offset=wall_offset,
              k_range=[float(k.min()), float(k.max())],
              omega_range=[float(omega.min()), float(omega.max())])
(folder/'initialization.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
