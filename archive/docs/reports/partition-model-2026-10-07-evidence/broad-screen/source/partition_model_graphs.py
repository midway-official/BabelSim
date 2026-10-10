"""Dual graphs of realizable face-based polyhedral meshes for partition studies.

Every graph comes from a conforming Cartesian submesh, triangular/tetrahedral
cell subdivision, or extrusion of a bounded convex Voronoi tessellation. Native
BABELSIM_MESH fixtures allow the geometry and adjacency to be checked by Mesh.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
import json
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import Voronoi

FAMILIES = ('hex_cube', 'thin_plate', 'long_duct', 'rotated_duct', 'obstructed_cube',
            'branching_duct', 'graded_channel', 'warped_prisms', 'tet_cube',
            'voronoi_prisms', 'graded_voronoi', 'disconnected', 'renumbered_cube')


@dataclass
class Graph:
    offsets: np.ndarray
    neighbours: np.ndarray
    centres: np.ndarray
    face_counts: np.ndarray
    metadata: dict
    recipe: dict | None = None

    @property
    def cells(self):
        return len(self.face_counts)

    def csr(self):
        return csr_matrix((np.ones(len(self.neighbours), dtype=np.uint8),
                           self.neighbours, self.offsets), shape=(self.cells, self.cells))

    def validate(self):
        assert self.offsets.dtype == self.neighbours.dtype == np.dtype('int32')
        assert len(self.offsets) == self.cells + 1 and self.offsets[0] == 0
        assert self.offsets[-1] == len(self.neighbours)
        assert np.all(np.diff(self.offsets) <= self.face_counts)
        assert np.all(self.face_counts >= 4) and np.isfinite(self.centres).all()
        assert not len(self.neighbours) or (self.neighbours.min() >= 0 and self.neighbours.max() < self.cells)
        matrix = self.csr()
        assert matrix.diagonal().sum() == 0 and (matrix != matrix.T).nnz == 0
        self.metadata.update(cells=self.cells, internalFaces=len(self.neighbours)//2,
            physicalFaces=int((self.face_counts-np.diff(self.offsets)).sum()),
            components=int(connected_components(matrix, directed=False, return_labels=False)),
            minCellFaces=int(self.face_counts.min()), maxCellFaces=int(self.face_counts.max()))

    def save(self, path):
        self.validate()
        np.savez(path, offsets=self.offsets, neighbours=self.neighbours,
            centres=self.centres, face_counts=self.face_counts,
            metadata=json.dumps(self.metadata, sort_keys=True))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as data:
            return cls(data['offsets'], data['neighbours'], data['centres'], data['face_counts'],
                       json.loads(str(data['metadata'])))


def graph_from_edges(cells, pairs, centres, face_counts, metadata, recipe=None):
    u, v = pairs[:, 0], pairs[:, 1]
    rows, cols = np.concatenate((u, v)), np.concatenate((v, u))
    matrix = csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cols)), shape=(cells, cells))
    matrix.sort_indices()
    assert np.all(matrix.data == 1), 'duplicate internal cell face'
    return Graph(matrix.indptr.astype(np.int32), matrix.indices.astype(np.int32),
                 np.asarray(centres, dtype=np.float64), np.asarray(face_counts, dtype=np.int32), metadata, recipe)


def rotation():
    a, b = 0.67, 0.38
    rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    return rz @ ry


def voxel_graph(shape, name, scale, *, mask=None, graded=False, rotate=False):
    nx, ny, nz = shape
    k, j, i = np.indices((nz, ny, nx))
    keep = np.ones((nz, ny, nx), bool) if mask is None else mask(i, j, k, nx, ny, nz)
    mapping = np.full((nz, ny, nx), -1, dtype=np.int32)
    mapping[keep] = np.arange(keep.sum(), dtype=np.int32)
    edges = []
    for axis in range(3):
        left, right = [slice(None)]*3, [slice(None)]*3
        left[axis], right[axis] = slice(None, -1), slice(1, None)
        u, v = mapping[tuple(left)].ravel(), mapping[tuple(right)].ravel()
        ok = (u >= 0) & (v >= 0)
        edges.append(np.column_stack((u[ok], v[ok])))
    physical = (8., 1., 1.) if name in ('long_duct', 'rotated_duct', 'graded_channel') else (1., 1., .01 if nz == 1 else 1.)
    axes = [np.linspace(0, length, n+1) for n, length in zip(shape, physical)]
    if graded:
        axes[0] = physical[0] * np.linspace(0, 1, nx+1)**3
    points = np.column_stack(tuple(((axis[:-1]+axis[1:])/2)[index[keep]] for axis, index in zip(axes, (i, j, k))))
    if rotate:
        points = points @ rotation().T
    recipe = dict(kind='voxel', axes=axes, kept=np.column_stack((i[keep], j[keep], k[keep])), rotate=rotate)
    return graph_from_edges(len(points), np.concatenate(edges), points, np.full(len(points), 6),
        dict(family=name, scale=scale, construction='conforming Cartesian submesh', shape=shape,
             graded=graded, rotated=rotate), recipe)


def prism_graph(n, scale):
    from generate_mpi_scaling_case import grid_triangles, warped_point
    triangles = list(grid_triangles(n))
    base = np.array([warped_point(i,j,n,.15) for j in range(n+1) for i in range(n+1)])
    edge_map, pairs = {}, []
    for owner, triangle in enumerate(triangles):
        for edge in range(3):
            key = tuple(sorted((triangle[edge], triangle[(edge+1)%3])))
            if key in edge_map: pairs.append((edge_map.pop(key), owner))
            else: edge_map[key] = owner
    base_pairs = np.array(pairs, dtype=np.int32)
    count = len(triangles)
    edges = [base_pairs+k*count for k in range(n)]
    ids = np.arange(count*n, dtype=np.int32).reshape(n,count)
    edges.append(np.column_stack((ids[:-1].ravel(),ids[1:].ravel())))
    centres = base[np.array(triangles)].mean(axis=1)
    points = np.column_stack((np.tile(centres,(n,1)),np.repeat((np.arange(n)+.5)/n,count)))
    return graph_from_edges(count*n, np.concatenate(edges), points, np.full(count*n,5),
        dict(family='warped_prisms',scale=scale,construction='warped checkerboard triangular-prism extrusion',n=n,warp=.15),
        dict(kind='prism',n=n))


TETS = tuple((0, 1 << p[0], (1 << p[0]) | (1 << p[1]), 7) for p in permutations(range(3)))
TET_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3))


def tet_graph(n, scale):
    ids = np.arange(6*n**3, dtype=np.int32).reshape(n, n, n, 6)
    lookup = {}
    for t, corners in enumerate(TETS):
        for face in TET_FACES:
            lookup.setdefault(tuple(sorted(corners[a] for a in face)), []).append(t)
    edges = []
    for ring, cells in lookup.items():
        if len(cells) == 2:
            edges.append(np.column_stack((ids[..., cells[0]].ravel(), ids[..., cells[1]].ravel())))
        else:
            for axis in range(3):
                if all(v & (1 << axis) for v in ring):
                    opposite = tuple(sorted(v ^ (1 << axis) for v in ring))
                    neighbour = lookup[opposite][0]
                    lower, upper = [slice(None)]*3, [slice(None)]*3
                    lower[2-axis], upper[2-axis] = slice(None, -1), slice(1, None)
                    edges.append(np.column_stack((ids[tuple(lower)+(cells[0],)].ravel(), ids[tuple(upper)+(neighbour,)].ravel())))
                    break
    k, j, i = np.indices((n, n, n))
    points = np.empty((n, n, n, 6, 3))
    for t, corners in enumerate(TETS):
        shift = np.mean([[v & 1, (v >> 1) & 1, (v >> 2) & 1] for v in corners], axis=0)
        points[..., t, :] = np.stack(((i+shift[0])/n, (j+shift[1])/n, (k+shift[2])/n), axis=-1)
    cells = 6*n**3
    return graph_from_edges(cells, np.concatenate(edges), points.reshape(cells, 3), np.full(cells, 4),
        dict(family='tet_cube', scale=scale, construction='conforming six-tetrahedra subdivision of each cube', n=n),
        dict(kind='tet', n=n))


def voronoi_graph(count, layers, scale, graded=False, seed=20261007):
    rng = np.random.default_rng(seed)
    sites = rng.random((count, 2))
    if graded:
        sites[:, 0] = sites[:, 0]**3
    mirrors = [sites, sites*np.array([-1, 1]), sites*np.array([-1, 1])+[2, 0],
               sites*np.array([1, -1]), sites*np.array([1, -1])+[0, 2]]
    vor = Voronoi(np.concatenate(mirrors))
    regions = [vor.regions[vor.point_region[i]] for i in range(count)]
    assert all(len(r) >= 3 and -1 not in r for r in regions)
    base_pairs = vor.ridge_points[(vor.ridge_points[:, 0] < count) & (vor.ridge_points[:, 1] < count)]
    base_pairs = np.sort(base_pairs, axis=1).astype(np.int32)
    centres = []
    for region in regions:
        p = vor.vertices[region]
        q = np.roll(p, -1, axis=0)
        cross = p[:, 0]*q[:, 1]-p[:, 1]*q[:, 0]
        centres.append(((p+q)*cross[:, None]).sum(axis=0)/(3*cross.sum()))
    edges = [base_pairs + k*count for k in range(layers)]
    ids = np.arange(count*layers, dtype=np.int32).reshape(layers, count)
    if layers > 1:
        edges.append(np.column_stack((ids[:-1].ravel(), ids[1:].ravel())))
    points = np.column_stack((np.tile(np.asarray(centres), (layers, 1)), np.repeat((np.arange(layers)+.5)/layers, count)))
    face_counts = np.tile(np.array([len(r)+2 for r in regions], dtype=np.int32), layers)
    name = 'graded_voronoi' if graded else 'voronoi_prisms'
    return graph_from_edges(count*layers, np.concatenate(edges), points, face_counts,
        dict(family=name, scale=scale, construction='bounded convex Voronoi-polygon extrusion',
             sites=count, layers=layers, graded=graded, seed=seed),
        dict(kind='voronoi', vertices=vor.vertices, regions=regions, layers=layers))


def make_graph(family, scale):
    sizes = {'fixture': (8, 4, 32, 8, 2, 32, 2), 'small': (24, 12, 128, 64, 8, 256, 8),
             'engineering': (80, 44, 1024, 384, 16, 16000, 32)}
    cube, tet, duct, plate, width, sites, layers = sizes[scale]
    if family == 'warped_prisms':
        return prism_graph({'fixture': 4, 'small': 16, 'engineering': 63}[scale], scale)
    if family == 'tet_cube':
        return tet_graph(tet, scale)
    if family in ('voronoi_prisms', 'graded_voronoi'):
        return voronoi_graph(sites, layers, scale, graded=family == 'graded_voronoi')
    if family == 'thin_plate':
        return voxel_graph((plate, plate, 1), family, scale)
    if family in ('long_duct', 'rotated_duct', 'graded_channel'):
        return voxel_graph((duct, width, width), family, scale, graded=family == 'graded_channel', rotate=family == 'rotated_duct')
    if family == 'obstructed_cube':
        return voxel_graph((cube, cube, cube), family, scale,
            mask=lambda i,j,k,nx,ny,nz: ((i+.5)/nx-.5)**2+((j+.5)/ny-.5)**2 > .18**2)
    if family == 'branching_duct':
        return voxel_graph((plate, plate, max(2, width//2)), family, scale,
            mask=lambda i,j,k,nx,ny,nz: (abs((j+.5)/ny-.65)<.10) | ((abs((i+.5)/nx-.5)<.10) & ((j+.5)/ny<=.65)))
    if family == 'disconnected':
        first = voxel_graph((cube, cube, max(1,cube//2)), family, scale)
        second = voxel_graph((max(2,cube//2),)*3, family, scale)
        matrix = csr_matrix((np.ones(len(first.neighbours), dtype=np.uint8), first.neighbours, first.offsets))
        from scipy.sparse import block_diag
        both = block_diag((matrix, second.csr()), format='csr')
        return Graph(both.indptr.astype(np.int32), both.indices.astype(np.int32),
            np.vstack((first.centres, second.centres+[2,0,0])), np.concatenate((first.face_counts, second.face_counts)),
            dict(family=family, scale=scale, construction='two disconnected conforming hex domains'),
            dict(kind='disconnected', components=[first.recipe, second.recipe]))
    graph = voxel_graph((cube,)*3, family, scale)
    if family == 'renumbered_cube':
        perm = np.random.default_rng(73021).permutation(graph.cells)
        matrix = graph.csr()[perm][:, perm]
        matrix.sort_indices()
        graph = Graph(matrix.indptr.astype(np.int32), matrix.indices.astype(np.int32), graph.centres[perm],
            graph.face_counts[perm], dict(family=family, scale=scale, construction='hex cube with scrambled cell IDs', seed=73021),
            dict(kind='renumbered', source=graph.recipe, perm=perm))
    return graph


def native_fixture(graph, path):
    """Produce actual face rings; the graph is checked independently from them."""
    from polyhedral_mesh import HEX_FACES, write_v3
    vertices, cells = [], []
    vertex_map = {}
    def vertex(point):
        key = tuple(float(x) for x in point)
        if key not in vertex_map:
            vertex_map[key] = len(vertices); vertices.append(key)
        return vertex_map[key]
    def voxel(recipe, translation=(0,0,0)):
        axes = recipe['axes']
        for i,j,k in recipe['kept']:
            rings = []
            corners = []
            for dx,dy,dz in ((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)):
                point = np.array([axes[0][i+dx], axes[1][j+dy], axes[2][k+dz]])
                if recipe['rotate']: point = rotation() @ point
                corners.append(vertex(point+translation))
            rings = [tuple(corners[a] for a in face) for face in HEX_FACES]
            cells.append(rings)
    recipe = graph.recipe
    if recipe['kind'] in ('voxel','disconnected','renumbered'):
        if recipe['kind'] == 'disconnected':
            voxel(recipe['components'][0]); voxel(recipe['components'][1], (2,0,0))
        elif recipe['kind'] == 'renumbered':
            voxel(recipe['source']); cells = [cells[int(i)] for i in recipe['perm']]
        else: voxel(recipe)
    elif recipe['kind'] == 'tet':
        n = recipe['n']
        for k in range(n):
            for j in range(n):
                for i in range(n):
                    for corners in TETS:
                        ids = [vertex(((i+(v&1))/n,(j+((v>>1)&1))/n,(k+((v>>2)&1))/n)) for v in corners]
                        cells.append([tuple(ids[a] for a in f) for f in TET_FACES])
    else:
        if recipe['kind'] == 'prism':
            n = recipe['n']; regions = []; base = []
            from generate_mpi_scaling_case import grid_triangles, warped_point
            for j in range(n+1):
                for i in range(n+1): base.append(warped_point(i,j,n,.15))
            regions = list(grid_triangles(n)); layers = n
        else:
            base, regions, layers = recipe['vertices'], recipe['regions'], recipe['layers']
        for k in range(layers):
            for region in regions:
                lower = [vertex((*base[v], k/layers)) for v in region]
                upper = [vertex((*base[v], (k+1)/layers)) for v in region]
                rings = [tuple(reversed(lower)), tuple(upper)]
                for a in range(len(region)):
                    b=(a+1)%len(region); rings.append((lower[a],lower[b],upper[b],upper[a]))
                cells.append(rings)
    records = {}
    for owner, rings in enumerate(cells):
        for ring in rings:
            key=tuple(sorted(ring))
            if key not in records: records[key]=[ring, owner, -1, 0]
            else:
                assert records[key][2] == -1, 'nonmanifold generated face'
                records[key][2], records[key][3] = owner, -1
    actual = np.array([(v[1],v[2]) for v in records.values() if v[2]>=0], dtype=np.int32).reshape(-1,2)
    checked = graph_from_edges(len(cells), actual, graph.centres, np.array([len(c) for c in cells]), {})
    assert np.array_equal(graph.offsets,checked.offsets) and np.array_equal(graph.neighbours,checked.neighbours)
    assert np.array_equal(graph.face_counts,checked.face_counts)
    write_v3(path, vertices, list(records.values()), [('boundary','generic')])
