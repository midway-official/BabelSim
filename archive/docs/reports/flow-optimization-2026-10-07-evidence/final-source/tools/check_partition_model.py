#!/usr/bin/env python3
"""Independent checks for partition witnesses, exact quotas and stencil halos."""
from collections import deque
from pathlib import Path
import argparse
import os

for setting in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[setting] = '1'

import numpy as np
from partition_model_graphs import FAMILIES, Graph
from partition_model_algorithms import (graph_growth, halo_masks, quotas, rebalance,
                                        removable, warm_kernels)


def reference_halo(graph, owners, part, layers):
    distance = {int(cell): 0 for cell in np.flatnonzero(owners == part)}
    frontier = deque(distance)
    while frontier:
        cell = frontier.popleft()
        if distance[cell] == layers:
            continue
        for neighbour in graph.neighbours[graph.offsets[cell]:graph.offsets[cell + 1]]:
            neighbour = int(neighbour)
            if neighbour not in distance:
                distance[neighbour] = distance[cell] + 1
                frontier.append(neighbour)
    return set(distance)


def reference_removable(graph, owners, cell):
    source = owners[cell]
    goals = {int(other) for other in graph.neighbours[graph.offsets[cell]:graph.offsets[cell + 1]]
             if owners[other] == source}
    if len(goals) <= 1:
        return True
    seed = min(goals)
    visited = {cell, seed}
    frontier = deque([seed])
    while frontier:
        current = frontier.popleft()
        for other in graph.neighbours[graph.offsets[current]:graph.offsets[current + 1]]:
            other = int(other)
            if owners[other] == source and other not in visited:
                visited.add(other)
                frontier.append(other)
    return goals <= visited


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph-root', type=Path, required=True)
    args = parser.parse_args()
    warm_kernels()
    # A surplus and deficit separated by a balanced intermediate partition.
    offsets = np.array([0, 1, 3, 5, 7, 9, 11, 13, 15, 16], np.int32)
    neighbours = np.array([1, 0, 2, 1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7], np.int32)
    owners, transfers, unsafe = rebalance(offsets, neighbours,
        np.array([0, 0, 0, 0, 1, 1, 1, 2, 2], np.int32), quotas(9, 3))
    assert np.array_equal(np.bincount(owners, minlength=3), [3, 3, 3])
    assert transfers == 2 and unsafe == 0
    # Valid disconnected cells have no graph edges. Check this distinct path.
    owners, transfers, unsafe = rebalance(np.zeros(11, np.int32), np.empty(0, np.int32),
        np.array([0] * 8 + [1, 2], np.int32), quotas(10, 3))
    assert np.array_equal(np.bincount(owners, minlength=3), [4, 3, 3])
    assert transfers == 4 and unsafe == 0
    print('PASS quota routing through a balanced part and edgeless components')
    for family in FAMILIES:
        graph = Graph.load(args.graph_root / f'{family}-fixture.npz')
        owners = graph_growth(graph.offsets, graph.neighbours, 4)
        masks, first = halo_masks(graph.offsets, graph.neighbours, owners, 3)
        for part in range(4):
            bit = np.uint64(1) << np.uint64(part)
            assert set(np.flatnonzero(masks & bit)) == reference_halo(graph, owners, part, 3)
            assert set(np.flatnonzero(first & bit)) == reference_halo(graph, owners, part, 1)
        visited = np.zeros(graph.cells, np.int64)
        goals = visited.copy()
        queue = np.empty(graph.cells, np.int32)
        sequence = np.zeros(1, np.int64)
        for cell in range(graph.cells):
            reference = reference_removable(graph, owners, cell)
            full = removable(graph.offsets, graph.neighbours, owners, cell,
                             visited, goals, queue, sequence, graph.cells)
            bounded = removable(graph.offsets, graph.neighbours, owners, cell,
                                visited, goals, queue, sequence, 64)
            assert full == reference
            assert not bounded or reference
        print(f'PASS {family}: independent 1/3-layer BFS and removal witnesses ({graph.cells} cells)')


if __name__ == '__main__':
    main()
