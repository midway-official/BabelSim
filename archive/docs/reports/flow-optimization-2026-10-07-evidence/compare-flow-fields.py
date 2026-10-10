#!/usr/bin/env python3
"""Validate large exported flow fields, geometry, time and global-id coverage."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def read(directory):
    ranks = sorted(directory.glob('rank-*'))
    assert ranks, directory
    fields, seen, manifests = {}, {}, []
    signature = None
    for rank, path in enumerate(ranks):
        metadata = path / 'metadata.bs'
        lines = [line.split() for line in metadata.read_text().splitlines()]
        values = {tokens[0]: tokens[1:] for tokens in lines if tokens[0] != 'field'}
        names = {tokens[1]: tokens[2:] for tokens in lines if tokens[0] == 'field'}
        assert names and len(names) == sum(tokens[0] == 'field' for tokens in lines)
        n = int(values['global_cell_count'][0])
        current = (float(values['time'][0]), n, tuple(sorted((k, tuple(v)) for k, v in names.items())))
        assert int(values['rank'][0]) == rank and int(values['ranks'][0]) == len(ranks)
        assert signature is None or current == signature
        signature = current
        manifests.append((str(metadata.relative_to(directory)), hashlib.sha256(metadata.read_bytes()).hexdigest()))
        for name in sorted(names):
            csv = path / (name + '.csv')
            with csv.open() as stream:
                header = stream.readline().strip().split(',')
            assert header[:4] == ['global_id', 'x', 'y', 'z']
            assert header[4:] and all(key.startswith('value') for key in header[4:])
            component_counts = {'scalar':1, 'vector':3, 'tensor':9}
            assert names[name][1:] == ['cell'], (name, names[name])
            assert len(header[4:]) == component_counts[names[name][0]], (name, header)
            data = np.loadtxt(csv, delimiter=',', skiprows=1, ndmin=2)
            assert data.shape == (int(values['owned_cells'][0]), len(header))
            assert np.isfinite(data).all()
            ids = data[:, 0].astype(np.int64)
            assert np.array_equal(ids, data[:, 0]) and ((ids >= 0) & (ids < n)).all()
            assert len(np.unique(ids)) == len(ids)
            if name not in fields:
                fields[name] = np.empty((n, len(header)-1))
                seen[name] = np.zeros(n, dtype=bool)
            assert fields[name].shape[1] == len(header)-1 and not seen[name][ids].any()
            fields[name][ids] = data[:, 1:]
            seen[name][ids] = True
            manifests.append((str(csv.relative_to(directory)), hashlib.sha256(csv.read_bytes()).hexdigest()))
    assert all(mask.all() for mask in seen.values())
    digest = hashlib.sha256(json.dumps(manifests, sort_keys=True).encode()).hexdigest()
    return signature, fields, digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairs', type=Path, required=True,
                        help='JSON array of label, reference, candidate directories')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--atol', type=float, default=5e-6)
    parser.add_argument('--rtol', type=float, default=5e-6)
    args = parser.parse_args()
    cache, results = {}, []
    for pair in json.loads(args.pairs.read_text()):
        items = []
        for key in ('reference', 'candidate'):
            directory = Path(pair[key]).resolve()
            if directory not in cache:
                cache[directory] = read(directory)
            items.append(cache[directory])
        (a, lhs, lhash), (b, rhs, rhash) = items
        assert a == b, (pair['label'], a, b)
        row = dict(pair, atol=args.atol, rtol=args.rtol, time=a[0], cells=a[1],
                   referenceFieldsSha256=lhash, candidateFieldsSha256=rhash,
                   fields={}, passed=True)
        for name in lhs:
            assert lhs[name].shape == rhs[name].shape
            geometry_error = np.max(np.abs(lhs[name][:, :3]-rhs[name][:, :3]))
            assert geometry_error <= 1e-11, (pair['label'], name, geometry_error)
            lvalues, rvalues = lhs[name][:, 3:], rhs[name][:, 3:]
            error = np.abs(lvalues-rvalues)
            limit = args.atol + args.rtol*np.maximum(np.abs(lvalues), np.abs(rvalues))
            count = int((error > limit).sum())
            row['fields'][name] = dict(components=lvalues.shape[1], maxAbsoluteDifference=float(error.max()),
                                     maxGeometryDifference=float(geometry_error), mismatchedComponents=count)
            row['passed'] &= count == 0
        results.append(row)
        args.output.write_text(json.dumps(results, indent=2)+'\n')
        print(json.dumps(row), flush=True)
        assert row['passed'], pair['label']


if __name__ == '__main__':
    main()
