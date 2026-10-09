#!/usr/bin/env python3
"""Validate measured runs against Hypre at one rank (or a specified reference rank)."""
import argparse
import json
import math
from pathlib import Path
import re

from compare_parallel_results import read_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('--reference-rank', type=int, default=1)
    args = parser.parse_args()
    meta = json.loads((args.study / 'study.json').read_text())
    references = {}
    comparisons = []
    for name, config in meta['cases'].items():
        family = name.split('-')[0]
        reference_case = meta['cases'][family + '-hypre']['path']
        reference_path = Path(reference_case) / f'results/benchmark-r{args.reference_rank}-n01'
        if family not in references:
            references[family] = read_result(reference_path)
        reference = references[family]
        atol = rtol = 1e-8 if family.startswith('heat') else 5e-6
        for rank in meta['ranks']:
            for repeat in range(1, meta['repeat'] + 1):
                run_dir = args.study / name / f'ranks-{rank}' / f'repeat-{repeat:02d}'
                run = json.loads((run_dir / 'run.json').read_text())
                if run['status'] != 'converged' or run['returncode'] != 0:
                    raise ValueError(f'unsuccessful solve: {run_dir}')
                if len(run['performance']) != rank:
                    raise ValueError(f'incomplete rank evidence: {run_dir}')
                expected_cells = run['performance'][0]['partition']['globalCells']
                candidate_path = Path(config['path']) / f'results/benchmark-r{rank}-n{repeat:02d}'
                candidate = read_result(candidate_path)
                if candidate.keys() != reference.keys():
                    raise ValueError(f'field sets differ: {run_dir}')
                result = dict(case=name, ranks=rank, repeat=repeat, reference=str(reference_path),
                              candidate=str(candidate_path), atol=atol, rtol=rtol, fields={})
                for field, expected in reference.items():
                    actual = candidate[field]
                    if expected.keys() != actual.keys() or len(actual) != expected_cells:
                        raise ValueError(f'global cell coverage differs: {run_dir}/{field}')
                    maximum = maximum_ratio = error2 = ref2 = 0.0
                    for cell, values in expected.items():
                        other = actual[cell]
                        if len(values) != len(other):
                            raise ValueError('component count differs')
                        for a, b in zip(values, other):
                            delta = abs(a - b)
                            maximum = max(maximum, delta)
                            maximum_ratio = max(maximum_ratio, delta / (atol + rtol * max(abs(a), abs(b))))
                            error2 += delta * delta
                            ref2 += a * a
                    result['fields'][field] = dict(cells=len(actual), maxAbs=maximum,
                                                   maxToleranceRatio=maximum_ratio,
                                                   relativeL2=math.sqrt(error2 / ref2) if ref2 else None)
                    if maximum_ratio > 1:
                        raise ValueError(f'field tolerance exceeded: {run_dir}/{field}: {maximum_ratio}')
                expected_time = 1.0 if family.startswith('heat') else 0.05
                for rank_dir in candidate_path.glob('rank-*'):
                    pairs = [line.split() for line in (rank_dir / 'metadata.bs').read_text().splitlines()]
                    saved_time = next(float(p[1]) for p in pairs if len(p) == 2 and p[0] == 'time')
                    if abs(saved_time - expected_time) > 1e-12:
                        raise ValueError(f'wrong final time: {rank_dir}')
                if family.startswith(('piso', 'naca')):
                    lines = [line for line in (run_dir / 'stdout.log').read_text().splitlines()
                             if re.match(r'PISO \d+ time=', line)]
                    if len(lines) != 5 or any('accepted=true' not in l or 'linear=ok' not in l for l in lines):
                        raise ValueError(f'invalid time-step acceptance: {run_dir}')
                    result['maxCoupling'] = max(float(re.search(r' rCoupling=([^ ]+)', l)[1]) for l in lines)
                    result['maxMass'] = max(float(re.search(r' mass=([^ ]+)', l)[1]) for l in lines)
                    result['correctors'] = [int(re.search(r' correctors=(\d+)', l)[1]) for l in lines]
                    mass_limit = 1e-8 if family == 'naca' else 1e-7
                    if result['maxCoupling'] > 1e-5 or result['maxMass'] > mass_limit:
                        raise ValueError(f'coupling or mass tolerance exceeded: {run_dir}')
                comparisons.append(result)
                print(f'PASS {name} ranks={rank} repeat={repeat}', flush=True)
    output = args.study / 'validation.json'
    output.write_text(json.dumps(dict(status='passed', measuredRuns=len(comparisons),
                                     comparisons=comparisons), indent=2) + '\n')
    print(f'Validated {len(comparisons)} measured runs: {output}')


if __name__ == '__main__':
    main()
