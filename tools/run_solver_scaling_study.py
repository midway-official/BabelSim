#!/usr/bin/env python3
"""Sequential, interleaved solver-only scaling experiments using the existing harness.

Example: --case heat=/tmp/heat --case candidate=/tmp/candidate --ranks 1,2,4,8
Only solverComputeSeconds determines speedup. Failed solves remain in raw evidence
and are excluded from performance summaries. Warmups are retained but not measured.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import random
import shutil
import statistics

import benchmark_backend as bench


def summarize(runs, output):
    groups = {}
    for name, run in runs:
        if run['warmup'] or run['status'] != 'converged':
            continue
        reports = run['performance']
        if len(reports) != run['ranks']:
            continue
        local_max = lambda key: max(r['local'][key] for r in reports)
        groups.setdefault((name, run['ranks']), []).append({
            'compute': max(r['solverComputeSeconds'] for r in reports),
            'setup': local_max('preconditionerSeconds'),
            'solve': local_max('linearSolveSeconds'),
            'iterations': local_max('krylovIterations'),
            'pcSetups': local_max('preconditionerSetups'),
            'linearSolves': local_max('linearSolves'),
        })
    rows = []
    for (name, rank), samples in sorted(groups.items()):
        times = [s['compute'] for s in samples]
        baseline = groups.get((name, 1))
        t = statistics.median(times)
        speedup = statistics.median(s['compute'] for s in baseline) / t if baseline else None
        rows.append(dict(case=name, ranks=rank, samples=len(samples),
                         computeMedian=t, computeMin=min(times), computeMax=max(times),
                         speedup=speedup, efficiency=speedup / rank if speedup else None,
                         **{k: statistics.median(s[k] for s in samples)
                            for k in ('setup', 'solve', 'iterations', 'pcSetups', 'linearSolves')}))
    (output / 'summary.json').write_text(json.dumps(rows, indent=2) + '\n')
    if rows:
        with (output / 'summary.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', action='append', required=True, help='name=/absolute/case/path')
    parser.add_argument('--ranks', type=bench.parse_ranks, default=[1, 2, 4, 8])
    parser.add_argument('--repeat', type=int, default=3)
    parser.add_argument('--warmup', type=int, default=1)
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--binary', type=Path, default=bench.ROOT / 'build/babelsim-solve')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=20260928)
    args = parser.parse_args()
    if args.repeat < 1 or args.warmup not in (0, 1):
        parser.error('repeat must be positive; warmup must be 0 or 1')
    args.output = args.output.resolve()
    args.binary = args.binary.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    if (args.output / 'study.json').exists():
        parser.error('choose a new evidence directory to preserve previous runs')
    cases = {}
    for item in args.case:
        name, path = item.split('=', 1)
        if not name or '/' in name or name in cases:
            parser.error('case names must be unique path-free labels')
        cases[name] = Path(path).resolve()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[key] = '1'
    # Keep small authoritative dictionaries even if exploratory cases later change.
    for name, case in cases.items():
        for source in case.rglob('*'):
            relative = source.relative_to(case)
            if ('results' not in relative.parts and source.is_file()
                    and source.suffix in ('.bs', '.field', '.json', '.dat')):
                target = args.output / 'inputs' / name / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
    # The harness records binding, commands, outputs, peak RSS, and each rank JSON.
    (args.output / 'study.json').write_text(json.dumps({
        'host': bench.host_metadata(), 'git': bench.git_metadata(),
        'binary': str(args.binary), 'binarySha256': bench.sha256(args.binary),
        'cases': {n: {'path': str(p), 'manifest': bench.case_manifest(p)} for n, p in cases.items()},
        'ranks': args.ranks, 'repeat': args.repeat, 'warmup': args.warmup,
        'seed': args.seed, 'metric': 'max_rank(solverComputeSeconds)',
    }, indent=2) + '\n')
    jobs = [(name, rank) for name in cases for rank in args.ranks]
    rng = random.Random(args.seed)
    runs = []
    for repeat in range(0 if args.warmup else 1, args.repeat + 1):
        order = jobs.copy()
        rng.shuffle(order)
        for name, rank in order:
            run_args = argparse.Namespace(case=cases[name], binary=args.binary,
                                          mode='complete', resume=False, timeout=args.timeout)
            print(f'RUN case={name} ranks={rank} repeat={repeat}', flush=True)
            result = bench.run_one(run_args, args.output / name, rank, repeat, warmup=repeat == 0)
            runs.append((name, result))
            summarize(runs, args.output)
            times = [r['solverComputeSeconds'] for r in result['performance']]
            print(f"DONE status={result['status']} compute={max(times) if times else None}", flush=True)
    failures = [r for _, r in runs if r['status'] != 'converged']
    print(f'Completed {len(runs)} runs; {len(failures)} did not converge. Evidence: {args.output}', flush=True)
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
