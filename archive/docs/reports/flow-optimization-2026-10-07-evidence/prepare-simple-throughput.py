#!/usr/bin/env python3
"""Prepare engineering SIMPLE fixed-work inputs without changing discretization.

Use an existing validated 500,094-cell mesh. Mesh construction is not repeated.
Run preparation before the numerical study, never alongside a timed solver.
The tight outer tolerances intentionally reach the iteration cap; this study
measures sustained SIMPLE throughput, not time to physical steady convergence.
"""
import json
from pathlib import Path
import shutil
import sys

REPOSITORY = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPOSITORY/'tools'))
import benchmark_backend as bench

SOURCE = Path('/tmp/babelsim-flow-opt-20261007/formal-inputs/simple128')
TARGET = Path('/tmp/babelsim-flow-opt-20261007/strong-inputs/simple384')


def copy_file(source, target):
    if Path(source).suffix == '.mesh':
        Path(target).hardlink_to(source)
    else:
        shutil.copy2(source, target)
    return str(target)


def main():
    assert not TARGET.exists(), TARGET
    shutil.copytree(SOURCE, TARGET, ignore=shutil.ignore_patterns('results','post'),
                    copy_function=copy_file)
    path = TARGET/'numerics/solution.bs'
    old = path.read_text()
    assert old.count('maxIterations 128\n') == 1
    new = old.replace('maxIterations 128\n', 'maxIterations 384\n', 1)
    path.write_text(new)
    assert new.replace('maxIterations 384\n','maxIterations 128\n',1) == old
    metadata = dict(case=str(TARGET), source=str(SOURCE), cells=500094, faces=1266111,
        outerIterations=384, measuredIntervals=383, status='expected maxIterations',
        purpose='SIMPLE throughput; not time to physical steady convergence',
        unchangedExceptOuterIterationCap=True,
        inputSha256=bench.case_manifest(TARGET)['sha256'],
        meshSha256=bench.sha256(TARGET/'mesh/flow.mesh'))
    output=Path(__file__).resolve().parent/'simple-throughput-input.json'
    output.write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata),flush=True)


if __name__ == '__main__':
    main()
