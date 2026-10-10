#!/usr/bin/env python3
"""Recompute only the final engineering fields lost with the old /tmp tree.

Keep the original performance evidence intact. Recovery runs use the exact same
binary bytes and input checksums, and execute under the benchmark advisory lock.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2]/'tools'))
import benchmark_backend as bench
from benchmark_solver_optimization import archive_fields

BASE = Path('/home/midway/babelsim-flow-opt-20261007')
SHAS = {'reference':'7a29412a574f912f6d67c7b3514b7aa15c60b854c222a455392ed977cb61a6b2',
        'optimized':'2114c30d7632f2fa500822cfafe2a134c69a5e280089bc312243ff3be2e7c53c'}
BINARIES = {'reference':BASE/'bin/babelsim-memory-partition-normalized',
            'optimized':BASE/'bin/babelsim-amg-krylov-policy-cf'}


def main():
    study = ROOT/'amg-strong-scaling'
    metadata = json.loads((study/'study.json').read_text())
    output = ROOT/'strong-field-recovery'
    output.mkdir(exist_ok=True)
    for variant, binary in BINARIES.items():
        assert bench.sha256(binary) == SHAS[variant]
    with open('/tmp/babelsim-solver-benchmark.lock','w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):
            os.environ[key] = '1'
        assert not os.environ.get('PETSC_OPTIONS') and not os.environ.get('LD_PRELOAD')
        plans, records = [], []
        for manifest in sorted(study.glob('*/*/ranks-*/repeat-01/run.json')):
            relative = manifest.parent.relative_to(study)
            index = manifest.with_name('fields.json')
            if index.exists() and json.loads(index.read_text())['status'] == 'retained':
                continue
            variant, case = relative.parts[:2]
            old = json.loads(manifest.read_text())
            assert old['status'] == 'converged'
            source = BASE/'strong-inputs'/f'{case}-piso10'
            assert bench.case_manifest(source) == metadata['cases'][case]['manifest']
            target = BASE/'field-recovery-work'/variant/case
            if not target.exists():
                for entry in metadata['cases'][case]['manifest']['files']:
                    rel = Path(entry['path']);dest = target/rel
                    dest.parent.mkdir(parents=True,exist_ok=True)
                    if rel.suffix in ('.mesh','.npz'):
                        dest.hardlink_to(source/rel)
                    else:
                        shutil.copy2(source/rel,dest)
            assert bench.case_manifest(target) == metadata['cases'][case]['manifest']
            plans.append(dict(variant=variant,case=case,ranks=old['ranks'],
                              originalRun=str(manifest),inputSha256=metadata['cases'][case]['manifest']['sha256']))
            (output/'plan.json').write_text(json.dumps(plans,indent=2)+'\n')
            print('RECOVER FIELDS',variant,case,old['ranks'],flush=True)
            args = argparse.Namespace(case=target,binary=BINARIES[variant],mode='complete',
                resume=True,timeout=7200,live_log=True,memory_trace=True)
            run = bench.run_one(args,output/variant/case,old['ranks'],1)
            run_directory = output/relative
            lines = re.findall(r'^PISO [0-9]+ .*accepted=true.*converged=true.*',
                               (run_directory/'stdout.log').read_text(),re.MULTILINE)
            assert run['returncode']==0 and run['status']=='converged'
            assert len(lines)==10 and all('linear=ok' in line for line in lines)
            assert len(run['performance'])==old['ranks']
            assert all(p['reportedIterationIntervals']==9 for p in run['performance'])
            run.update(variant=variant,caseLabel=case,originalRun=str(manifest),
                workUnits=10,validWork=True,measurementValid=True,
                computeSeconds=max(p['iterationComputeSeconds'] for p in run['performance']))
            run['sustainedSample'] = run['computeSeconds']>=60
            run['fieldArchive'] = archive_fields(run,BASE/'field-archive/strong-field-recovery'/relative,
                                                 run_directory/'fields.json')
            assert run['fieldArchive']['status']=='retained'
            records.append(run)
            (output/'records.json').write_text(json.dumps(records,indent=2)+'\n')
            print('RECOVERED',variant,case,old['ranks'],run['computeSeconds'],flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
