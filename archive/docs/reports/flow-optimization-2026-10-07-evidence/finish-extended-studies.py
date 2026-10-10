#!/usr/bin/env python3
"""Finish fixed-work scaling and repeat the serial efficiency denominator.

Prism uses 10 accepted PISO steps; Voronoi uses 20 at every rank/version.
The nine lost field exports were recovered in an earlier sequential campaign.
The short Voronoi-10 point stays in that campaign, never in this final curve.
"""
from pathlib import Path
import copy
from datetime import datetime, timezone
import fcntl
import json
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[2]
BASE = Path('/home/midway/babelsim-flow-opt-20261007')
DRIVER = REPOSITORY/'tools/benchmark_solver_optimization.py'
STEPS = {'prism':10, 'voronoi':20}


def execute(label, command):
    print('START',label,flush=True)
    state = dict(phase=label, command=command,
                 updatedUtc=datetime.now(timezone.utc).isoformat())
    (ROOT/'extended-campaign-state.json').write_text(json.dumps(state,indent=2)+'\n')
    with (ROOT/(label+'-driver.log')).open('a') as log:
        result = subprocess.run(command,cwd=REPOSITORY,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'{label} exited {result.returncode}; inspect retained evidence')
    print('FINISHED',label,flush=True)


def study(label, case, ranks, repeat, *, throughput=False):
    input_name = 'simple384' if throughput else f'{case}-piso{STEPS[case]}'
    command = [sys.executable,'-u',str(DRIVER),
        '--variant','reference='+str(BASE/'bin/babelsim-memory-partition-normalized'),
        '--variant','optimized='+str(BASE/'bin/babelsim-amg-krylov-policy-cf'),
        '--case',case+'='+str(BASE/'strong-inputs'/input_name),
        '--throughput' if throughput else '--accepted-steps',
        case+'='+str(384 if throughput else STEPS[case]),
        '--ranks',ranks,'--repeat',str(repeat),'--warmup','0',
        '--timeout','7200','--minimum-compute','60','--seed','20261007',
        '--output',str(ROOT/label),'--work-root',str(BASE/(label+'-work')),
        '--field-archive-root',str(BASE/'field-archive'/label)]
    if (ROOT/label/'study.json').exists():
        command += ['--resume']
    execute(label,command)
    records = json.loads((ROOT/label/'records.json').read_text())
    assert len(records)==2*len(ranks.split(','))*repeat
    assert all(r['validWork'] and r['measurementValid'] for r in records)
    assert all(r['sustainedSample'] for r in records), \
        'short sample retained; increase work uniformly before publishing'


def prism_aggregate():
    original = json.loads((ROOT/'amg-strong-scaling/records.json').read_text())
    recovery = json.loads((ROOT/'strong-field-recovery/records.json').read_text())
    replacements = {(r['variant'],r['caseLabel'],r['ranks']):r for r in recovery}
    fresh=[]
    for old in original:
        if old['caseLabel']!='prism':
            continue
        key = old['variant'],old['caseLabel'],old['ranks']
        row = copy.deepcopy(replacements.get(key,old))
        assert row['fieldArchive']['status']=='retained'
        row['originStudy'] = 'strong-field-recovery' if key in replacements else 'amg-strong-scaling'
        fresh.append(row)
    assert len(fresh)==8 and all(r['validWork'] and r['sustainedSample'] for r in fresh)
    metadata = json.loads((ROOT/'amg-strong-scaling/study.json').read_text())
    metadata['cases']={'prism':metadata['cases']['prism']}
    metadata['acceptedSteps']={'prism':10}
    metadata.update(selection='all eight post-restart prism measurements selected by field availability, never by speed',
        originalStudy='amg-strong-scaling',fieldRecoveryStudy='strong-field-recovery',
        excludedCase='Voronoi uses a new fixed-work 20-step curve because a 10-step run was shorter than 60 seconds',
        sourceInputAndBinaryHashesUnchanged=True,
        recoveredBuildRecipe=str(ROOT/'recovered-build-recipe.json'))
    aggregate=ROOT/'amg-strong-post-restart'
    aggregate.mkdir(exist_ok=True)
    (aggregate/'study.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (aggregate/'records.json').write_text(json.dumps(fresh,indent=2)+'\n')


def field_plan():
    pairs=[]
    for case, study_name in (('prism','amg-strong-post-restart'),
                             ('voronoi','amg-voronoi20-strong-scaling')):
        records=json.loads((ROOT/study_name/'records.json').read_text())
        directories={}
        final_time='0.1' if case=='prism' else '0.2'
        for row in records:
            fields=row['fieldArchive']
            assert fields['status']=='retained' and fields['files']
            times={Path(entry['path']).parent.parent for entry in fields['files']}
            # Ignore the production latest-state alias below the run root.
            named_times={p for p in times if p.name==final_time}
            assert len(named_times)==1
            key=row['variant'],row['ranks']
            assert key not in directories
            directories[key]=str(Path(fields['root'])/named_times.pop())
        assert len(directories)==8
        for ranks in (1,4,8,12):
            pairs.append(dict(label=f'{case}-versions-{ranks}',reference=directories['reference',ranks],
                              candidate=directories['optimized',ranks]))
        for ranks in (4,8,12):
            pairs.append(dict(label=f'{case}-optimized-1-{ranks}',reference=directories['optimized',1],
                              candidate=directories['optimized',ranks]))
    plan=ROOT/'strong-field-pairs.json'
    plan.write_text(json.dumps(pairs,indent=2)+'\n')
    return plan


def main():
    for key in ('LD_PRELOAD','PETSC_OPTIONS'):
        os.environ.pop(key,None)
    os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with open('/tmp/babelsim-solver-benchmark.lock','w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
    prism_aggregate()
    # Every group runs sequentially; repetitions shuffle both binaries and ranks.
    # Repeat T1 as well as T12, so a slow isolated serial run cannot inflate E12.
    study('amg-voronoi20-strong-scaling','voronoi','1,4,8,12',1)
    plan=field_plan()
    execute('strong-field-comparison',[sys.executable,'-u',str(ROOT/'compare-flow-fields.py'),
        '--pairs',str(plan),'--output',str(ROOT/'strong-field-comparisons.json')])
    study('simple-strong-scaling','simple','1,4,8,12',1,throughput=True)
    study('amg-prism-repeat','prism','1,12',2)
    study('amg-voronoi20-repeat','voronoi','1,12',2)
    study('simple-repeat','simple','1,12',2,throughput=True)
    for script in ('analyze-strong-scaling.py','analyze-simple-throughput.py','render-strong-report.py'):
        execute(script.removesuffix('.py'),[sys.executable,str(ROOT/script)])
    (ROOT/'extended-campaign-state.json').write_text(json.dumps(dict(
        phase='completed',updatedUtc=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
    print('ALL FIXED-WORK STUDIES AND REPORT GENERATION COMPLETED',flush=True)


if __name__=='__main__':
    main()
