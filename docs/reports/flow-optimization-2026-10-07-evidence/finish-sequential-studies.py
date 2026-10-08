#!/usr/bin/env python3
"""Finish the remaining studies in order, with one numerical task at a time.

Wait for the resumed PISO study to release its lock before launching anything.
Never continue to another study after incomplete/invalid work. Repetition and
field-recovery evidence remain distinct from the main strong-scaling curve.
"""
from pathlib import Path
import fcntl
import json
import os
import copy
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[2]
BASE = Path('/home/midway/babelsim-flow-opt-20261007')
DRIVER = REPOSITORY/'tools/benchmark_solver_optimization.py'


def execute(label, command):
    print('START', label, flush=True)
    (ROOT/'sequential-campaign-state.json').write_text(json.dumps(dict(phase=label,command=command),indent=2)+'\n')
    with (ROOT/(label+'-driver.log')).open('a') as log:
        result = subprocess.run(command,cwd=REPOSITORY,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'{label} exited {result.returncode}; inspect preserved logs')
    print('FINISHED', label, flush=True)


def study(label, cases, ranks, repeat, *, throughput=False):
    command = [sys.executable,'-u',str(DRIVER),
        '--variant','reference='+str(BASE/'bin/babelsim-memory-partition-normalized'),
        '--variant','optimized='+str(BASE/'bin/babelsim-amg-krylov-policy-cf')]
    for case in cases:
        input_name = 'simple384' if throughput else case+'-piso10'
        command += ['--case',case+'='+str(BASE/'strong-inputs'/input_name),
                    '--throughput' if throughput else '--accepted-steps',case+'='+('384' if throughput else '10')]
    command += ['--ranks',ranks,'--repeat',str(repeat),'--warmup','0','--timeout','7200',
        '--output',str(ROOT/label),'--work-root',str(BASE/(label+'-work')),
        '--field-archive-root',str(BASE/'field-archive'/label)]
    if (ROOT/label/'study.json').exists():
        command += ['--resume']
    execute(label,command)
    records = json.loads((ROOT/label/'records.json').read_text())
    assert all(r['validWork'] and r['measurementValid'] for r in records)
    assert all(r['sustainedSample'] for r in records), 'short sample retained; increase fixed work for every rank/version before publishing'


def main():
    for key in ('LD_PRELOAD','PETSC_OPTIONS'):
        os.environ.pop(key,None)
    with open('/tmp/babelsim-solver-benchmark.lock','w') as lock:
        print('WAITING FOR CURRENT PISO STUDY',flush=True)
        fcntl.flock(lock,fcntl.LOCK_EX)
        records = json.loads((ROOT/'amg-strong-scaling/records.json').read_text())
        assert len(records)==16 and all(r['validWork'] and r['sustainedSample'] for r in records)
    execute('strong-field-recovery',[sys.executable,'-u',str(ROOT/'recover-strong-fields.py')])
    # Preselect the complete post-restart epoch by availability, never by speed.
    # The nine old fields were lost, so all nine associated runs are repeated.
    # Newly completed slots and these recovery runs form one full 16-point set.
    main_records = json.loads((ROOT/'amg-strong-scaling/records.json').read_text())
    recovery = json.loads((ROOT/'strong-field-recovery/records.json').read_text())
    replacements = {(r['variant'],r['caseLabel'],r['ranks']):r for r in recovery}
    fresh=[]
    for old in main_records:
        key = old['variant'],old['caseLabel'],old['ranks']
        row = copy.deepcopy(replacements.get(key,old))
        assert row['fieldArchive']['status']=='retained'
        row['originStudy'] = 'strong-field-recovery' if key in replacements else 'amg-strong-scaling'
        fresh.append(row)
    assert len(fresh)==16 and all(r['sustainedSample'] for r in fresh)
    metadata = json.loads((ROOT/'amg-strong-scaling/study.json').read_text())
    metadata.update(selection='all 16 measurements after the workspace restart; selected before observing recovery timings',
        originalStudy='amg-strong-scaling',fieldRecoveryStudy='strong-field-recovery',
        sourceInputAndBinaryHashesUnchanged=True)
    aggregate=ROOT/'amg-strong-post-restart'
    aggregate.mkdir(exist_ok=True)
    (aggregate/'study.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (aggregate/'records.json').write_text(json.dumps(fresh,indent=2)+'\n')
    study('amg-strong-repeat-12',['prism','voronoi'],'12',2)
    study('simple-strong-scaling',['simple'],'1,4,8,12',1,throughput=True)
    study('simple-repeat-12',['simple'],'12',2,throughput=True)
    pairs=[]
    for case in ('prism','voronoi'):
        directories={}
        for variant in ('reference','optimized'):
            for ranks in (1,4,8,12):
                relative = Path(variant)/case/f'ranks-{ranks}'/'repeat-01'
                index = ROOT/'amg-strong-scaling'/relative/'fields.json'
                if not index.exists() or json.loads(index.read_text())['status']!='retained':
                    index = ROOT/'strong-field-recovery'/relative/'fields.json'
                fields = json.loads(index.read_text())
                assert fields['status']=='retained' and fields['files']
                times = {Path(entry['path']).parent.parent for entry in fields['files']}
                # Production export also maintains a latest-state alias below
                # the run root. Select the named accepted physical time set.
                named_times={p for p in times if p.name=='0.1'}
                assert len(named_times)==1
                final_time=named_times.pop()
                directories[variant,ranks]=str(Path(fields['root'])/final_time)
        for ranks in (1,4,8,12):
            pairs.append(dict(label=f'{case}-versions-{ranks}',reference=directories['reference',ranks],
                              candidate=directories['optimized',ranks]))
        for ranks in (4,8,12):
            pairs.append(dict(label=f'{case}-optimized-1-{ranks}',reference=directories['optimized',1],
                              candidate=directories['optimized',ranks]))
    plan=ROOT/'strong-field-pairs.json'
    plan.write_text(json.dumps(pairs,indent=2)+'\n')
    execute('strong-field-comparison',[sys.executable,'-u',str(ROOT/'compare-flow-fields.py'),
        '--pairs',str(plan),'--output',str(ROOT/'strong-field-comparisons.json')])
    for script in ('analyze-strong-scaling.py','analyze-simple-throughput.py','render-strong-report.py'):
        execute(script.removesuffix('.py'),[sys.executable,str(ROOT/script)])
    (ROOT/'sequential-campaign-state.json').write_text(json.dumps(dict(phase='completed'),indent=2)+'\n')
    print('ALL NUMERICAL STUDIES AND REPORT GENERATION COMPLETED',flush=True)


if __name__ == '__main__':
    main()
