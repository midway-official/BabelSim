#!/usr/bin/env python3
"""Audit paired screening; never substitute these for a production curve."""
from pathlib import Path
import hashlib
import json
import math
import statistics
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
STUDY=ROOT/'cg-single-reduction-screen'


def time_directory(run):
    archive=run['fieldArchive']
    assert archive['status']=='retained'
    directories={Path(item['path']).parent.parent for item in archive['files']}
    named={p for p in directories if p.name=='0.1'}
    assert len(named)==1
    return str(Path(archive['root'])/named.pop())


def main():
    records=json.loads((STUDY/'records.json').read_text())
    assert len(records)==6
    by_key={(r['variant'],r['repeat']):r for r in records}
    assert len(by_key)==6
    for r in records:
        assert r['status']=='converged' and r['returncode']==0
        assert r['validWork'] and r['measurementValid'] and r['sustainedSample']
        assert r['workUnits']==10 and r['ranks']==12
        assert len(r['performance'])==12
        assert all(p['reportedIterationIntervals']==9 for p in r['performance'])
    pairs=[]
    comparisons=[]
    for repeat in (1,2,3):
        baseline=by_key['reference',repeat]
        candidate=by_key['optimized',repeat]
        pairs.append(dict(label=f'cg-single-reduction-{repeat}',
            reference=time_directory(baseline),candidate=time_directory(candidate)))
        comparisons.append(dict(repeat=repeat,baselineSeconds=baseline['computeSeconds'],
            singleReductionSeconds=candidate['computeSeconds'],
            speedup=baseline['computeSeconds']/candidate['computeSeconds'],
            baselineKrylov=baseline['performance'][0]['krylovIterations'],
            singleReductionKrylov=candidate['performance'][0]['krylovIterations'],
            baselinePeakPssKiB=baseline['memory']['peakTreePssKiB'],
            singleReductionPeakPssKiB=candidate['memory']['peakTreePssKiB']))
    plan=STUDY/'field-pairs.json'
    plan.write_text(json.dumps(pairs,indent=2)+'\n')
    output=STUDY/'field-comparisons.json'
    subprocess.run([sys.executable,str(ROOT/'compare-flow-fields.py'),
        '--pairs',str(plan),'--output',str(output)],check=True)
    checks=json.loads(output.read_text())
    assert len(checks)==3 and all(c['passed'] for c in checks)
    speedups=[c['speedup'] for c in comparisons]
    geo=math.exp(statistics.mean(math.log(x) for x in speedups))
    data=dict(purpose='paired diagnostic screening; not a production scaling curve',
        baseline='unchanged production v12',candidate='same v12 with public CG single-reduction overlay',
        paired=comparisons,pairedGeometricMeanSpeedup=geo,
        baselineMedianSeconds=statistics.median(c['baselineSeconds'] for c in comparisons),
        singleReductionMedianSeconds=statistics.median(c['singleReductionSeconds'] for c in comparisons),
        memoryMetric='whole-lifetime process-tree PSS; name-matched peakSolverRssKiB excluded because the diagnostic wrapper execs a differently named binary',
        allFieldsPassed=True,
        promising=geo>=1.05 and min(speedups)>1.0,
        screenCriterion='at least 5 percent geometric paired gain and positive gain in all three pairs; inspect variability before production adoption',
        overlayManifestSha256=hashlib.sha256((ROOT/'cg-single-reduction-screen-overlay.json').read_bytes()).hexdigest())
    (ROOT/'cg-single-reduction-screen-summary.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data),flush=True)


if __name__=='__main__':
    main()
