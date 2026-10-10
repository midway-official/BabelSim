#!/usr/bin/env python3
"""Verify measured binaries, tested source, final data and delivery artifacts.

Run after all numerical work and report generation. Existing stage snapshots
remain immutable; the delivery snapshot includes later measurement fixes.
"""
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess

ROOT=Path(__file__).resolve().parent
REPOSITORY=ROOT.parents[2]
BASE=Path('/home/midway/babelsim-flow-opt-20261007')
BASE_COMMIT='83d7b72757f7f9af2db1a01a097389f86c296117'
SHAS=dict(reference='7a29412a574f912f6d67c7b3514b7aa15c60b854c222a455392ed977cb61a6b2',
    optimized='2114c30d7632f2fa500822cfafe2a134c69a5e280089bc312243ff3be2e7c53c')


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1048576),b''):
            h.update(chunk)
    return h.hexdigest()


def item(path):
    return dict(path=str(path.relative_to(REPOSITORY)),sha256=sha(path),bytes=path.stat().st_size)


def verify_scaling(summary, studies):
    """Independently rebuild median denominators from distinct raw runs."""
    runs=[]
    identities=set()
    case_hashes=summary['caseHashes'] if 'caseHashes' in summary else {'simple':summary['caseHash']}
    for study in studies:
        metadata=json.loads((ROOT/study/'study.json').read_text())
        assert {name:value['sha256'] for name,value in metadata['variants'].items()}==SHAS
        for case,value in metadata['cases'].items():
            assert value['manifest']['sha256']==case_hashes[case],(study,case)
        records=json.loads((ROOT/study/'records.json').read_text())
        for run in records:
            if run['warmup']:
                continue
            assert run['validWork'] and run['measurementValid'] and run['sustainedSample']
            assert math.isfinite(run['computeSeconds']) and run['computeSeconds']>=60
            assert run['returncode']==(2 if run['caseLabel']=='simple' else 0)
            assert run['status']==('maxIterations' if run['caseLabel']=='simple' else 'converged')
            reports=run['performance']
            assert len(reports)==run['ranks']
            identity=str(Path(reports[0]['path']).parent.resolve())
            assert identity not in identities, ('duplicate measured run',identity)
            identities.add(identity)
            assert abs(run['computeSeconds']-max(p['iterationComputeSeconds'] for p in reports))<1e-9
            for report in reports:
                current=json.loads(Path(report['path']).read_text())
                assert current['iterationComputeSeconds']==report['iterationComputeSeconds']
                assert current['reportedIterationIntervals']==report['reportedIterationIntervals']
                assert current['status']==run['status']
            runs.append(run)
    groups={}
    for run in runs:
        groups.setdefault((run['caseLabel'],run['variant'],run['ranks']),[]).append(run)
    rows=summary['scaling']
    expected={(row.get('case','simple'),row['variant'],row['ranks']) for row in rows}
    assert len(expected)==len(rows) and set(groups)==expected
    for row in rows:
        case=row.get('case','simple')
        group=groups[case,row['variant'],row['ranks']]
        times=[run['computeSeconds'] for run in group]
        serial=statistics.median(run['computeSeconds'] for run in groups[case,row['variant'],1])
        median=statistics.median(times)
        assert len(times)==row['samples']==(3 if row['ranks'] in (1,12) else 1)
        values=dict(computeMedian=median,computeMin=min(times),computeMax=max(times),
                    serialMedian=serial,speedup=serial/median,
                    efficiency=serial/(row['ranks']*median))
        for name,value in values.items():
            assert math.isclose(row[name],value,rel_tol=1e-12,abs_tol=1e-12),(case,row['variant'],row['ranks'],name)
    return dict(distinctRuns=len(runs),groups=len(groups),medianDenominatorsRecomputed=True,
                currentPerformanceFilesMatch=True)


def fields_digest(directory):
    """Hash current archived contents in the comparison's manifest order."""
    directory=Path(directory)
    ranks=sorted(directory.glob('rank-*'))
    assert ranks,directory
    entries=[]
    for rank in ranks:
        metadata=rank/'metadata.bs'
        entries.append((str(metadata.relative_to(directory)),sha(metadata)))
        names=sorted(tokens[1] for line in metadata.read_text().splitlines()
                     if (tokens:=line.split()) and tokens[0]=='field')
        assert names and len(names)==len(set(names))
        for name in names:
            assert Path(name).name==name
            path=rank/(name+'.csv')
            entries.append((str(path.relative_to(directory)),sha(path)))
    return hashlib.sha256(json.dumps(entries,sort_keys=True).encode()).hexdigest()


def main():
    previous=json.loads((ROOT/'final-source-amg/manifest.json').read_text())
    assert previous['binarySha256']==SHAS['optimized']
    binaries={'reference':BASE/'bin/babelsim-memory-partition-normalized',
              'optimized':BASE/'bin/babelsim-amg-krylov-policy-cf'}
    for variant,path in binaries.items():
        assert sha(path)==SHAS[variant],path
    assert sha(REPOSITORY/'build-petsc/babelsim-solve')==SHAS['optimized']
    # Production and test source must still be exactly the source already tested.
    for entry in previous['files']:
        name=entry['path']
        if name=='Makefile' or name.startswith(('include/','src/','tests/')):
            assert sha(REPOSITORY/name)==entry['sha256'],name
    test=previous['fullTest']
    assert test['exitCode']==0 and sha(REPOSITORY/test['log'])==test['logSha256']
    physics=subprocess.check_output(['git','diff','--name-only',BASE_COMMIT,'--','src/physics'],
                                   cwd=REPOSITORY,text=True)
    assert not physics.strip()
    scaling_checks={}
    studies={'strong-summary.json':('amg-strong-post-restart','amg-voronoi20-strong-scaling',
                                   'amg-prism-repeat','amg-voronoi20-repeat'),
             'simple-summary.json':('simple-strong-scaling','simple-repeat')}
    for name in studies:
        data=json.loads((ROOT/name).read_text())
        assert data['binarySha256']==SHAS
        for row in data['scaling']:
            assert row['samples']==(3 if row['ranks'] in (1,12) else 1)
            assert row['computeMin']>=60
            assert abs(row['efficiency']-row['serialMedian']/(row['ranks']*row['computeMedian']))<1e-12
        scaling_checks[name]=verify_scaling(data,studies[name])
    checks=json.loads((ROOT/'strong-field-comparisons.json').read_text())
    assert len(checks)==14 and all(c['passed'] for c in checks)
    assert sum(c['label'].startswith('prism-') and c['cells']==1024000 and c['time']==.1 for c in checks)==7
    assert sum(c['label'].startswith('voronoi-') and c['cells']==512000 and c['time']==.2 for c in checks)==7
    field_hashes={}
    for check in checks:
        for side in ('reference','candidate'):
            directory=check[side]
            if directory not in field_hashes:
                field_hashes[directory]=fields_digest(directory)
            assert field_hashes[directory]==check[side+'FieldsSha256'],(check['label'],side)
    assert json.loads((ROOT/'extended-campaign-state.json').read_text())['phase']=='completed'
    target=ROOT/'delivery-source'
    target.mkdir(exist_ok=True)
    paths={REPOSITORY/e['path'] for e in previous['files']}
    paths.update(ROOT.glob('*.py'))
    paths.add(ROOT/'cg-single-reduction-preload.c')
    for path in sorted(paths):
        destination=target/path.relative_to(REPOSITORY)
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,destination)
        assert sha(path)==sha(destination)
    diff=subprocess.check_output(['git','diff',BASE_COMMIT,'--',
        *[str(p.relative_to(REPOSITORY)) for p in sorted(paths)]],cwd=REPOSITORY)
    (target/'diff-from-baseline.patch').write_bytes(diff)
    artifact_names=['strong-summary.json','simple-summary.json','strong-field-comparisons.json',
        'strong-field-pairs.json','strong-scaling.csv','strong-comparison.csv','simple-scaling.csv',
        'simple-comparison.csv','voronoi20-input.json','simple-throughput-input.json',
        'recovery-2026-10-08.json','recovered-build-recipe.json',
        'checkpoint-git-push.json',
        'framework-scope-review.md',
        'benchmark-exit-contract-audit.json',
        'cg-single-reduction-screen-summary.json','cg-single-reduction-decision.json',
        'cg-single-reduction-screen-overlay.json','amg-strong-profiles/v12-current-profile-summary.json']
    artifacts=[ROOT/name for name in artifact_names]
    artifacts.append(ROOT.parent/'flow-optimization-2026-10-07.md')
    for prefix in ('strong-scaling','strong-memory','simple-scaling'):
        artifacts.extend(ROOT/(prefix+'.'+suffix) for suffix in ('png','svg','pdf'))
    for study in ('amg-strong-post-restart','amg-voronoi20-strong-scaling','amg-prism-repeat',
                  'amg-voronoi20-repeat','simple-strong-scaling','simple-repeat'):
        artifacts.extend(ROOT/study/name for name in ('study.json','records.json'))
    manifest=dict(createdUtc=datetime.now(timezone.utc).isoformat(),
        git=dict(branch=subprocess.check_output(['git','branch','--show-current'],cwd=REPOSITORY,text=True).strip(),
                 commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPOSITORY,text=True).strip(),
                 sourceBaselineCommit=BASE_COMMIT,
                 trackedPaths=[str(p.relative_to(REPOSITORY)) for p in sorted(paths)]),
        binarySha256=SHAS,binaries={k:str(v) for k,v in binaries.items()},
        physicsDiffEmpty=True,productionAndTestSourceMatchPassedFullRegression=True,
        fullTest=test,build=json.loads((ROOT/'recovered-build-recipe.json').read_text()),
        scalingDataAudit=scaling_checks,currentFieldDirectoriesRehashed=len(field_hashes),
        files=[item(p) for p in sorted(paths)],artifacts=[item(p) for p in artifacts],
        finalFieldPairsPassed=14,previousImmutableSnapshot=str(ROOT/'final-source-amg/manifest.json'))
    (target/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status='passed',files=len(paths),artifacts=len(artifacts),
        manifest=str(target/'manifest.json'),binarySha256=SHAS),indent=2))


if __name__=='__main__':
    main()
