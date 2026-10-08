#!/usr/bin/env python3
"""Verify measured binaries, tested source, final data and delivery artifacts.

Run after all numerical work and report generation. Existing stage snapshots
remain immutable; the delivery snapshot includes later measurement fixes.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
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
    for name in ('strong-summary.json','simple-summary.json'):
        data=json.loads((ROOT/name).read_text())
        assert data['binarySha256']==SHAS
        for row in data['scaling']:
            assert row['samples']==(3 if row['ranks'] in (1,12) else 1)
            assert row['computeMin']>=60
            assert abs(row['efficiency']-row['serialMedian']/(row['ranks']*row['computeMedian']))<1e-12
    checks=json.loads((ROOT/'strong-field-comparisons.json').read_text())
    assert len(checks)==14 and all(c['passed'] for c in checks)
    assert sum(c['label'].startswith('prism-') and c['cells']==1024000 and c['time']==.1 for c in checks)==7
    assert sum(c['label'].startswith('voronoi-') and c['cells']==512000 and c['time']==.2 for c in checks)==7
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
        files=[item(p) for p in sorted(paths)],artifacts=[item(p) for p in artifacts],
        finalFieldPairsPassed=14,previousImmutableSnapshot=str(ROOT/'final-source-amg/manifest.json'))
    (target/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status='passed',files=len(paths),artifacts=len(artifacts),
        manifest=str(target/'manifest.json'),binarySha256=SHAS),indent=2))


if __name__=='__main__':
    main()
