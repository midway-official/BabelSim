#!/usr/bin/env python3
"""Increase fixed Voronoi work for every variant/rank after a short sample."""
from pathlib import Path
import json
import shutil
import sys

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'tools'))
import benchmark_backend as bench
BASE=Path('/home/midway/babelsim-flow-opt-20261007')


def main():
    source=BASE/'strong-inputs/voronoi-piso10'
    target=BASE/'strong-inputs/voronoi-piso20'
    original=bench.case_manifest(source)
    assert not target.exists(),target
    for entry in original['files']:
        rel=Path(entry['path']);dest=target/rel
        dest.parent.mkdir(parents=True,exist_ok=True)
        if rel.suffix in ('.mesh','.npz'):
            dest.hardlink_to(source/rel)
        else:
            shutil.copy2(source/rel,dest)
    assert (target/'control.bs').read_text()=='startTime 0\nendTime 0.1\ndeltaT 0.01\n'
    (target/'control.bs').write_text('startTime 0\nendTime 0.2\ndeltaT 0.01\n')
    output=target/'output.bs';old=output.read_text()
    assert old.count('writeInterval 10\n')==1
    output.write_text(old.replace('writeInterval 10\n','writeInterval 20\n',1))
    generation=target/'generation.json';data=json.loads(generation.read_text())
    data.update(acceptedSteps=20,measuredIntervals=19,parentCase=str(source),
        parentInputSha256=original['sha256'],
        reason='55.395-second optimized 12-rank sample at 10 steps; fixed work increased for all ranks and both versions')
    generation.write_text(json.dumps(data,indent=2)+'\n')
    assert bench.sha256(target/'numerics/solution.bs')==bench.sha256(source/'numerics/solution.bs')
    assert bench.sha256(target/'numerics/methods.bs')==bench.sha256(source/'numerics/methods.bs')
    manifest=bench.case_manifest(target)
    metadata=dict(case=str(target),cells=512000,faces=2064032,acceptedSteps=20,
        measuredIntervals=19,deltaT=.01,finalTime=.2,input=manifest,
        unchangedMeshAndNumericalControls=True,allVariantsAndRanksUseThisFixedWork=True)
    (ROOT/'voronoi20-input.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata),flush=True)


if __name__=='__main__':
    main()
