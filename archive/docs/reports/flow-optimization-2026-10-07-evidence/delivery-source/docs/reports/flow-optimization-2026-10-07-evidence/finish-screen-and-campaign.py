#!/usr/bin/env python3
"""Wait for the current diagnostic driver, audit it, then run the campaign."""
from pathlib import Path
import fcntl
import json
import subprocess
import sys

ROOT=Path(__file__).resolve().parent


def run(script):
    print('START',script,flush=True)
    with (ROOT/(Path(script).stem+'-launcher.log')).open('a') as log:
        subprocess.run([sys.executable,'-u',str(ROOT/script)],
                       cwd=ROOT.parents[2],stdout=log,stderr=subprocess.STDOUT,check=True)
    print('FINISHED',script,flush=True)


def main():
    print('WAITING FOR THE SINGLE CURRENT DIAGNOSTIC TASK',flush=True)
    with open('/tmp/babelsim-solver-benchmark.lock','w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
    run('analyze-cg-single-reduction-screen.py')
    screen=json.loads((ROOT/'cg-single-reduction-screen-summary.json').read_text())
    # If the observed result changes, stop for source integration and regression;
    # never silently mix an overlay with the frozen production scaling binary.
    assert not screen['promising'], 'promising alternative requires production integration and regression before the final curve'
    slowdown=1/screen['pairedGeometricMeanSpeedup']-1
    reason=(f"保留已通过完整回归的 v12；单归约在三组工程规模配对中平均慢 {100*slowdown:.2f}%，没有支持接入的速度收益"
            if slowdown>=0 else
            f"保留已通过完整回归的 v12；单归约配对加速几何平均为 {screen['pairedGeometricMeanSpeedup']:.4f} 倍，但配对结果未满足稳定收益门槛")
    decision=dict(adopted=False,productionBinaryUnchanged=True,
        reason=reason,
        fieldConsistencyPassed=screen['allFieldsPassed'],
        pairedGeometricMeanSpeedup=screen['pairedGeometricMeanSpeedup'],
        scope='this hardware and fixed-work engineering screen; not a claim that single-reduction CG is universally slower',
        productionBinarySha256='2114c30d7632f2fa500822cfafe2a134c69a5e280089bc312243ff3be2e7c53c')
    (ROOT/'cg-single-reduction-decision.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(decision,ensure_ascii=False),flush=True)
    run('finish-extended-studies.py')
    run('audit-delivery.py')
    print('NUMERICAL CAMPAIGN AND DELIVERY AUDIT COMPLETED',flush=True)


if __name__=='__main__':
    main()
