#!/usr/bin/env python3
"""Run a bounded full-mesh check without modifying production dictionaries."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--steps',type=int,default=5)
p.add_argument('--ranks',type=int,default=4)
p.add_argument('--dt',type=float,default=0.01)
p.add_argument('--solver',type=Path,default=REPO/'build-petsc/babelsim-solve')
p.add_argument('--runs-directory',type=Path,default=ROOT/'results')
p.add_argument('--performance',action='store_true',help='write per-rank PETSc timing JSON')
p.add_argument('--pressure-ksp',choices=('cg','bcgs','gmres','fgmres'),help='override pressure-correction KSP in this isolated run')
p.add_argument('--pressure-pc',choices=('bjacobi','hypre','gamg'),help='override pressure-correction PC in this isolated run')
a=p.parse_args()
if a.steps<1 or a.ranks<1: p.error('steps and ranks must be positive')
if a.dt<0.005: p.error('deltaT must be >= 0.005 (user-required lower bound)')
a.runs_directory.mkdir(parents=True,exist_ok=True)
if not a.solver.is_file(): p.error(f'solver executable does not exist: {a.solver}')
run=a.runs_directory.resolve()/('smoke-'+time.strftime('%Y%m%d-%H%M%S'))
run.mkdir(parents=True)
for name in ('case.bs','output.bs'): shutil.copy2(ROOT/name,run/name)
for name in ('fields','numerics','physics'): shutil.copytree(ROOT/name,run/name)
(run/'mesh').symlink_to(ROOT/'mesh',target_is_directory=True)
if a.pressure_ksp or a.pressure_pc:
    solution=run/'numerics/solution.bs'
    text=solution.read_text()
    for key,value in (('kspType',a.pressure_ksp),('pcType',a.pressure_pc)):
        if value is None: continue
        text,count=re.subn(
            rf'^(equation\.pressureCorrection\.{key}\s+).+$',
            rf'\g<1>{value}',text,flags=re.M)
        if count != 1: p.error(f'expected one pressureCorrection.{key} setting, found {count}')
    solution.write_text(text)
(run/'control.bs').write_text(f'startTime 0\nendTime {a.steps*a.dt:.12g}\ndeltaT {a.dt:.12g}\n')
(run/'output.bs').write_text('directory results\ntimeName final\nwriteInterval 1\n')
cmd=['mpirun','-np',str(a.ranks),str(a.solver.resolve()),'-case',str(run)]
if a.performance: cmd.extend(['-performance',str(run/'performance')])
env=os.environ.copy(); env['OMP_NUM_THREADS']='1'
mesh_hash=hashlib.sha256((ROOT/'mesh/naca0012.mesh').read_bytes()).hexdigest()
start=time.monotonic()
with (run/'solver.log').open('w') as log:
    result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
report={'mesh_sha256':mesh_hash,'command':cmd,'returncode':result.returncode,'steps_requested':a.steps,'deltaT':a.dt,'pressure_ksp_override':a.pressure_ksp,'pressure_pc_override':a.pressure_pc,'elapsed_seconds':time.monotonic()-start,'run_directory':str(run)}
(run/'run_summary.json').write_text(json.dumps(report,indent=2)+'\n')
from summarize_run import summarize
report=summarize(run)
print(json.dumps({k:v for k,v in report.items() if k!='saved_steps'},indent=2))
raise SystemExit(0 if report['smoke_passed'] else 1)
