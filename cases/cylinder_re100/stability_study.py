#!/usr/bin/env python3
"""Sequential two-core, one-factor-at-a-time startup stability experiments."""
from pathlib import Path
import json,os,re,shutil,subprocess,time
import numpy as np
from postprocess import Geometry,field
CASE=Path(__file__).resolve().parent
REPO=CASE.parents[1]
EXPERIMENTS=[
 ('bdf2_dt005',.005,{},{}),
 ('bdf2_dt005_one_corrector',.005,{}, {'nCorrectors':'1'}),
 ('euler_dt005',.005,{'time':'euler'},{}),
 ('bdf2_dt005_limited',.005,{'equation.momentum.diffusion':'limitedCorrected','equation.pressureCorrection.diffusion':'limitedCorrected'},{}),
 ('bdf2_dt005_upwind',.005,{'equation.momentum.convection':'upwind'},{}),
 ('bdf2_dt002',.002,{},{}),
]

def replace(s,changes):
 for k,v in changes.items():s=re.sub(r'^'+re.escape(k)+r'\s+.*$',k+' '+v,s,flags=re.M)
 return s

def main(experiments=EXPERIMENTS, report_name="stability_study"):
 root=CASE/'results/stability';root.mkdir(parents=True,exist_ok=True);geo=Geometry();report=[]
 env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
 for name,dt,method,solution in experiments:
  path=root/name;path.mkdir(exist_ok=True)
  for item in ['mesh','physics','fields']:
   if not (path/item).exists():(path/item).symlink_to(CASE/item,target_is_directory=True)
  shutil.copyfile(CASE/'case.bs',path/'case.bs');(path/'numerics').mkdir(exist_ok=True)
  for item,changes in [('methods',method),('solution',solution)]:
   (path/f'numerics/{item}.bs').write_text(replace((CASE/f'numerics/{item}.bs').read_text(),changes))
  (path/'control.bs').write_text(f'startTime 0\nendTime 0.4\ndeltaT {dt}\n')
  (path/'output.bs').write_text('directory results\ntimeName run\nwriteInterval 5\nwriteFields U p\n')
  cmd=['mpirun','-np','2','--bind-to','core','--map-by','core',str(REPO/'build-petsc/babelsim-solve'),'-case',str(path),'-time','run']
  started=time.monotonic()
  with (path/'solver.log').open('w') as out:p=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT,env=env)
  lines=(path/'solver.log').read_text().splitlines();steps=[x for x in lines if x.startswith('PISO ')]
  samples=[]
  for d in (path/'results/run').glob('*'):
   try:t=float(d.name)
   except ValueError:continue
   try:u=field(d,'U',geo.n);pressure=field(d,'p',geo.n)
   except (ValueError,OSError):continue
   _,co,_=geo.diagnostics(u,pressure,dt)
   imax=int(np.argmax(np.linalg.norm(u,axis=1)));ip=int(np.argmax(abs(pressure[:,0])))
   samples.append({'t':t,'max_speed':float(np.linalg.norm(u,axis=1).max()),'max_abs_p':float(abs(pressure).max()),
     'max_Co_interior':float(co.max()),'speed_max_xy':geo.c[imax].tolist(),'pressure_max_xy':geo.c[ip].tolist()})
  samples.sort(key=lambda x:x['t'])
  entry={'name':name,'dt':dt,'methods_changes':method,'solution_changes':solution,'exit_code':p.returncode,
    'elapsed_seconds':time.monotonic()-started,'reported_steps':len(steps),'last_record':steps[-1] if steps else None,'samples':samples}
  report.append(entry);(CASE/f'validation/{report_name}.json').write_text(json.dumps(report,indent=2)+'\n')
  print(json.dumps({k:v for k,v in entry.items() if k!='samples'}),flush=True)
  print('last sample: '+str(samples[-1] if samples else None),flush=True)
if __name__=='__main__':main()
