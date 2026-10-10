from pathlib import Path
import subprocess,os,time,json,re,signal
root=Path(__file__).resolve().parent
solver=Path('/home/midway/BabelSim/build-petsc/babelsim-solve')
for name in ['upwind_probe','pressure_corrected','correctors8']:
 case=root/name; start=time.time();reason=None
 with (case/'solver.log').open('w') as log:
  p=subprocess.Popen(['mpirun','-np','4',str(solver),'-case',str(case)],stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'OMP_NUM_THREADS':'1'},start_new_session=True)
  while p.poll() is None:
   time.sleep(2)
   text=(case/'solver.log').read_text();v=re.findall(r'Umax=(\S+)',text)
   if v and (not float(v[-1])<20):
    reason='stopped: Umax exceeds 20 m/s diagnosis threshold';os.killpg(p.pid,signal.SIGTERM);break
  code=p.wait()
 (case/'status.json').write_text(json.dumps({'exit_code':code,'wall_seconds':time.time()-start,'reason':reason}))
 print(name,code,reason,flush=True)
