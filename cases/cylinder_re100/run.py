#!/usr/bin/env python3
"""Launch a detached, two-core PISO run with range guards and final postprocessing.
No scheduled service is installed. A local worker owns the mpirun process group.
"""
from pathlib import Path
import argparse,datetime,json,os,signal,subprocess,sys,time
import numpy as np
CASE=Path(__file__).resolve().parent
REPO=CASE.parents[1]
RUN=CASE/'results/run'
STATUS=CASE/'results/run_status.json'

def write_status(state):
    tmp=STATUS.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(STATUS)

def complete_snapshot(path,n):
    arrays=[]
    for name in ['U','p']:
        files=sorted(path.glob(f'rank-*/{name}.csv'))
        if len(files)!=2:return None
        try:a=np.concatenate([np.loadtxt(p,delimiter=',',skiprows=1,ndmin=2) for p in files])
        except (ValueError,OSError):return None
        if len(a)!=n or len(np.unique(a[:,0]))!=n:return None
        arrays.append(a[:,4:])
    u,p=arrays
    return {'finite':bool(np.isfinite(u).all() and np.isfinite(p).all()),
        'max_speed':float(np.linalg.norm(u,axis=1).max()),'max_abs_p':float(abs(p).max())}

def worker():
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    n=json.loads((CASE/'mesh/mesh_quality.json').read_text())['cells']
    cmd=['mpirun','-np','2','--bind-to','core','--map-by','core',str(REPO/'build/babelsim-solve'),'-case',str(CASE),'-time','run']
    started=time.monotonic();state={'status':'running','started_at':datetime.datetime.now().astimezone().isoformat(),
        'worker_pid':os.getpid(),'mpi_ranks':2,'threads_per_rank':1,'command':cmd,
        'control':(CASE/'control.bs').read_text(),'last_complete_time':None,
        'guard_note':'Range checks detect gross blowup; they are not physical validation.'}
    with (CASE/'results/solver.log').open('w') as log:
        process=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        state['mpirun_pid']=process.pid;write_status(state)
        stopped=False
        while process.poll() is None:
            snapshots=[]
            if RUN.exists():
                for p in RUN.iterdir():
                    if not p.is_dir():continue
                    try:t=float(p.name)
                    except ValueError:continue
                    if t>(state['last_complete_time'] or -1):snapshots.append((t,p))
            for t,p in sorted(snapshots):
                metrics=complete_snapshot(p,n)
                if metrics is None:continue
                state.update(last_complete_time=t,latest=metrics,elapsed_seconds=time.monotonic()-started)
                if not metrics['finite'] or metrics['max_speed']>10 or metrics['max_abs_p']>100:
                    state['status']='stopped_by_range_guard';stopped=True
                    os.killpg(process.pid,signal.SIGINT);write_status(state);break
            write_status(state)
            if stopped:break
            time.sleep(10)
        code=process.wait();state['solver_exit_code']=code;state['elapsed_seconds']=time.monotonic()-started
    if not stopped:state['status']='solver_completed' if code==0 else 'solver_failed'
    write_status(state)
    if code==0:
        with (CASE/'results/postprocess.log').open('w') as out:
            post=subprocess.run([sys.executable,str(CASE/'postprocess.py'),'--dt','.001','--gif'],stdout=out,stderr=subprocess.STDOUT,env=env)
        state['postprocess_exit_code']=post.returncode
        state['status']='completed' if post.returncode==0 else 'postprocess_failed';write_status(state)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--worker',action='store_true');args=parser.parse_args()
    (CASE/'results').mkdir(exist_ok=True)
    if args.worker:worker()
    else:
        if STATUS.exists():
            previous=json.loads(STATUS.read_text())
            if previous.get('status')=='running':raise SystemExit('A run is already marked running; inspect run_status.json before launching again.')
        if RUN.exists():raise SystemExit('Refusing to overwrite existing production results/run')
        with (CASE/'results/worker.log').open('w') as out:
            p=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--worker'],stdout=out,stderr=subprocess.STDOUT,start_new_session=True,
                env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
        print(f'Started worker PID {p.pid}; status: {STATUS}')
