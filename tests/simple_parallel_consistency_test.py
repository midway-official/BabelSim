#!/usr/bin/env python3
"""Full SIMPLE 1/2/4-rank comparisons, including every transient output step."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from compare_parallel_results import read_result


def generated_mesh(path, shape, warp):
    nx, ny, nz = shape
    def vid(i,j,k): return i+(nx+1)*(j+(ny+1)*k)
    points=[]
    for k in range(nz+1):
        for j in range(ny+1):
            for i in range(nx+1):
                x,y,z=i/nx,j/ny,k/nz
                if warp: x += .18*y + .08*y*(1-y)*math.sin(math.pi*z)
                points.append((x,y,z))
    cells=[]; faces=[[] for _ in range(6)]
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c=[vid(i,j,k),vid(i+1,j,k),vid(i+1,j+1,k),vid(i,j+1,k),
                   vid(i,j,k+1),vid(i+1,j,k+1),vid(i+1,j+1,k+1),vid(i,j+1,k+1)]
                cells.append(c)
                for side,yes,ids in [(0,i==0,[0,4,7,3]),(1,i==nx-1,[1,2,6,5]),
                    (2,j==0,[0,1,5,4]),(3,j==ny-1,[3,7,6,2]),
                    (4,k==0,[0,3,2,1]),(5,k==nz-1,[4,5,6,7])]:
                    if yes: faces[side].append([c[t] for t in ids])
    with path.open('w') as f:
        f.write(f'BABELSIM_MESH 2\nvertices {len(points)}\n')
        for p in points: f.write(' '.join(format(x,'.17g') for x in p)+'\n')
        f.write(f'cells {len(cells)}\n')
        for c in cells: f.write(' '.join(map(str,c))+'\n')
        f.write('patches 6\n')
        for name,patch in zip(['cavity_left','cavity_right','cavity_bottom','lid','front','back'],faces):
            f.write(f'patch {name} wall {len(patch)}\n')
            for face in patch: f.write(' '.join(map(str,face))+'\n')
        f.write('end\n')


def prepare(base, name, mode, ghost=3):
    source=ROOT/'cases'/('poiseuille' if name=='channel' else 'cavity')
    target=base/(name+'-'+mode)
    shutil.copytree(source,target,ignore=shutil.ignore_patterns('results','post','validation','benchmark*'))
    entry=target/'case.bs'
    entry.write_text(entry.read_text().replace('solver simple','solver '+('simple' if mode=='steady' else 'transientSimple'))+f'\nghostLayers {ghost}\n')
    if name.startswith(('cube','warped')):
        warp=name.startswith('warped')
        generated_mesh(target/'mesh/cavity.mesh',(10,8,6) if warp else (8,8,6),warp)
        for field in ('U','p'):
            p=target/f'fields/initial/{field}.field'
            s=p.read_text().replace('front { type symmetry }','front { type '+('fixedValue value (0 0 0)' if field=='U' else 'zeroGradient')+' }')
            s=s.replace('back { type symmetry }','back { type '+('fixedValue value (0 0 0)' if field=='U' else 'zeroGradient')+' }')
            p.write_text(s)
        (target/'numerics/methods.bs').write_text('interpolation '+('corrected' if warp else 'linear')+'\ngradient '+('greenGauss' if warp else 'leastSquares')+'\nconvection '+('linearUpwind' if warp else 'upwind')+'\ndiffusion '+('corrected' if warp else 'orthogonal')+'\ntime steady\n')
        (target/'numerics/solution.bs').write_text('maxIterations 8000\nnonOrthogonalCorrections 2\nvelocityRelaxation 0.5\npressureRelaxation 0.3\ncontinuityTolerance 1e-9\nvelocityTolerance 1e-8\npressureCorrectionTolerance 1e-8\nmomentumTolerance 1e-8\nvectorSolver bicgstab ilut 1e-15 1e-11 2000\nscalarSolver cg incompleteCholesky 1e-15 1e-11 2000\n')
    if mode!='steady':
        methods=target/'numerics/methods.bs'
        methods.write_text(methods.read_text().replace('time steady','time '+mode))
        (target/'control.bs').write_text('startTime 0\nendTime 0.05\ndeltaT 0.01\n')
        # Stricter inner/outer accuracy for comparing transient pressure at every step.
        solution=target/'numerics/solution.bs'
        s=solution.read_text()
        for key,value in [('maxIterations','8000'),('velocityTolerance','1e-8'),('pressureCorrectionTolerance','1e-8'),('momentumTolerance','1e-8')]:
            s=re.sub(r'^'+key+r' .+$',key+' '+value,s,flags=re.M) if re.search(r'^'+key+r' ',s,re.M) else s+'\n'+key+' '+value+'\n'
        s=re.sub(r'^vectorSolver .+$','vectorSolver bicgstab ilut 1e-15 1e-11 2000',s,flags=re.M)
        s=re.sub(r'^scalarSolver .+$','scalarSolver cg incompleteCholesky 1e-15 1e-11 2000',s,flags=re.M)
        solution.write_text(s)
    (target/'output.bs').write_text('directory results\ntimeName final\nwriteInterval 1\n')
    mesh=next((target/'mesh').glob('*.mesh'))
    count=int(re.search(r'^cells (\d+)$',mesh.read_text(),re.M)[1])
    return target,count


def snapshot(path,ranks,count,expected_time):
    directories=sorted(path.glob('rank-*'))
    assert len(directories)==ranks,(path,'rank count')
    geometry={}
    for rank,directory in enumerate(directories):
        lines=[s.split() for s in (directory/'metadata.bs').read_text().splitlines()]
        metadata={s[0]:s[1:] for s in lines if s[0]!='field'}
        assert metadata['format']==['babelsim_result','2']
        assert int(metadata['rank'][0])==rank and int(metadata['ranks'][0])==ranks
        assert int(metadata['global_cell_count'][0])==count
        assert abs(float(metadata['time'][0])-expected_time)<1e-12
        rows=(directory/'mesh.geometry').read_text().splitlines()
        assert len(rows)==int(metadata['owned_cells'][0])
        for row in rows:
            t=row.split(','); gid=int(t[0]); assert gid not in geometry
            geometry[gid]=tuple(map(float,t[1:]))
    fields=read_result(path)
    assert set(fields)=={'U','p'},(path,fields.keys())
    assert set(geometry)==set(range(count))
    for field in fields.values(): assert set(field)==set(range(count))
    return fields,geometry


def compare(reference,candidate,atol,rtol):
    a,ag=reference; b,bg=candidate
    assert ag==bg,'geometry differs across partitions'
    out={}
    for name in a:
        errors=[]; values=[]; worst=0.0
        for gid in a[name]:
            assert len(a[name][gid])==len(b[name][gid])
            for x,y in zip(a[name][gid],b[name][gid]):
                error=abs(x-y); errors.append(error); values.append(x*x)
                worst=max(worst,error/(atol+rtol*max(abs(x),abs(y))))
        out[name]={'max_abs':max(errors),'rms':math.sqrt(sum(e*e for e in errors)/len(errors)),
                   'relative_l2':math.sqrt(sum(e*e for e in errors)/sum(values)) if sum(values)>0 else max(errors),
                   'tolerance_ratio':worst}
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--cases',nargs='+',default=['cavity','channel','cube','warped','warped4'])
    parser.add_argument('--modes',nargs='+',default=['steady','euler','bdf2'])
    args=parser.parse_args()
    base=args.output or Path(tempfile.mkdtemp(prefix='babelsim-simple-mpi-'))
    base.mkdir(parents=True,exist_ok=True)
    summary={'directory':str(base),'baseline':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
             'binary_sha256':hashlib.sha256((ROOT/'build/babelsim-solve').read_bytes()).hexdigest(),
             'atol':5e-6,'rtol':5e-6,'runs':[],'comparisons':[],'passed':False}
    print('Evidence:',base,flush=True)
    for name in args.cases:
        for mode in args.modes:
            case,count=prepare(base,name,mode,4 if name=='warped4' else 3)
            reference={}
            for ranks in (1,2,4):
                label=f'np{ranks}'
                command=['mpirun','-np',str(ranks),str(ROOT/'build/babelsim-solve'),'-case',str(case),'-time',label]
                begin=time.monotonic()
                result=subprocess.run(command,cwd=ROOT,env=dict(os.environ,TMPDIR='/tmp'),text=True,
                                      stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=600)
                (base/f'{name}-{mode}-{label}.log').write_text(result.stdout)
                converged=[line for line in result.stdout.splitlines() if 'converged=true' in line]
                run={'case':name,'mode':mode,'ranks':ranks,'cells':count,'ghost_layers':4 if name=='warped4' else 3,
                     'exit':result.returncode,'seconds':time.monotonic()-begin,'converged_steps':len(converged),'terminal':converged}
                summary['runs'].append(run)
                (base/'summary.json').write_text(json.dumps(summary,indent=2))
                assert result.returncode==0,(name,mode,ranks,result.stdout[-2000:])
                assert len(converged)==(1 if mode=='steady' else 5),(name,mode,ranks,'missing converged time step')
                times=[('final',0.0)] if mode=='steady' else [(format(i*.01,'.2f'),i*.01) for i in range(1,6)]
                root=case/'results'/label
                if mode!='steady':
                    numeric=sorted(float(p.name) for p in root.iterdir() if p.is_dir() and re.fullmatch(r'\d+(\.\d+)?',p.name))
                    assert numeric==[t for _,t in times],(name,mode,ranks,numeric)
                for timestamp,physical_time in times:
                    current=snapshot(root if timestamp=='final' else root/timestamp,ranks,count,physical_time)
                    if ranks==1: reference[timestamp]=current
                    else:
                        diff=compare(reference[timestamp],current,summary['atol'],summary['rtol'])
                        summary['comparisons'].append({'case':name,'mode':mode,'ranks':ranks,'time':physical_time,'fields':diff})
                        (base/'summary.json').write_text(json.dumps(summary,indent=2))
                        assert all(v['tolerance_ratio']<=1 for v in diff.values()),(name,mode,ranks,timestamp,diff)
                if mode!='steady':
                    final=snapshot(root,ranks,count,.05)
                    last=snapshot(root/'0.05',ranks,count,.05)
                    assert final==last,'final alias differs from final physical state'
                print(f'{name}/{mode}/np{ranks}: {count} cells, {len(converged)} converged step(s), {run["seconds"]:.2f}s',flush=True)
    summary['passed']=True
    (base/'summary.json').write_text(json.dumps(summary,indent=2))
    print('simple_parallel_consistency_test: all states, fields, geometry, and physical times passed',flush=True)

if __name__=='__main__': main()
