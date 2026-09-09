from pathlib import Path
import shutil, subprocess, csv, json
root=Path('/home/midway/BabelSim')
base=Path('/tmp/babelsim-review-20260909')
patches=['hot','cold','lower','upper','front','back']
def field(name,kind,value,bc='zeroGradient'):
 return f'field {name}\n{{\ntype {kind}\nlocation cell\ninternal uniform ({value})\nboundary\n{{\n'+''.join(f'{p} {{ type {bc}'+(f' value ({value})' if bc=='fixedValue' else '')+' }\n' for p in patches)+'}\n}\n'
def case(label,model,time,clip=False):
 d=base/label;shutil.copytree(root/'cases/heat',d,ignore=shutil.ignore_patterns('results','post'),dirs_exist_ok=True)
 (d/'case.bs').write_text((d/'case.bs').read_text().replace('solver heat','solver '+('simple' if time=='steady' else 'transientSimple')))
 (d/'numerics/methods.bs').write_text(f'interpolation linear\ngradient leastSquares\nconvection upwind\ndiffusion orthogonal\ntime {time}\n')
 (d/'control.bs').write_text('startTime 0\nendTime '+('1' if clip else '0.001')+'\ndeltaT '+('1' if clip else '0.001')+'\n')
 (d/'numerics/solution.bs').write_text('scalarSolver bicgstab ilut 1e-14 1e-10 1000\nvectorSolver bicgstab ilut 1e-14 1e-10 1000\nmaxIterations 2000\nvelocityTolerance 1e-8\ncontinuityTolerance 1e-8\npressureCorrectionTolerance 1e-6\n')
 phy='density 1\ndynamicViscosity 0.01\nturbulenceModel '+model+'\n'
 if model!='none':phy+='turbulenceRelaxation 0.5\nturbulenceTolerance 1e-8\n'
 if clip:phy+='kMin 1\nepsilonMin 1\n'
 (d/'physics/thermal.bs').write_text(phy)
 for n,k,v in [('U','vector','0 0 0'),('p','scalar','1'),('k','scalar','1'),('epsilon','scalar','1'),('omega','scalar','1'),('nuTilda','scalar','0.03'),('wallDistance','scalar','0.5')]:
  (d/'fields/initial'/f'{n}.field').write_text(field(n,k,v,'fixedValue' if n=='U' or (time=='steady' and n not in ['p','wallDistance']) else 'zeroGradient'))
 return d
rows=[]
for time in ['steady','euler']:
 for model in ['none','SA','kOmega','kEpsilon']:
  d=case(f'smoke-{time}-{model}',model,time)
  for n in [1,2]:
   p=subprocess.run(['mpirun','-np',str(n),str(root/'build/babelsim-solve'),'-case',str(d),'-time',f'np{n}'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=90,env={**__import__('os').environ,'TMPDIR':'/tmp'})
   (base/f'{d.name}-np{n}.log').write_text(p.stdout)
   row={'case':d.name,'ranks':n,'exit':p.returncode}
   for f in ['k','omega','epsilon','nuTilda','mut']:
    path=d/'results'/f'np{n}'/'rank-0000'/f'{f}.csv'
    if path.exists():
     data=list(csv.DictReader(path.open()));row[f]=data[0]
   rows.append(row);print(json.dumps(row),flush=True)
d=case('clipping-false-convergence','kEpsilon','euler',True)
for n in [1,2]:
 p=subprocess.run(['mpirun','-np',str(n),str(root/'build/babelsim-solve'),'-case',str(d),'-time',f'np{n}'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=90,env={**__import__('os').environ,'TMPDIR':'/tmp'})
 (base/f'clipping-np{n}.log').write_text(p.stdout)
 print('CLIPPING',n,p.returncode,p.stdout,flush=True)
(base/'rans-summary.json').write_text(json.dumps(rows,indent=2))
