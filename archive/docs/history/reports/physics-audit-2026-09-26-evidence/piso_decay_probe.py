"""Reuse existing case generation only; compare production PISO with analytic decay."""
from pathlib import Path
import sys, json
ROOT=Path(__file__).resolve().parents[3]
source=ROOT/'tests/rans_validation_test.py'
namespace={'__file__':str(source),'__name__':'audit_fixture'}
sys.argv=[str(source),'--solver',str(ROOT/'build-petsc/babelsim-solve'),'--equations',str(ROOT/'build-petsc/rans_equations_test')]
exec(compile(source.read_text().split('summary = {')[0],str(source),'exec'),namespace)
results=[]
for model in ['kOmega','kEpsilon','SA']:
 for alpha in [0.7,1.0]:
  for dt in ([0.01,0.005,0.0025,0.0005,0.0001] if model=="kOmega" and alpha==0.7 else [0.01,0.005,0.0025]):
   fixture=namespace['case'](f'piso-{model}-{alpha}-{dt}',model,'euler',dt,0.1)
   p=fixture/'case.bs';p.write_text(p.read_text().replace('solver transientSimple','solver piso'))
   p=fixture/'numerics/solution.bs';p.write_text(p.read_text().replace('maxIterations 3000','maxIterations 1').replace('turbulenceRelaxation 0.7',f'turbulenceRelaxation {alpha}'))
   namespace['run'](fixture,1,'audit')
   expected=namespace['exact'](model,0.1)
   fields={name:list(namespace['values'](fixture,'audit',name).values()) for name in expected}
   actual={name:values[0] for name,values in fields.items()}
   ranges={name:[min(values),max(values)] for name,values in fields.items()}
   record=dict(model=model,alpha=alpha,dt=dt,actual=actual,ranges=ranges,expected=expected,max_error=max(abs(actual[k]-v) for k,v in expected.items()))
   results.append(record);print(json.dumps(record),flush=True)
print('fixture_directory='+str(namespace['BASE']))
(ROOT/'docs/reports/physics-audit-2026-09-26-evidence/piso_decay.json').write_text(json.dumps(results,indent=2)+'\n')
