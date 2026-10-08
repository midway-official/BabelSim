from pathlib import Path
import sys,json,subprocess
ROOT=Path(__file__).resolve().parents[3]
source=ROOT/'tests/rans_validation_test.py'
n={'__file__':str(source),'__name__':'audit_fixture'}
sys.argv=[str(source),'--solver',str(ROOT/'build-petsc/babelsim-solve'),'--equations',str(ROOT/'build-petsc/rans_equations_test')]
exec(compile(source.read_text().split('summary = {')[0],str(source),'exec'),n)
f=n['case']('piso-failed-turbulence','kOmega','euler',0.01,0.01,fixed=True)
p=f/'case.bs';p.write_text(p.read_text().replace('solver transientSimple','solver piso'))
p=f/'numerics/solution.bs';p.write_text(p.read_text().replace('maxIterations 3000','maxIterations 1').replace('equation.kTransport.maxIterations 2000','equation.kTransport.maxIterations 1').replace('equation.kTransport.pcType bjacobi','equation.kTransport.pcType jacobi'))
result=subprocess.run([str(ROOT/'build-petsc/babelsim-solve'),'-case',str(f),'-time','audit','-performance',str(ROOT/'docs/reports/physics-audit-2026-09-26-evidence/piso-failure-performance')],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
print('exit_code='+str(result.returncode));print('fixture='+str(f));print(result.stdout)
