"""PISO must distinguish a conservative flux from a settled momentum correction."""
import argparse
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--solver', type=Path, required=True)
solver = parser.parse_args().solver.resolve()

def invoke(case, ranks):
    result = subprocess.run(['mpirun', '-np', str(ranks), str(solver), '-case', str(case)],
                            env={**os.environ, 'OMP_NUM_THREADS': '1'},
                            capture_output=True, text=True, timeout=120)
    rows = [dict(re.findall(r'(\w+)=([^ ]+)', line)) for line in result.stdout.splitlines()
            if line.startswith('PISO 1 ')]
    assert rows, result.stdout + result.stderr
    return result, rows[-1]

with tempfile.TemporaryDirectory(prefix='babelsim-piso-coupling-') as temporary:
    root = Path(temporary)
    baseline = root / 'baseline'
    shutil.copytree(ROOT/'cases/cavity', baseline,
                    ignore=shutil.ignore_patterns('results', 'post'))
    path = baseline/'case.bs'
    path.write_text(path.read_text().replace('solver simple', 'solver piso') + '\nghostLayers 3\n')
    (baseline/'control.bs').write_text('startTime 0\nendTime 0.01\ndeltaT 0.01\n')
    (baseline/'output.bs').write_text('directory results\ntimeName final\nwriteInterval 1\n')
    path = baseline/'numerics/methods.bs'
    path.write_text(path.read_text().replace('time steady', 'time bdf2'))
    path = baseline/'numerics/solution.bs'
    settings = '\n'.join(line for line in path.read_text().splitlines()
                         if not line.startswith(('maxIterations ', 'velocityRelaxation ', 'pressureRelaxation ')))
    settings += '\nmaxIterations 1\nvelocityRelaxation 1\nnCorrectors 1\n'
    path.write_text(settings)
    process, row = invoke(baseline, 1)
    assert process.returncode == 0, process.stdout + process.stderr
    tolerance = float(row['rCoupling']) * 0.5
    assert tolerance > 1e-10, row
    reference = None
    for ranks in (1, 2, 4):
        for adaptive in (False, True):
            case = root/f'np{ranks}-adaptive{adaptive}'
            shutil.copytree(baseline, case, ignore=shutil.ignore_patterns('results'))
            (case/'numerics/solution.bs').write_text(settings +
                f'couplingTolerance {tolerance:.17g}\nmaxCorrectors {30 if adaptive else 1}\n')
            process, row = invoke(case, ranks)
            if not adaptive:
                assert process.returncode == 2, process.stdout + process.stderr
                assert row['accepted'] == 'false' and row['converged'] == 'false', row
                assert float(row['mass']) < 1e-8 and row['linear'] == 'ok', row
                assert not (case/'results').exists()
            else:
                assert process.returncode == 0, process.stdout + process.stderr
                assert row['accepted'] == 'true' and int(row['correctors']) > 1, row
                assert float(row['rCoupling']) <= tolerance and float(row['CoMax']) > 0, row
                data = {}
                for path in (case/'results/final').glob('rank-*/U.csv'):
                    for value in csv.DictReader(path.open()):
                        data[int(value['global_id'])] = tuple(float(value[f'value{i}']) for i in range(3))
                assert data
                if reference is None:
                    reference = data
                else:
                    assert data.keys() == reference.keys()
                    assert max(abs(a-b) for i in data for a,b in zip(data[i], reference[i])) < 1e-7
    print('piso_coupling_test: conservative-but-unsettled rejection, adaptive correction, 1/2/4 ranks passed')
