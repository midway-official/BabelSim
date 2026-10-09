#!/usr/bin/env python3
"""Run an immutable case snapshot; monitor progress and postprocess completion."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from zoneinfo import ZoneInfo

import numpy as np

CASE = Path(__file__).resolve().parent
REPO = CASE.parents[1]
ENV = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
       'MKL_NUM_THREADS': '1', 'TMPDIR': '/tmp'}


def now():
    return datetime.datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()


def save_status(path, state):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(state, indent=2)+'\n')
    temporary.replace(path)


def snapshot_metrics(path, n, ranks, viscosity):
    values = {}
    for name in ['U', 'p', 'k', 'omega', 'mut']:
        files = sorted(path.glob(f'rank-*/{name}.csv'))
        if len(files) != ranks:
            return None
        try:
            array = np.concatenate([np.loadtxt(p, delimiter=',', skiprows=1, ndmin=2) for p in files])
        except (ValueError, OSError):
            return None
        if len(array) != n or not np.array_equal(np.sort(array[:, 0]), np.arange(n)):
            return None
        values[name] = array[:, 4:]
    return {'finite': bool(all(np.isfinite(a).all() for a in values.values())),
            'max_speed': float(np.linalg.norm(values['U'], axis=1).max()),
            'max_abs_p': float(abs(values['p']).max()),
            'k_range': [float(values['k'].min()), float(values['k'].max())],
            'omega_range': [float(values['omega'].min()), float(values['omega'].max())],
            'mut_range': [float(values['mut'].min()), float(values['mut'].max())],
            'max_mut_over_mu': float(values['mut'].max()/viscosity)}


def worker(manifest_path):
    manifest = json.loads(manifest_path.read_text())
    label = manifest['label']
    directory = CASE/'results'
    status_path = directory/f'{label}_status.json'
    log_path = directory/('solver.log' if label == 'run' else f'{label}_solver.log')
    results = directory/label
    parameters = manifest['parameters']
    n = manifest['cells']
    started = time.monotonic()
    binary = Path(manifest.get('binary_path', str(REPO/'build-petsc/babelsim-solve')))
    command = ['mpirun', '-np', str(manifest['ranks']), '--bind-to', 'core', '--map-by', 'core',
               str(binary), '-case', manifest['input_case'], '-time', label,
               '-performance', str(directory/(label+'_performance'))]
    state = {**manifest, 'status': 'running', 'started_at': now(), 'worker_pid': os.getpid(),
             'command': command, 'solver_log': str(log_path), 'last_complete_time': None,
             'reported_steps': 0, 'reported_outer_iterations': 0, 'reported_outer_records': 0,
             'last_solver_time': 0., 'last_iteration_time': 0., 'max_mass': 0., 'max_CoMax': 0.,
             'max_rCoupling': 0., 'max_rTurb': 0., 'all_steps_accepted': True,
             'all_linear_ok': True, 'all_turbulence_settled': True}
    try:
        if hashlib.sha256(binary.read_bytes()).hexdigest() != manifest['binary_sha256']:
            raise ValueError('Solver binary changed after the run manifest was written')
    except (OSError, ValueError) as error:
        state.update(status='worker_failed', error=str(error), updated_at=now())
        save_status(status_path, state)
        return 1
    stopped = False
    position = 0
    previous_time, previous_index = None, 0
    with log_path.open('w') as output:
        process = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT,
                                   env=ENV, start_new_session=True)
        state['mpirun_pid'] = process.pid
        save_status(status_path, state)
        try:
            while True:
                with log_path.open(errors='replace') as source:
                    source.seek(position)
                    while True:
                        line_position = source.tell()
                        line = source.readline()
                        if not line or not line.endswith('\n'):
                            position = line_position
                            break
                        if re.match(r'^PISO \d+ time=', line):
                            record = dict(re.findall(r'(\w+)=(\S+)', line))
                            index = int(line.split()[1])
                            current_time = float(record['time'])
                            state['reported_outer_records'] += 1
                            state['reported_outer_iterations'] += index - previous_index if current_time == previous_time else index
                            previous_time, previous_index = current_time, index
                            state['last_iteration_time'] = current_time
                            for key in ['mass', 'CoMax', 'rCoupling']:
                                state['max_'+key] = max(state['max_'+key], float(record[key]))
                            state['all_steps_accepted'] &= record.get('accepted') == 'true'
                            state['all_linear_ok'] &= record.get('linear') == 'ok'
                            # Intermediate Picard iterations are expected to be
                            # unsettled. Count a physical step only at convergence.
                            if record.get('converged') == 'true':
                                state['reported_steps'] += 1
                                state['last_solver_time'] = float(record['time'])
                                state['max_rTurb'] = max(state['max_rTurb'], float(record['rTurb']))
                                state['all_turbulence_settled'] &= record.get('turbulenceSettled') == 'true'
                state['elapsed_seconds'] = time.monotonic()-started
                state['updated_at'] = now()
                if state['last_solver_time'] > 0:
                    state['eta_seconds_estimate'] = state['elapsed_seconds']*max(
                        0., manifest['end_time']/state['last_solver_time']-1)
                if results.exists():
                    snapshots = []
                    for path in results.iterdir():
                        try:
                            t = float(path.name)
                        except ValueError:
                            continue
                        if path.is_dir() and t > (state['last_complete_time'] if state['last_complete_time'] is not None else -1):
                            snapshots.append((t, path))
                    for t, path in sorted(snapshots):
                        metrics = snapshot_metrics(path, n, manifest['ranks'], parameters['dynamicViscosity'])
                        if metrics is None:
                            continue
                        state.update(last_complete_time=t, latest_fields=metrics)
                        if (not metrics['finite'] or metrics['max_speed'] > 10 or metrics['max_abs_p'] > 100
                                or metrics['k_range'][0] < 0 or metrics['omega_range'][0] <= 0
                                or metrics['mut_range'][0] < 0):
                            state['status'] = 'stopped_by_field_guard'
                            stopped = True
                            if process.poll() is None:
                                os.killpg(process.pid, signal.SIGINT)
                            break
                save_status(status_path, state)
                if process.poll() is not None or stopped:
                    break
                time.sleep(5)
            code = process.wait(timeout=30)
        except BaseException as exc:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
            state.update(status='worker_failed', error=str(exc), updated_at=now())
            save_status(status_path, state)
            raise
    state.update(solver_exit_code=code, elapsed_seconds=time.monotonic()-started, updated_at=now())
    if not stopped:
        state['status'] = 'solver_completed' if code == 0 else 'solver_failed'
    if code == 0 and state['last_solver_time'] < manifest['end_time']-1e-7:
        state['status'] = 'solver_ended_before_target'
    save_status(status_path, state)
    if state['status'] == 'solver_completed':
        state['status'] = 'postprocessing'
        save_status(status_path, state)
        analysis = CASE/'validation/production' if label == 'run' else CASE/'validation'/label/'analysis'
        with (directory/f'{label}_postprocess.log').open('w') as output:
            post = subprocess.run([sys.executable, str(CASE/'postprocess.py'), '--results', str(results),
                                   '--out', str(analysis), '--log', str(log_path),
                                   '--case', manifest['input_case'], '--require-complete', '--gif'],
                                  stdout=output, stderr=subprocess.STDOUT, env=ENV)
        state.update(postprocess_exit_code=post.returncode, analysis_directory=str(analysis),
                     status='completed' if post.returncode == 0 else 'postprocess_failed', updated_at=now())
        save_status(status_path, state)
    return 0 if state['status'] == 'completed' else 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--label', default='run')
    parser.add_argument('--ranks', type=int, default=8)
    parser.add_argument('--binary', type=Path, default=REPO/'build-petsc/babelsim-solve',
                        help='solver binary to freeze in the run manifest')
    parser.add_argument('--end-time', type=float)
    parser.add_argument('--dt', type=float)
    parser.add_argument('--write-interval', type=int)
    parser.add_argument('--foreground', action='store_true')
    parser.add_argument('--worker-manifest', type=Path)
    args = parser.parse_args()
    if args.worker_manifest:
        raise SystemExit(worker(args.worker_manifest))
    if not re.fullmatch(r'[a-zA-Z0-9_-]+', args.label) or args.ranks < 1:
        raise SystemExit('Invalid label or rank count')
    if any(value is not None and (not np.isfinite(value) or value <= 0)
           for value in [args.end_time, args.dt, args.write_interval]):
        raise SystemExit('Time and write interval overrides must be positive and finite')
    binary = args.binary.resolve()
    if not binary.is_file():
        raise SystemExit('Solver binary does not exist')
    directory = CASE/'results'
    inputs = directory/'inputs'/args.label
    status_path = directory/f'{args.label}_status.json'
    if inputs.exists() or (directory/args.label).exists() or status_path.exists():
        raise SystemExit('Refusing to overwrite an existing run; choose a new --label')
    inputs.mkdir(parents=True)
    for name in ['numerics', 'physics', 'fields']:
        shutil.copytree(CASE/name, inputs/name)
    (inputs/'mesh').mkdir()
    for name in ['cylinder_re3900.mesh', 'layout.npz', 'mesh_quality.json', 'reader_quality.json']:
        shutil.copy2(CASE/'mesh'/name, inputs/'mesh'/name)
    for name in ['case.bs', 'control.bs', 'output.bs', 'case_parameters.json']:
        shutil.copy2(CASE/name, inputs/name)
    text = (inputs/'output.bs').read_text()
    text = re.sub(r'(?m)^directory \S+$', 'directory results', text)
    text = re.sub(r'(?m)^timeName \S+$', 'timeName '+args.label, text)
    if args.write_interval is not None:
        text = re.sub(r'(?m)^writeInterval \S+$', 'writeInterval '+str(args.write_interval), text)
    (inputs/'output.bs').write_text(text)
    text = (inputs/'control.bs').read_text()
    for name, value in [('endTime', args.end_time), ('deltaT', args.dt)]:
        if value is not None:
            text = re.sub(r'(?m)^'+name+r' \S+$', name+' '+str(value), text)
    (inputs/'control.bs').write_text(text)
    # The solver requires a relative output directory without parent components.
    # Expose its local output at the public case result path with a one-way link.
    (directory/args.label).symlink_to(inputs/'results'/args.label, target_is_directory=True)
    control = dict(re.findall(r'(?m)^(\w+)\s+(\S+)$', (inputs/'control.bs').read_text()))
    parameters = json.loads((inputs/'case_parameters.json').read_text())
    if hashlib.sha256((inputs/'mesh/cylinder_re3900.mesh').read_bytes()).hexdigest() != parameters['mesh_sha256']:
        raise SystemExit('Mesh hash disagrees with generated parameters')
    settings = dict(re.findall(r'(?m)^([A-Za-z][A-Za-z0-9_.]*)\s+(\S+)$', (inputs/'numerics/solution.bs').read_text()))
    methods = dict(re.findall(r'(?m)^([A-Za-z][A-Za-z0-9_.]*)\s+(\S+)$', (inputs/'numerics/methods.bs').read_text()))
    manifest = {'label': args.label, 'ranks': args.ranks, 'threads_per_rank': 1,
                'input_case': str(inputs), 'parameters': parameters,
                'cells': json.loads((inputs/'mesh/mesh_quality.json').read_text())['cells'],
                'delta_t': float(control['deltaT']), 'end_time': float(control['endTime']),
                'time_method': methods.get('time'), 'numerical_methods': methods, 'solution_controls': settings,
                'git_head': subprocess.check_output(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip(),
                'binary_path': str(binary), 'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
                'guard_note': 'Finite fields and range guards are numerical diagnostics, not physical validation.'}
    manifest_path = inputs/'run_manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    if args.foreground:
        raise SystemExit(worker(manifest_path))
    else:
        with (directory/f'{args.label}_worker.log').open('w') as output:
            process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
                                        '--worker-manifest', str(manifest_path)],
                                       stdout=output, stderr=subprocess.STDOUT, env=ENV, start_new_session=True)
        print(json.dumps({'worker_pid': process.pid, 'status': str(status_path), 'input_case': str(inputs)}, indent=2))


if __name__ == '__main__':
    main()
