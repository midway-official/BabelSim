#!/usr/bin/env python3
"""Compare finite startup snapshots; do not claim time-step independence."""
from pathlib import Path
from datetime import datetime, timezone
import json
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from postprocess import CASE, Geometry, NAMES, field, log_metrics

LABELS = {'reference_dt0p001': 'run', 'dt0p005': 'trial_dt0p005', 'dt0p01': 'trial_dt0p01_200'}
OUT = CASE/'validation/timestep_study'


def times(directory):
    result = {}
    if directory.exists():
        for path in directory.iterdir():
            try:
                t = float(path.name)
            except ValueError:
                continue
            if path.is_dir():
                result[round(t, 8)] = path
    return result


def norm(value, weight):
    return float(np.sqrt(np.sum(weight[:, None]*value**2)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    geo = Geometry(CASE/'results/inputs/run')
    area = np.load(CASE/'results/inputs/run/mesh/layout.npz')['areas']
    parameters = json.loads((CASE/'case_parameters.json').read_text())
    wake = (geo.c[:, 0] > -.6) & (geo.c[:, 0] < 5) & (abs(geo.c[:, 1]) < 2)
    series = {key: times(CASE/'results'/label) for key, label in LABELS.items()}
    common = sorted(set.intersection(*(set(paths) for paths in series.values())))
    reports = []
    force_history = {key: [] for key in LABELS}
    for t in common:
        try:
            arrays = {key: {name: field(paths[t], name, geo.n) for name in NAMES} for key, paths in series.items()}
        except (FileNotFoundError, ValueError):
            continue
        reference = arrays['reference_dt0p001']
        row = {'time': t, 'cases': {}}
        for key, values in arrays.items():
            forces, yplus = geo.wall_metrics(values['U'], values['p'], parameters)
            force_history[key].append([t, *forces])
            results = {'force_estimate': forces.tolist(), 'max_yplus_estimate': float(yplus.max())}
            if key != 'reference_dt0p001':
                results['differences'] = {}
                for name in NAMES:
                    difference = values[name]-reference[name]
                    results['differences'][name] = {
                        'relative_L2_cells_domain': float(np.linalg.norm(difference)/np.linalg.norm(reference[name])),
                        'relative_L2_cells_near_wake': float(np.linalg.norm(difference[wake])/np.linalg.norm(reference[name][wake])),
                        'relative_L2_area_weighted_domain': norm(difference, area)/norm(reference[name], area),
                        'relative_L2_area_weighted_near_wake': norm(difference[wake], area[wake])/norm(reference[name][wake], area[wake])}
            row['cases'][key] = results
        reports.append(row)
    if not reports:
        raise RuntimeError('No complete common time snapshots yet')
    # Compare field errors only at common times. Sample wall-force histories at
    # every saved time up to that same endpoint, so isolated pressure excursions
    # are not hidden by the coarser intersection of output intervals.
    end = reports[-1]['time']
    force_history = {key: [] for key in LABELS}
    for key, paths in series.items():
        for t, path in sorted(paths.items()):
            if t > end:
                break
            try:
                u, p = (field(path, name, geo.n) for name in ('U', 'p'))
            except (FileNotFoundError, ValueError):
                continue
            forces, _ = geo.wall_metrics(u, p, parameters)
            force_history[key].append([t, *forces])
    states = {}
    for key, label in LABELS.items():
        status = json.loads((CASE/'results'/(label+'_status.json')).read_text())
        log = CASE/'results'/('solver.log' if label == 'run' else label+'_solver.log')
        records = [dict(re.findall(r'(\w+)=(\S+)', line)) for line in log.read_text().splitlines()
                   if line.startswith('kOmega iteration=')]
        states[key] = {'status': status['status'], 'delta_t': status['delta_t'], 'ranks': status['ranks'],
                       'solver_exit_code': status.get('solver_exit_code'), 'postprocess_exit_code': status.get('postprocess_exit_code'),
                       'solver_log': log_metrics(log),
                       'min_k_raw': min(float(r['kMinRaw']) for r in records),
                       'min_omega_raw': min(float(r['omegaMinRaw']) for r in records),
                       'max_omega_bounded_cells': max(float(r['omegaBoundedCells']) for r in records)}
    report = {'updated_at': datetime.now(timezone.utc).isoformat(),
              'simulation': parameters['simulation'], 'Re_D': parameters['Re_D'],
              'near_wake_region': '-0.6 < x/D < 5; abs(y/D) < 2',
              'normalization': 'relative L2 against dt=0.001 at identical saved times; area weighting also reported',
              'evidence_limit': 'startup comparison; dt=0.001 is a reference, not an exact solution; no developed wake statistics or grid independence established',
              'rank_note': 'reference uses 8 ranks; trials use 4 ranks, same geometry and initial fields; trials allow 200 turbulence iterations at the same 1e-6 tolerance',
              'runs': states, 'times': reports,
              'saved_force_histories': force_history}
    (OUT/'comparison.json').write_text(json.dumps(report, indent=2)+'\n')
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), constrained_layout=True)
    colors = {'reference_dt0p001': 'black', 'dt0p005': 'tab:blue', 'dt0p01': 'tab:orange'}
    labels = {'reference_dt0p001': 'dt=0.001 reference', 'dt0p005': 'dt=0.005', 'dt0p01': 'dt=0.01'}
    for key in ('dt0p005', 'dt0p01'):
        ts = [r['time'] for r in reports]
        for ax, name in zip(axes[:2], ['U', 'p']):
            ax.plot(ts, [100*r['cases'][key]['differences'][name]['relative_L2_area_weighted_near_wake'] for r in reports],
                    marker='o', ms=3, label=labels[key], color=colors[key])
    axes[0].set(ylabel='Velocity difference (%)', title='Near-wake relative L2, cell-area weighted')
    axes[1].set(ylabel='Pressure difference (%)')
    for key, points in force_history.items():
        a = np.asarray(points)
        axes[2].plot(a[:, 0], a[:, 1], marker='o', ms=3, label=labels[key], color=colors[key])
    axes[2].set(ylabel='Cd estimate', title='Wall pressure + first-cell molecular shear; not a validated force integral')
    for ax in axes:
        ax.set_xlabel('t U/D')
        ax.grid(alpha=.2)
        ax.legend()
    fig.suptitle('Re=3900, 2D URANS/BDF2 | startup time-step comparison')
    fig.savefig(OUT/'comparison.png', dpi=160)
    plt.close(fig)
    # The same wall locations and Cp definition expose circumferential pressure
    # oscillations directly, without interpreting convergence flags as accuracy.
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), constrained_layout=True)
    theta = np.degrees(np.arctan2(geo.wc[:, 1], geo.wc[:, 0]))
    order = np.argsort(theta)
    normalization = .5*parameters['density']*parameters['U_infinity']**2
    latest_pressures = {}
    for key, paths in series.items():
        p = field(paths[end], 'p', geo.n)[geo.wo, 0]/normalization
        latest_pressures[key] = p
        axes[0].plot(theta[order], p[order], label=labels[key], color=colors[key], lw=1)
    for key in ('dt0p005', 'dt0p01'):
        axes[1].plot(theta[order], (latest_pressures[key]-latest_pressures['reference_dt0p001'])[order],
                     label=labels[key], color=colors[key], lw=1)
    axes[0].set(ylabel='Cp (wall-adjacent pressure)', title=f'Identical saved time: t U/D={end:g}')
    axes[1].set(ylabel='Cp difference from dt=0.001')
    for ax in axes:
        ax.set(xlabel='Cylinder angle from +x (degrees)', xlim=(-180, 180))
        ax.grid(alpha=.2)
        ax.legend()
    fig.suptitle('Re=3900, 2D URANS | wall pressure estimate')
    fig.savefig(OUT/'wall_pressure.png', dpi=160)
    plt.close(fig)
    print(json.dumps({'latest_common_time': reports[-1]['time'], 'latest_comparison': reports[-1], 'runs': states}, indent=2))


if __name__ == '__main__':
    main()
