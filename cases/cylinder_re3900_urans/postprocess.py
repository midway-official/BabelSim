#!/usr/bin/env python3
"""Inspect complete owned-cell snapshots of the Re=3900 2D URANS case.

Signed vorticity uses independent least-squares gradients. Wall forces and y+
use a first-cell tangential-gradient estimate, not a mesh-converged wall stress.
"""
from pathlib import Path
import argparse
import json
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation, PillowWriter

CASE = Path(__file__).resolve().parent
NAMES = ('U', 'p', 'k', 'omega', 'mut')


def field(path, name, count):
    files = sorted(path.glob(f'rank-*/{name}.csv'))
    if not files:
        files = sorted(path.glob(f'{name}.csv'))
    if not files:
        raise FileNotFoundError(f'Missing {name} at {path}')
    values = np.concatenate([np.loadtxt(p, delimiter=',', skiprows=1, ndmin=2) for p in files])
    values = values[np.argsort(values[:, 0])]
    if len(values) != count or not np.array_equal(values[:, 0], np.arange(count)):
        raise ValueError(f'Incomplete or duplicate global IDs for {name} at {path}')
    if not np.isfinite(values).all():
        raise ValueError(f'Nonfinite {name} at {path}')
    return values[:, 4:]


class Geometry:
    def __init__(self, case=CASE):
        data = np.load(case/'mesh/layout.npz')
        self.xy, self.tri, self.c = data['xy'], data['tri'], data['centres']
        self.n = len(self.c)
        edges = data['edge_cells']
        self.o, self.nb = edges[:, 0], edges[:, 1]
        wall = data['wall']
        self.wo = wall[:, 2]
        self.wc = .5*(self.xy[wall[:, 0]]+self.xy[wall[:, 1]])
        d = self.xy[wall[:, 1]]-self.xy[wall[:, 0]]
        self.wlen = np.linalg.norm(d, axis=1)
        self.wn = np.column_stack([d[:, 1], -d[:, 0]])/self.wlen[:, None]
        self.wdist = np.sum((self.wc-self.c[self.wo])*self.wn, axis=1)
        self.i = np.r_[self.o, self.nb, self.wo]
        self.j = np.r_[self.nb, self.o]
        self.delta = np.vstack([self.c[self.nb]-self.c[self.o],
                                self.c[self.o]-self.c[self.nb], self.wc-self.c[self.wo]])
        self.w = 1/np.sum(self.delta*self.delta, axis=1)
        self.matrix = np.zeros((self.n, 2, 2))
        np.add.at(self.matrix, self.i, self.w[:, None, None]*self.delta[:, :, None]*self.delta[:, None, :])
        if not np.all(np.linalg.det(self.matrix) > 0):
            raise ValueError('Singular independent gradient stencil')

    def gradient(self, u, wall_values=None):
        if wall_values is None:
            wall_values = np.zeros((len(self.wo), u.shape[1]))
        diff = np.vstack([u[self.j]-u[self.i[:len(self.j)]], wall_values-u[self.wo]])
        rhs = np.zeros((self.n, 2, u.shape[1]))
        np.add.at(rhs, self.i, self.w[:, None, None]*self.delta[:, :, None]*diff[:, None, :])
        return np.linalg.solve(self.matrix, rhs)

    def check_affine(self):
        a = np.array([[1.2, -.4], [.7, 2.1]])
        error = float(np.max(abs(self.gradient(self.c@a+.3, self.wc@a+.3)-a)))
        if error > 1e-10:
            raise ValueError(f'Independent gradient affine error {error}')
        return error

    def wall_metrics(self, u, p, parameters):
        nu, rho, speed, diameter = (parameters[k] for k in ['nu', 'density', 'U_infinity', 'diameter'])
        v = u[self.wo, :2]
        tangential = v-np.sum(v*self.wn, axis=1)[:, None]*self.wn
        shear = rho*nu*np.linalg.norm(tangential, axis=1)/self.wdist
        yplus = self.wdist*np.sqrt(shear/rho)/nu
        traction = p[self.wo, 0, None]*self.wn+rho*nu*tangential/self.wdist[:, None]
        coefficient = np.sum(traction*self.wlen[:, None], axis=0)/(.5*rho*speed**2*diameter)
        return coefficient, yplus


def log_metrics(log):
    records = []
    if log.exists():
        for line in log.read_text(errors='replace').splitlines():
            if re.match(r'^PISO \d+ time=', line):
                record = dict(re.findall(r'(\w+)=(\S+)', line))
                record['outer_index'] = int(line.split()[1])
                records.append(record)
    if not records:
        return {'reported_steps': 0}
    completed = [r for r in records if r.get('converged') == 'true']
    last_index = {r['time']: r['outer_index'] for r in records}
    result = {'reported_steps': len(completed), 'reported_outer_iterations': sum(last_index.values()),
              'reported_outer_records': len(records),
              'last_time': float(completed[-1]['time']) if completed else 0.,
              'last_iteration_time': float(records[-1]['time']),
              'last_iteration_converged': records[-1].get('converged') == 'true',
              'all_steps_accepted': all(r.get('accepted') == 'true' for r in records),
              'all_linear_ok': all(r.get('linear') == 'ok' for r in records),
              'all_turbulence_settled': bool(completed) and all(r.get('turbulenceSettled') == 'true' for r in completed),
              'all_completed_steps_converged': bool(completed)}
    for key in ['mass', 'CoMax', 'rCoupling']:
        result['max_'+key] = max(float(r[key]) for r in records)
    for key in ['rTurb', 'dTurb']:
        result['max_'+key] = max((float(r[key]) for r in completed), default=None)
        result['max_'+key+'_before_convergence'] = max(float(r[key]) for r in records)
    return result


def main():
    active_label = re.search(r'(?m)^timeName\s+(\S+)$', (CASE/'output.bs').read_text()).group(1)
    default_results = CASE/'results'/active_label
    default_analysis = CASE/'validation'/active_label/'analysis' if active_label != 'run' else CASE/'validation/production'
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, default=default_results)
    parser.add_argument('--out', type=Path, default=default_analysis)
    parser.add_argument('--log', type=Path, default=CASE/'results'/('solver.log' if active_label == 'run' else active_label+'_solver.log'))
    parser.add_argument('--case', type=Path, default=CASE)
    parser.add_argument('--latest-only', action='store_true')
    parser.add_argument('--require-complete', action='store_true')
    parser.add_argument('--gif', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    geometry_case = args.case
    if geometry_case == CASE and args.results == default_results and (CASE/'results/inputs'/active_label).exists():
        geometry_case = CASE/'results/inputs'/active_label
    geo = Geometry(geometry_case)
    parameters = json.loads((geometry_case/'case_parameters.json').read_text())
    paths = []
    for path in args.results.iterdir():
        try:
            t = float(path.name)
        except ValueError:
            continue
        if path.is_dir():
            paths.append((t, path))
    paths.sort()
    skipped = []
    if args.latest_only:
        candidates = list(reversed(paths))
    else:
        candidates = paths
    probes_xy = np.array([[2., .5], [4., .5], [6., 0.], [10., .5]])
    probes = np.argmin(np.sum((geo.c[:, None, :]-probes_xy[None, :, :])**2, axis=2), axis=0)
    rows, complete = [], []
    latest = None
    yplus_max = 0.
    ranges = {name: [float('inf'), float('-inf')] for name in ['k', 'omega', 'mut']}
    for t, path in candidates:
        try:
            fields = {name: field(path, name, geo.n) for name in NAMES}
        except (FileNotFoundError, ValueError) as exc:
            skipped.append({'time': t, 'reason': str(exc)})
            continue
        u, p = fields['U'], fields['p']
        coefficient, yplus = geo.wall_metrics(u, p, parameters)
        yplus_max = max(yplus_max, float(yplus.max()))
        rows.append([t, *u[probes, 1], *coefficient, float(np.linalg.norm(u, axis=1).max()),
                     float(abs(p).max()), float(yplus.max()), float(fields['mut'].max()/parameters['dynamicViscosity'])])
        for name in ranges:
            ranges[name][0] = min(ranges[name][0], float(fields[name].min()))
            ranges[name][1] = max(ranges[name][1], float(fields[name].max()))
        latest = (t, fields, yplus)
        complete.append((t, path))
        if args.latest_only:
            break
    if latest is None:
        raise RuntimeError('No complete finite U/p/k/omega/mut snapshots')
    if skipped and args.require_complete:
        raise RuntimeError(f'Incomplete or invalid finished-run snapshots: {skipped}')
    rows = np.asarray(rows)
    np.savetxt(args.out/'probes_forces.csv', rows, delimiter=',', comments='',
               header='time,Uy_2_0.5,Uy_4_0.5,Uy_6_0,Uy_10_0.5,Cd_estimate,Cl_estimate,max_speed,max_abs_p,max_yplus_estimate,max_mut_over_mu')
    t, fields, yplus = latest
    gradient = geo.gradient(fields['U'][:, :2])
    vorticity = gradient[:, 0, 1]-gradient[:, 1, 0]
    fig, axes = plt.subplots(4, 1, figsize=(13, 12), constrained_layout=True)
    variables = [(vorticity, 'Signed vorticity: dUy/dx - dUx/dy', 'RdBu_r', {'vmin': -4, 'vmax': 4}),
                 (np.linalg.norm(fields['U'], axis=1), 'Speed / U_infinity', 'viridis', {'vmin': 0, 'vmax': 2}),
                 (np.maximum(fields['k'][:, 0], 1e-8), 'Modelled turbulent kinetic energy k', 'magma', {'norm': LogNorm(1e-8, .1)}),
                 (np.maximum(fields['mut'][:, 0]/parameters['dynamicViscosity'], 1e-4), 'Turbulent / molecular viscosity', 'viridis', {'norm': LogNorm(1e-4, 1e3)})]
    for ax, (values, title, cmap, limits) in zip(axes, variables):
        artist = ax.tripcolor(geo.xy[:, 0], geo.xy[:, 1], geo.tri, facecolors=values, cmap=cmap, **limits)
        ax.add_patch(plt.Circle((0, 0), .5, color='.65', zorder=3))
        ax.set(xlim=(-2, 16), ylim=(-3, 3), aspect='equal', xlabel='x/D', ylabel='y/D', title=title)
        fig.colorbar(artist, ax=ax, extend='both')
    fig.suptitle(f'Re=3900 | 2D URANS, Wilcox1988m k-omega | t U/D={t:g}')
    fig.savefig(args.out/'flow_latest.png', dpi=160)
    plt.close(fig)
    angle = np.degrees(np.arctan2(geo.wc[:, 1], geo.wc[:, 0]))
    order = np.argsort(angle)
    np.savetxt(args.out/'wall_yplus.csv', np.column_stack([angle[order], yplus[order]]),
               delimiter=',', header='theta_deg,yplus_first_cell_estimate', comments='')
    spectrum = {'status': 'insufficient developed-flow sampling'}
    if not args.latest_only:
        fig, axes = plt.subplots(3, 1, figsize=(11, 9), constrained_layout=True)
        for i, xy in enumerate(probes_xy):
            axes[0].plot(rows[:, 0], rows[:, i+1], label=f'{tuple(xy)}')
        axes[0].set(ylabel='Uy / U_infinity', title='Wake probes')
        axes[0].legend(ncol=4)
        axes[1].plot(rows[:, 0], rows[:, 5], label='Cd estimate')
        axes[1].plot(rows[:, 0], rows[:, 6], label='Cl estimate')
        axes[1].legend()
        axes[1].set(ylabel='Force coefficient', title='Wall pressure + first-cell molecular shear estimate')
        sample = rows[rows[:, 0] >= max(40., .5*rows[-1, 0])]
        if len(sample) >= 64 and sample[-1, 0]-sample[0, 0] >= 25:
            ts, values = sample[:, 0], sample[:, 2]
            interval = float(np.median(np.diff(ts)))
            if np.allclose(np.diff(ts), interval, atol=1e-8):
                signal = values-np.polyval(np.polyfit(ts, values, 1), ts)
                power = abs(np.fft.rfft(signal*np.hanning(len(signal))))**2
                frequency = np.fft.rfftfreq(len(signal), interval)
                valid = (frequency > .05) & (frequency < .5)
                peak = np.flatnonzero(valid)[np.argmax(power[valid])]
                axes[2].plot(frequency, power)
                axes[2].set(xlim=(0, .6), xlabel='f D/U_infinity', ylabel='Windowed spectral power')
                spectrum = {'status': 'finite-window estimate, requires time/mesh sensitivity for validation',
                            'St_peak': float(frequency[peak]), 'bin_width': float(frequency[1]),
                            'window': [float(ts[0]), float(ts[-1])], 'probe_Uy_rms': float(np.std(values)),
                            'Cd_mean_estimate': float(sample[:, 5].mean()), 'Cl_rms_estimate': float(sample[:, 6].std())}
        else:
            axes[2].text(.5, .5, 'More developed-flow samples needed for shedding spectrum',
                         ha='center', transform=axes[2].transAxes)
        for ax in axes:
            ax.grid(alpha=.2)
        axes[0].set_xlabel('t U/D')
        axes[1].set_xlabel('t U/D')
        fig.savefig(args.out/'probes_spectrum.png', dpi=150)
        plt.close(fig)
    report = {'cells': geo.n, 'simulation': parameters['simulation'], 'model': parameters['turbulence_model'],
              'Re_D': parameters['Re_D'], 'complete_finite_snapshots': len(complete),
              'last_saved_time': t, 'skipped_snapshots': skipped, 'field_ranges_over_inspected_snapshots': ranges,
              'latest_max_speed': float(np.linalg.norm(fields['U'], axis=1).max()),
              'latest_pressure_range': [float(fields['p'].min()), float(fields['p'].max())],
              'latest_signed_vorticity_range': [float(vorticity.min()), float(vorticity.max())],
              'max_yplus_over_inspected_snapshots': yplus_max, 'latest_yplus_p95': float(np.percentile(yplus, 95)),
              'wall_shear_method': 'first-cell tangential velocity / wall-normal distance; molecular viscosity',
              'affine_gradient_max_error': geo.check_affine(), 'actual_probe_xy': geo.c[probes].tolist(),
              'solver_log': log_metrics(args.log), 'spectrum': spectrum,
              'physical_validation': '2D URANS model calculation; 3D turbulent wake and mesh/time independence are not established'}
    (args.out/'diagnostics.json').write_text(json.dumps(report, indent=2)+'\n')
    if args.gif:
        chosen = complete[-100:][::2]
        fig, ax = plt.subplots(figsize=(11, 4), constrained_layout=True)
        artist = ax.tripcolor(geo.xy[:, 0], geo.xy[:, 1], geo.tri, facecolors=vorticity, cmap='RdBu_r', vmin=-4, vmax=4)
        ax.add_patch(plt.Circle((0, 0), .5, color='.65', zorder=3))
        ax.set(xlim=(-2, 16), ylim=(-3, 3), aspect='equal', xlabel='x/D', ylabel='y/D')
        fig.colorbar(artist, ax=ax, extend='both')
        def update(i):
            snapshot_time, path = chosen[i]
            grad = geo.gradient(field(path, 'U', geo.n)[:, :2])
            artist.set_array(grad[:, 0, 1]-grad[:, 1, 0])
            ax.set_title(f'Re=3900, 2D URANS signed vorticity, t U/D={snapshot_time:g}')
            return artist,
        FuncAnimation(fig, update, frames=len(chosen)).save(args.out/'vorticity.gif', writer=PillowWriter(fps=8), dpi=110)
        plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
