#!/usr/bin/env python3
"""Audit fixed-work engineering runs and reproduce scaling tables and figures.

Run only after the sequential numerical study has finished. The efficiency
baseline belongs to the same binary as its parallel points. Never substitute
an older/slower serial version. Calibration and profiling are kept separate.
"""
import csv
import json
from pathlib import Path
import re
import statistics

ROOT = Path(__file__).resolve().parent
SHAS = {
    'reference': '7a29412a574f912f6d67c7b3514b7aa15c60b854c222a455392ed977cb61a6b2',
    'optimized': '2114c30d7632f2fa500822cfafe2a134c69a5e280089bc312243ff3be2e7c53c',
}
CELLS = {'prism': 1024000, 'voronoi': 512000}
STEPS = {'prism': 10, 'voronoi': 20}
LABELS = {'prism': 'Warped triangular prisms (1,024,000 cells)',
          'voronoi': 'Voronoi prisms (512,000 cells)'}


def load_study(name):
    directory = ROOT / name
    metadata = json.loads((directory / 'study.json').read_text())
    for variant, item in metadata['variants'].items():
        assert item['sha256'] == SHAS[variant], (name, variant, item['sha256'])
    assert metadata['acceptedSteps'] == {case:STEPS[case] for case in metadata['cases']}
    records = json.loads((directory / 'records.json').read_text())
    result = []
    for run in records:
        if run['warmup']:
            continue
        assert run['validWork'] and run['measurementValid'] and run['sustainedSample']
        steps=STEPS[run['caseLabel']]
        assert run['status'] == 'converged' and run['returncode'] == 0 and run['workUnits'] == steps
        assert len(run['performance']) == run['ranks']
        assert all(p['reportedIterationIntervals'] == steps-1 for p in run['performance'])
        assert all(p['partition']['globalCells'] == CELLS[run['caseLabel']]
                   for p in run['performance'])
        assert abs(run['computeSeconds'] - max(p['iterationComputeSeconds']
                                              for p in run['performance'])) < 1e-9
        log = Path(run['performance'][0]['path']).parent.parent / 'stdout.log'
        lines = [line for line in log.read_text().splitlines()
                 if re.match(r'^PISO [0-9]+ .*accepted=true', line)]
        assert len(lines) == steps
        assert all('linear=ok' in line and 'converged=true' in line for line in lines)
        times = [float(re.search(r' time=([^ ]+)', line)[1]) for line in lines]
        assert all(abs(time - (i + 1)*.01) < 1e-12 for i, time in enumerate(times))
        result.append(dict(run, sourceStudy=name,
            correctors=[int(re.search(r' correctors=([0-9]+)', line)[1]) for line in lines],
            maximumMassError=max(float(re.search(r' mass=([^ ]+)', line)[1]) for line in lines),
            maximumCouplingResidual=max(float(re.search(r' rCoupling=([^ ]+)', line)[1])
                                        for line in lines)))
    return result, {name: item['manifest']['sha256'] for name, item in metadata['cases'].items()}


def med(runs, function):
    return statistics.median(function(run) for run in runs)


def table(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(records):
    rows = []
    for case in CELLS:
        for variant in SHAS:
            serial = [r for r in records if (r['caseLabel'], r['variant'], r['ranks'])
                      == (case, variant, 1)]
            assert serial
            t1 = med(serial, lambda r:r['computeSeconds'])
            for ranks in (1, 4, 8, 12):
                runs = [r for r in records if (r['caseLabel'], r['variant'], r['ranks'])
                        == (case, variant, ranks)]
                assert runs, (case, variant, ranks)
                times = [r['computeSeconds'] for r in runs]
                time = statistics.median(times)
                parts = runs[0]['performance']
                owned = [p['partition']['ownedCells'] for p in parts]
                faces = [p['partition']['ownedFaces'] for p in parts]
                local_faces = [p['partition']['localFaces'] for p in parts]
                critical = lambda r:max(r['performance'], key=lambda p:p['solverComputeSeconds'])
                pressure = lambda r:next(e for e in critical(r)['local']['equationSystems']
                                        if e['identity'].startswith('scalar:pressureCorrection:'))
                rows.append(dict(case=case, variant=variant, ranks=ranks, samples=len(runs),
                    acceptedSteps=STEPS[case],measuredIntervals=STEPS[case]-1,
                    cells=CELLS[case], computeMedian=time, computeMin=min(times), computeMax=max(times),
                    speedup=t1/time, efficiency=t1/(ranks*time), serialMedian=t1,
                    peakRssGiB=med(runs, lambda r:r['memory']['peakTreeRssKiB'])/1048576,
                    peakPssGiB=med(runs, lambda r:r['memory']['peakTreePssKiB'])/1048576,
                    peakSwapMiB=max(r['memory']['peakTreeSwapKiB'] for r in runs)/1024,
                    ownedMin=min(owned), ownedMax=max(owned),
                    ownedFaceImbalance=max(faces)/statistics.mean(faces),
                    localFaceImbalance=max(local_faces)/statistics.mean(local_faces),
                    totalGhostCells=sum(p['partition']['ghostCells'] for p in parts),
                    maxNeighbourRanks=max(p['partition']['neighbourRanks'] for p in parts),
                    krylovIterations=med(runs, lambda r:r['performance'][0]['krylovIterations']),
                    preconditionerSetups=med(runs, lambda r:r['performance'][0]['preconditionerSetups']),
                    linearSolveSeconds=med(runs, lambda r:critical(r)['local']['linearSolveSeconds']),
                    preconditionerSeconds=med(runs, lambda r:critical(r)['local']['preconditionerSeconds']),
                    pressureSolveSeconds=med(runs, lambda r:pressure(r)['solveSeconds']),
                    pressureSetupSeconds=med(runs, lambda r:pressure(r)['preconditionerSeconds']),
                    pressureShareOfFullCallback=med(runs, lambda r:
                        (pressure(r)['solveSeconds']+pressure(r)['preconditionerSeconds']) /
                         critical(r)['solverComputeSeconds']),
                    haloSeconds=med(runs, lambda r:critical(r)['local']['haloSeconds']),
                    reductionSeconds=med(runs, lambda r:critical(r)['local']['globalReductionSeconds']),
                    caseSetupSeconds=med(runs, lambda r:max(p['caseSetupSeconds'] for p in r['performance'])),
                    maximumMassError=max(r['maximumMassError'] for r in runs),
                    maximumCouplingResidual=max(r['maximumCouplingResidual'] for r in runs),
                    correctorSequences=json.dumps([r['correctors'] for r in runs])))
    comparisons = []
    for new in (r for r in rows if r['variant'] == 'optimized'):
        old = next(r for r in rows if (r['case'],r['variant'],r['ranks'])
                   == (new['case'],'reference',new['ranks']))
        comparisons.append(dict(case=new['case'], ranks=new['ranks'],
            referenceSamples=old['samples'], optimizedSamples=new['samples'],
            referenceSeconds=old['computeMedian'], optimizedSeconds=new['computeMedian'],
            versionSpeedup=old['computeMedian']/new['computeMedian'],
            timeChange=new['computeMedian']/old['computeMedian']-1,
            referenceEfficiency=old['efficiency'], optimizedEfficiency=new['efficiency'],
            referencePssGiB=old['peakPssGiB'], optimizedPssGiB=new['peakPssGiB'],
            referenceRssGiB=old['peakRssGiB'], optimizedRssGiB=new['peakRssGiB'],
            referenceSwapMiB=old['peakSwapMiB'], optimizedSwapMiB=new['peakSwapMiB']))
    return rows, comparisons


def figures(rows, comparisons):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = {'reference':'#9e5a19','optimized':'#126e82'}
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.9), layout='constrained')
    for column, case in enumerate(CELLS):
        for variant in SHAS:
            selected = [r for r in rows if r['case'] == case and r['variant'] == variant]
            x = [r['ranks'] for r in selected]
            y = [r['computeMedian'] for r in selected]
            errors = [[r['computeMedian']-r['computeMin'] for r in selected],
                      [r['computeMax']-r['computeMedian'] for r in selected]]
            axes[0,column].errorbar(x, y, yerr=errors, fmt='o-', capsize=4,
                color=colors[variant], label=variant.capitalize())
            axes[1,column].plot(x, [100*r['efficiency'] for r in selected], 'o-',
                color=colors[variant], label=variant.capitalize())
            if variant == 'optimized':
                axes[0,column].plot(x, [y[0]/p for p in x], ':', color=colors[variant],
                                   label='Optimized ideal T1/P')
        axes[0,column].set(title=LABELS[case]+f'\n{STEPS[case]} accepted steps, timing {STEPS[case]-1}', ylabel='Compute time (s)', yscale='log')
        axes[0,column].legend()
        axes[1,column].set(xlabel='MPI processes / physical cores', ylabel='Parallel efficiency (%)')
        axes[1,column].axhline(100, color='gray', ls=':', lw=1)
        for axis in axes[:,column]:
            axis.set_xticks([1,4,8,12]); axis.grid(alpha=.25)
    fig.suptitle('Fixed mesh and case-specific work at every process count\n'
                 'Each version uses its own serial time; setup and output excluded', fontsize=12)
    for suffix in ('png','svg','pdf'):
        fig.savefig(ROOT/('strong-scaling.'+suffix), dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), layout='constrained')
    for column, case in enumerate(CELLS):
        for variant in SHAS:
            selected = [r for r in rows if r['case'] == case and r['variant'] == variant]
            axes[column].plot([r['ranks'] for r in selected], [r['peakPssGiB'] for r in selected],
                              'o-', color=colors[variant], label=variant.capitalize())
        axes[column].set(title=LABELS[case], xlabel='MPI processes / physical cores',
                         ylabel='Whole-run process-tree peak PSS (GiB)')
        axes[column].set_xticks([1,4,8,12]); axes[column].grid(alpha=.25); axes[column].legend()
    fig.suptitle('PSS apportions shared pages; peaks include initialization and output')
    for suffix in ('png','svg','pdf'):
        fig.savefig(ROOT/('strong-memory.'+suffix), dpi=180)
    plt.close(fig)


def main():
    records, hashes = load_study('amg-strong-post-restart')
    for name in ('amg-voronoi20-strong-scaling','amg-prism-repeat','amg-voronoi20-repeat'):
        extra, extra_hashes = load_study(name)
        assert all(case not in hashes or hashes[case]==value for case,value in extra_hashes.items())
        hashes.update(extra_hashes)
        records.extend(extra)
    rows, comparisons = summarize(records)
    (ROOT/'strong-summary.json').write_text(json.dumps(dict(binarySha256=SHAS, caseHashes=hashes,
        metric='max_rank(iterationComputeSeconds)', acceptedSteps=STEPS,
        reportedIntervals={k:v-1 for k,v in STEPS.items()},
        memoryMetric='whole-lifetime simultaneous process-tree peak',
        scaling=rows, comparisons=comparisons),indent=2)+'\n')
    table(ROOT/'strong-scaling.csv', rows)
    table(ROOT/'strong-comparison.csv', comparisons)
    figures(rows, comparisons)


if __name__ == '__main__':
    main()
