#!/usr/bin/env python3
"""Reproduce final-version scaling, reference comparisons and exported figures."""
from pathlib import Path
import csv
import json
import statistics

ROOT = Path(__file__).resolve().parent
FINAL_SHA = '7a29412a574f912f6d67c7b3514b7aa15c60b854c222a455392ed977cb61a6b2'
REFERENCE_SHA = 'd2b572ecd1446fd309af35b02918f4c741f3e0ef1ddd90161cf6b5541a1c8a14'


def collect(directory, variant, expected_sha):
    meta = json.loads((ROOT / directory / 'study.json').read_text())
    assert meta['variants'][variant]['sha256'] == expected_sha
    records = json.loads((ROOT / directory / 'records.json').read_text())
    result = []
    for record in records:
        if record['variant'] != variant or record['warmup']:
            continue
        assert record['validWork'] and record['measurementValid'] and record['sustainedSample']
        expected = 128 if record['caseLabel'] == 'simple' else 6
        assert record['workUnits'] == expected
        assert all(p['reportedIterationIntervals'] == expected - 1 for p in record['performance'])
        record = dict(record, sourceStudy=directory)
        result.append(record)
    return result, {key: value['manifest']['sha256'] for key, value in meta['cases'].items()}


def median(runs, function):
    return statistics.median(function(run) for run in runs)


def main():
    final, hashes = collect('final-strong-scaling', 'final', FINAL_SHA)
    reference, reference_hashes = collect('memory-partition-ab', 'rhs-budget', REFERENCE_SHA)
    assert hashes == reference_hashes
    repeat_path = ROOT / 'final-repeat-12/records.json'
    if repeat_path.exists():
        extra, extra_hashes = collect('final-repeat-12', 'final', FINAL_SHA)
        old_extra, old_hashes = collect('final-repeat-12', 'rhs-budget', REFERENCE_SHA)
        assert hashes == extra_hashes == old_hashes
        final.extend(extra)
        reference.extend(old_extra)
    rows, comparisons = [], []
    for case in ('simple', 'piso'):
        serial = [run for run in final if run['caseLabel'] == case and run['ranks'] == 1]
        t1 = median(serial, lambda r: r['computeSeconds'])
        for ranks in (1, 4, 8, 12):
            runs = [r for r in final if r['caseLabel'] == case and r['ranks'] == ranks]
            assert runs
            times = [r['computeSeconds'] for r in runs]
            time = statistics.median(times)
            parts = runs[0]['performance']
            faces = [p['partition']['localFaces'] for p in parts]
            def pressure_time(run):
                critical = max(run['performance'], key=lambda p:p['solverComputeSeconds'])
                equation = next(e for e in critical['local']['equationSystems']
                                if e['identity'].startswith('scalar:pressureCorrection:'))
                time = equation['preconditionerSeconds'] + equation['solveSeconds']
                return time, time / critical['solverComputeSeconds']
            row = dict(case=case, ranks=ranks, samples=len(runs), computeMedian=time,
                computeMin=min(times), computeMax=max(times), speedup=t1/time,
                efficiency=t1/(time*ranks),
                peakRssGiB=median(runs, lambda r: r['memory']['peakTreeRssKiB'])/1048576,
                peakPssGiB=median(runs, lambda r: r['memory']['peakTreePssKiB'])/1048576,
                peakSwapMiB=max(r['memory']['peakTreeSwapKiB'] for r in runs)/1024,
                totalGhostCells=sum(p['partition']['ghostCells'] for p in parts),
                localFaceImbalance=max(faces)/statistics.mean(faces),
                maxNeighbourRanks=max(p['partition']['neighbourRanks'] for p in parts),
                kspIterations=median(runs, lambda r: r['performance'][0]['krylovIterations']),
                preconditionerSetups=median(runs, lambda r: r['performance'][0]['preconditionerSetups']),
                preconditionerSeconds=median(runs, lambda r: max(p['local']['preconditionerSeconds'] for p in r['performance'])),
                linearSolveSeconds=median(runs, lambda r: max(p['local']['linearSolveSeconds'] for p in r['performance'])),
                pressureSetupAndSolveSeconds=median(runs, lambda r: pressure_time(r)[0]),
                pressureShareOfFullCallback=median(runs, lambda r: pressure_time(r)[1]),
                caseSetupSeconds=median(runs, lambda r: max(p['caseSetupSeconds'] for p in r['performance'])))
            rows.append(row)
            old = [r for r in reference if r['caseLabel'] == case and r['ranks'] == ranks]
            if old:
                old_time = median(old, lambda r: r['computeSeconds'])
                old_rss = median(old, lambda r: r['memory']['peakTreeRssKiB'])/1048576
                old_pss = median(old, lambda r: r['memory']['peakTreePssKiB'])/1048576
                comparisons.append(dict(case=case, ranks=ranks, referenceSamples=len(old),
                    finalSamples=len(runs), referenceSeconds=old_time, finalSeconds=time,
                    computeChange=time/old_time-1, referenceRssGiB=old_rss,
                    finalRssGiB=row['peakRssGiB'], rssChange=row['peakRssGiB']/old_rss-1,
                    referencePssGiB=old_pss, finalPssGiB=row['peakPssGiB'],
                    pssChange=row['peakPssGiB']/old_pss-1,
                    referenceSwapMiB=max(r['memory']['peakTreeSwapKiB'] for r in old)/1024,
                    finalSwapMiB=row['peakSwapMiB']))
    data = dict(finalBinarySha256=FINAL_SHA, referenceBinarySha256=REFERENCE_SHA,
                caseHashes=hashes, metric='max_rank(iterationComputeSeconds)',
                memoryMetric='whole-process-lifetime simultaneous process-tree peak',
                scaling=rows, comparisons=comparisons)
    (ROOT/'final-summary.json').write_text(json.dumps(data, indent=2)+'\n')
    for name, table in [('final-scaling', rows), ('final-comparison', comparisons)]:
        with (ROOT/(name+'.csv')).open('w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=list(table[0]))
            writer.writeheader(); writer.writerows(table)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), layout='constrained')
    for column, case in enumerate(('simple', 'piso')):
        selected = [r for r in rows if r['case'] == case]
        x = [r['ranks'] for r in selected]
        y = [r['computeMedian'] for r in selected]
        axes[0,column].plot(x, y, 'o-', label='Measured')
        axes[0,column].plot(x, [y[0]/p for p in x], '--', color='gray', label='Ideal T1/P')
        axes[0,column].set(title=case.upper(), ylabel='Compute time (s)')
        axes[0,column].legend()
        axes[1,column].plot(x, [100*r['efficiency'] for r in selected], 'o-')
        axes[1,column].set(xlabel='MPI ranks (physical cores)', ylabel='Parallel efficiency (%)', ylim=(0,105))
        for axis in axes[:,column]:
            axis.set_xticks(x); axis.grid(alpha=.25)
    fig.suptitle('500,094 cells; same final binary and fixed work at all rank counts')
    for suffix in ('png','svg','pdf'): fig.savefig(ROOT/('final-scaling.'+suffix), dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), layout='constrained')
    labels = [r['case'].upper()+' / '+str(r['ranks']) for r in comparisons]
    positions = list(range(len(comparisons)))
    for axis, metric, title in zip(axes, ('Rss','Pss'), ('RSS (includes repeated shared pages)', 'PSS (apportions shared pages)')):
        axis.bar([p-.18 for p in positions], [r['reference'+metric+'GiB'] for r in comparisons], .36, label='Before memory/partition')
        axis.bar([p+.18 for p in positions], [r['final'+metric+'GiB'] for r in comparisons], .36, label='Final')
        axis.set_xticks(positions, labels); axis.set(ylabel='Whole-run peak (GiB)', title=title)
        axis.grid(axis='y', alpha=.25); axis.legend()
    for suffix in ('png','svg','pdf'): fig.savefig(ROOT/('final-memory.'+suffix), dpi=160)
    plt.close(fig)
    print(json.dumps(data, indent=2))


if __name__ == '__main__':
    main()
