#!/usr/bin/env python3
"""Audit sustained SIMPLE throughput, separately from physical convergence."""
import csv
import json
from pathlib import Path
import re
import statistics

ROOT=Path(__file__).resolve().parent
SHAS={
    'reference':'7a29412a574f912f6d67c7b3514b7aa15c60b854c222a455392ed977cb61a6b2',
    'optimized':'2114c30d7632f2fa500822cfafe2a134c69a5e280089bc312243ff3be2e7c53c',
}


def load(name):
    folder=ROOT/name
    metadata=json.loads((folder/'study.json').read_text())
    assert metadata['throughput']=={'simple':384} and not metadata['acceptedSteps']
    for variant,item in metadata['variants'].items():
        assert item['sha256']==SHAS[variant]
    records=json.loads((folder/'records.json').read_text())
    result=[]
    for run in records:
        if run['warmup']:
            continue
        assert run['caseLabel']=='simple' and run['workUnits']==384
        assert run['status']=='maxIterations' and run['returncode']==2 and run['validWork']
        assert run['measurementValid'] and run['sustainedSample']
        assert len(run['performance'])==run['ranks']
        assert all(p['reportedIterationIntervals']==383 and
                   p['partition']['globalCells']==500094 for p in run['performance'])
        assert abs(run['computeSeconds']-max(p['iterationComputeSeconds']
                                            for p in run['performance']))<1e-9
        log=Path(run['performance'][0]['path']).parent.parent/'stdout.log'
        iterations=[line for line in log.read_text().splitlines()
                    if re.match(r'^SIMPLE [0-9]+ ',line)]
        assert iterations and max(int(line.split()[1]) for line in iterations)==384
        assert all('linear=ok' in line for line in iterations)
        result.append(dict(run,sourceStudy=name))
    return result,metadata['cases']['simple']['manifest']['sha256']


def med(runs,function):
    return statistics.median(function(run) for run in runs)


def main():
    records,case_hash=load('simple-strong-scaling')
    if (ROOT/'simple-repeat').exists():
        extra,extra_hash=load('simple-repeat')
        assert case_hash==extra_hash
        records.extend(extra)
    rows=[]
    for variant in SHAS:
        serial=[r for r in records if r['variant']==variant and r['ranks']==1]
        assert serial
        t1=med(serial,lambda r:r['computeSeconds'])
        for ranks in (1,4,8,12):
            runs=[r for r in records if r['variant']==variant and r['ranks']==ranks]
            assert runs
            values=[r['computeSeconds'] for r in runs]; time=statistics.median(values)
            parts=runs[0]['performance'];faces=[p['partition']['ownedFaces'] for p in parts]
            local_faces=[p['partition']['localFaces'] for p in parts]
            rows.append(dict(variant=variant,ranks=ranks,samples=len(runs),cells=500094,
                outerIterations=384,measuredIntervals=383,computeMedian=time,
                computeMin=min(values),computeMax=max(values),serialMedian=t1,
                speedup=t1/time,efficiency=t1/(time*ranks),
                peakPssGiB=med(runs,lambda r:r['memory']['peakTreePssKiB'])/1048576,
                peakRssGiB=med(runs,lambda r:r['memory']['peakTreeRssKiB'])/1048576,
                peakSwapMiB=max(r['memory']['peakTreeSwapKiB'] for r in runs)/1024,
                totalGhostCells=sum(p['partition']['ghostCells'] for p in parts),
                ownedFaceImbalance=max(faces)/statistics.mean(faces),
                localFaceImbalance=max(local_faces)/statistics.mean(local_faces),
                maxNeighbourRanks=max(p['partition']['neighbourRanks'] for p in parts),
                krylovIterations=med(runs,lambda r:r['performance'][0]['krylovIterations']),
                preconditionerSetups=med(runs,lambda r:r['performance'][0]['preconditionerSetups']),
                preconditionerSeconds=med(runs,lambda r:max(p['local']['preconditionerSeconds']
                                                          for p in r['performance'])),
                linearSolveSeconds=med(runs,lambda r:max(p['local']['linearSolveSeconds']
                                                      for p in r['performance']))))
    comparisons=[]
    for new in (r for r in rows if r['variant']=='optimized'):
        old=next(r for r in rows if r['variant']=='reference' and r['ranks']==new['ranks'])
        comparisons.append(dict(ranks=new['ranks'],referenceSeconds=old['computeMedian'],
            optimizedSeconds=new['computeMedian'],versionSpeedup=old['computeMedian']/new['computeMedian'],
            referenceEfficiency=old['efficiency'],optimizedEfficiency=new['efficiency']))
    data=dict(binarySha256=SHAS,caseHash=case_hash,metric='max_rank(iterationComputeSeconds)',
        purpose='SIMPLE throughput, not time to physical steady convergence',
        scaling=rows,comparisons=comparisons)
    (ROOT/'simple-summary.json').write_text(json.dumps(data,indent=2)+'\n')
    for filename,table in [('simple-scaling.csv',rows),('simple-comparison.csv',comparisons)]:
        with (ROOT/filename).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11.2,4.7),layout='constrained')
    colors={'reference':'#9e5a19','optimized':'#126e82'}
    for variant in SHAS:
        chosen=[r for r in rows if r['variant']==variant]; x=[r['ranks'] for r in chosen]
        y=[r['computeMedian'] for r in chosen]
        error=[[r['computeMedian']-r['computeMin'] for r in chosen],
               [r['computeMax']-r['computeMedian'] for r in chosen]]
        axes[0].errorbar(x,y,yerr=error,fmt='o-',capsize=4,color=colors[variant],label=variant.capitalize())
        axes[1].plot(x,[100*r['efficiency'] for r in chosen],'o-',color=colors[variant],label=variant.capitalize())
        if variant=='optimized':
            axes[0].plot(x,[y[0]/p for p in x],':',color=colors[variant],label='Optimized ideal T1/P')
    axes[0].set(ylabel='Compute time (s)',yscale='log')
    axes[1].set(ylabel='Parallel efficiency (%)');axes[1].axhline(100,color='gray',ls=':',lw=1)
    for axis in axes:
        axis.set_xticks([1,4,8,12]);axis.set_xlabel('MPI processes / physical cores');axis.grid(alpha=.25);axis.legend()
    fig.suptitle('SIMPLE throughput: 500,094 cells, 384 outer iterations (timing 383)\n'
                 'Fixed work, setup and output excluded; not time to steady convergence',fontsize=12)
    for suffix in ('png','svg','pdf'):
        fig.savefig(ROOT/('simple-scaling.'+suffix),dpi=180)
    plt.close(fig)


if __name__=='__main__':
    main()
