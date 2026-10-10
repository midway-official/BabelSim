#!/usr/bin/env python3
"""Sequential Python study of FVM mesh partition quality and load balance.

Graph creation, I/O, Numba compilation, quality evaluation and plotting are
outside partition timers. No solver is run. Solver speed must be measured after
an algorithm has been selected and integrated through the existing mesh API.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import gc
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import shutil
import threading
import time
import numpy as np
import psutil
import scipy
import numba

from partition_model_graphs import FAMILIES, Graph, make_graph, native_fixture
from partition_model_algorithms import (ALGORITHMS, Metis, components_per_part,
    geometric, graph_growth, halo_masks, quotas, rebalance, warm_kernels)

ROOT=Path(__file__).resolve().parents[1]
SOURCES={
 'METIS manual':'https://karypis.github.io/glaros/files/sw/metis/manual.pdf',
 'Graph partitioning research':'https://papers.karypis.org/glaros/projects/gp.html',
 'Zoltan RCB':'https://sandialabs.github.io/Zoltan/ug_html/ug_alg_rcb.html',
 'Zoltan RIB':'https://sandialabs.github.io/Zoltan/ug_html/ug_alg_rib.html',
 'Zoltan geometry guide':'https://sandialabs.github.io/Zoltan/ug_html/ug_alg_geom.html'}


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1<<20),b''): digest.update(block)
    return digest.hexdigest()


class MemorySampler:
    def __init__(self):
        self.process=psutil.Process(); self.before=self.process.memory_info().rss
        self.peak=self.before; self.stop_event=threading.Event()
        self.thread=threading.Thread(target=self.sample,daemon=True)
    def sample(self):
        while not self.stop_event.is_set():
            self.peak=max(self.peak,self.process.memory_info().rss)
            self.stop_event.wait(.02)
    def start(self): self.thread.start()
    def stop(self):
        self.peak=max(self.peak,self.process.memory_info().rss)
        self.stop_event.set(); self.thread.join()
        return dict(beforeRSSBytes=self.before,peakRSSBytes=self.peak,
                    sampledRSSIncrementBytes=max(0,self.peak-self.before))


def metrics(graph,owners,parts,layers=3):
    assert len(owners)==graph.cells and owners.min()>=0 and owners.max()<parts
    expected=quotas(graph.cells,parts)
    owned=np.bincount(owners,minlength=parts)
    assert np.array_equal(owned,expected), 'exact owned-cell quotas failed'
    rows=np.repeat(np.arange(graph.cells,dtype=np.int32),np.diff(graph.offsets))
    columns=graph.neighbours
    unique=rows<columns; u,v=rows[unique],columns[unique]
    cross=owners[u]!=owners[v]
    cut=int(cross.sum())
    cutloads=np.bincount(np.concatenate((owners[u[cross]],owners[v[cross]])),minlength=parts)
    masks,first=halo_masks(graph.offsets,graph.neighbours,owners,layers)
    physical=graph.face_counts-np.diff(graph.offsets)
    ghost=[]; faces=[]; firstghost=[]; neighbourcounts=[]
    work=np.bincount(owners,weights=graph.face_counts,minlength=parts).astype(np.int64)
    for part in range(parts):
        selected=(masks & (np.uint64(1)<<np.uint64(part)))!=0
        ghosts=selected & (owners!=part)
        ghost.append(int(ghosts.sum()))
        firstghost.append(int((((first & (np.uint64(1)<<np.uint64(part)))!=0) & (owners!=part)).sum()))
        faces.append(int((selected[u]|selected[v]).sum()+physical[selected].sum()))
        neighbourcounts.append(int(len(np.unique(owners[ghosts]))))
    localcells=owned+np.array(ghost)
    globalfaces=len(u)+int(physical.sum())
    result=dict(cutFaces=cut,cutFacesPerPart=cutloads.tolist(),ghostLayers=layers,
        ownedCellsPerPart=owned.tolist(),ghostCellsPerPart=ghost,firstLayerGhostsPerPart=firstghost,
        localCellsPerPart=localcells.tolist(),localFacesPerPart=faces,neighbourRanksPerPart=neighbourcounts,
        ownedFaceIncidencesPerPart=work.tolist(),
        ownedCellImbalance=float(owned.max()*parts/owned.sum()),
        ownedFaceWorkImbalance=float(work.max()*parts/work.sum()),
        localFaceImbalance=float(max(faces)*parts/sum(faces)),
        localCellImbalance=float(localcells.max()*parts/localcells.sum()),
        totalGhostCells=sum(ghost),maxGhostCells=max(ghost),maxNeighbourRanks=max(neighbourcounts),
        firstLayerCommunicationVolume=sum(firstghost),threeLayerCommunicationVolume=sum(ghost),
        cellStorageDuplication=float(localcells.sum()/graph.cells),
        faceStorageDuplication=float(sum(faces)/globalfaces),
        componentsPerPart=components_per_part(graph.offsets,graph.neighbours,owners,parts).tolist())
    result['ownershipSHA256']=hashlib.sha256(owners.tobytes()).hexdigest()
    return result


def run_job(graph,algorithm,parts,seed,metis):
    gc.collect()
    sampler=MemorySampler(); sampler.start()
    start=time.perf_counter()
    if algorithm=='graph-growth': raw=graph_growth(graph.offsets,graph.neighbours,parts); objective=None
    elif algorithm in ('rcb','rib','morton'): raw=geometric(graph,parts,algorithm); objective=None
    else: raw,objective=metis.partition(graph,parts,algorithm,seed)
    initial_end=time.perf_counter()
    raw_count=np.bincount(raw,minlength=parts)
    owners,transfers,unsafe=rebalance(graph.offsets,graph.neighbours,raw,quotas(graph.cells,parts))
    end=time.perf_counter()
    memory=sampler.stop()
    result=dict(algorithm=algorithm,parts=parts,seed=seed,partitionSeconds=end-start,
        initialPartitionSeconds=initial_end-start,exactBalanceSeconds=end-initial_end,
        quotaTransfers=int(transfers),unsafeQuotaTransfers=int(unsafe),rawCellImbalance=float(raw_count.max()*parts/raw_count.sum()),
        nativeObjective=objective,memory=memory)
    result.update(metrics(graph,owners,parts))
    return result


def summarize(output):
    records=[]
    for path in sorted((output/'runs').glob('*.json')):
        data=json.loads(path.read_text())
        if data.get('status')=='ok': records.append(data)
    groups={}
    for record in records:
        key=(record['graph'],record['algorithm'],record['parts'])
        groups.setdefault(key,[]).append(record)
    summaries=[]
    for (graph,algorithm,parts),values in groups.items():
        item=dict(graph=graph,algorithm=algorithm,parts=parts,samples=len(values))
        for metric in ('partitionSeconds','cutFaces','totalGhostCells','maxGhostCells',
                       'maxNeighbourRanks','ownedFaceWorkImbalance','localFaceImbalance',
                       'cellStorageDuplication','faceStorageDuplication','firstLayerCommunicationVolume','quotaTransfers'):
            measures=[v[metric] for v in values]
            item[metric]=statistics.median(measures)
            item[metric+'Min']=min(measures); item[metric+'Max']=max(measures)
        item['maxComponentsPerPart']=max(max(v['componentsPerPart']) for v in values)
        summaries.append(item)
    (output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    if summaries:
        with (output/'summary.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(summaries[0])); writer.writeheader(); writer.writerows(summaries)
    ranking=[]
    for algorithm in sorted({s['algorithm'] for s in summaries}):
        ratios=[]; work=[]; components=[]
        for result in summaries:
            if result['algorithm']!=algorithm: continue
            baseline=next((s for s in summaries if s['graph']==result['graph'] and s['parts']==result['parts'] and s['algorithm']=='graph-growth'),None)
            if baseline is None: continue
            ratios.append(result['totalGhostCells']/max(baseline['totalGhostCells'],1))
            work.append(result['ownedFaceWorkImbalanceMax']); components.append(result['maxComponentsPerPart'])
        if ratios:
            ranking.append(dict(algorithm=algorithm,comparisons=len(ratios),
                geometricMeanGhostRatio=float(np.exp(np.log(np.maximum(ratios,1e-12)).mean())),
                worstGhostRatio=max(ratios),maxOwnedFaceWorkImbalance=max(work),maxComponentsPerPart=max(components)))
    ranking.sort(key=lambda r:r['geometricMeanGhostRatio'])
    (output/'ranking.json').write_text(json.dumps(ranking,indent=2)+'\n')
    return summaries,ranking


def plot(output,summaries):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    algorithms=[a for a in ALGORITHMS if any(s['algorithm']==a for s in summaries)]
    ranks=max(s['parts'] for s in summaries)
    graphs=sorted({s['graph'] for s in summaries})
    matrix=np.full((len(graphs),len(algorithms)),np.nan)
    for y,graph in enumerate(graphs):
        base=next((s for s in summaries if s['graph']==graph and s['parts']==ranks and s['algorithm']=='graph-growth'),None)
        if base is None: continue
        for x,algorithm in enumerate(algorithms):
            result=next((s for s in summaries if s['graph']==graph and s['parts']==ranks and s['algorithm']==algorithm),None)
            if result: matrix[y,x]=result['totalGhostCells']/max(base['totalGhostCells'],1)
    fig,ax=plt.subplots(figsize=(13,max(5,len(graphs)*.40)))
    graphic=ax.imshow(matrix,vmin=0,vmax=max(1.5,float(np.nanmax(matrix))),cmap='RdYlGn_r',aspect='auto')
    ax.set_xticks(range(len(algorithms)),algorithms,rotation=30,ha='right'); ax.set_yticks(range(len(graphs)),graphs)
    ax.set_title(f'Total 3-layer ghost cells / graph-growth baseline, {ranks} partitions (lower is better)')
    for y in range(len(graphs)):
        for x in range(len(algorithms)):
            if np.isfinite(matrix[y,x]): ax.text(x,y,f'{matrix[y,x]:.2f}',ha='center',va='center',fontsize=8)
    fig.colorbar(graphic,ax=ax,label='Ratio'); fig.tight_layout()
    fig.savefig(output/'ghost-comparison.png',dpi=170); fig.savefig(output/'ghost-comparison.svg'); plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--families',nargs='+',choices=FAMILIES,default=list(FAMILIES))
    parser.add_argument('--scales',nargs='+',choices=('fixture','small','engineering'),default=['small','engineering'])
    parser.add_argument('--algorithms',nargs='+',choices=ALGORITHMS,default=list(ALGORITHMS))
    parser.add_argument('--ranks',nargs='+',type=int,default=[4,12])
    parser.add_argument('--seeds',nargs='+',type=int,default=[0,7,29])
    parser.add_argument('--metis-header',type=Path,default=Path.home()/'opt/metis-5.1.0/include/metis.h')
    parser.add_argument('--generate',action='store_true'); parser.add_argument('--run',action='store_true')
    parser.add_argument('--analyze',action='store_true'); parser.add_argument('--fixtures',action='store_true')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if any(p<2 or p>63 for p in args.ranks): parser.error('ranks must be between 2 and 63')
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[key]='1'
    args.graph_root.mkdir(parents=True,exist_ok=True); args.output.mkdir(parents=True,exist_ok=True)
    with Path('/tmp/babelsim-solver-benchmark.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.generate:
            metadata=[]
            for scale in args.scales:
                for family in args.families:
                    path=args.graph_root/f'{family}-{scale}.npz'
                    if not path.exists():
                        print(f'GENERATE {path.name}',flush=True)
                        graph=make_graph(family,scale); graph.save(path)
                        if args.fixtures and scale=='fixture': native_fixture(graph,path.with_suffix('.mesh'))
                        del graph; gc.collect()
                    graph=Graph.load(path); metadata.append(dict(path=str(path),sha256=sha(path),**graph.metadata)); del graph
            (args.output/'graphs.json').write_text(json.dumps(metadata,indent=2)+'\n')
        if args.run:
            versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__,
                host=platform.platform(),cpuAffinity=sorted(os.sched_getaffinity(0)),environment={k:os.getenv(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS')})
            sources={name:sha(ROOT/'tools'/name) for name in ('study_mesh_partition.py','partition_model_graphs.py','partition_model_algorithms.py')}
            config=dict(families=args.families,scales=args.scales,algorithms=args.algorithms,ranks=args.ranks,seeds=args.seeds,
                sources=sources,versions=versions,metisHeader=str(args.metis_header),metisHeaderSHA256=sha(args.metis_header),
                sourcesReferences=SOURCES,ghostLayers=3,exactCellQuotas=True)
            configfile=args.output/'study.json'
            if configfile.exists() and (not args.resume or json.loads(configfile.read_text())!=config):
                parser.error('study exists or sources changed; choose a new evidence directory')
            configfile.write_text(json.dumps(config,indent=2)+'\n')
            (args.output/'source').mkdir(exist_ok=True)
            for name in sources: shutil.copy2(ROOT/'tools'/name,args.output/'source'/name)
            (args.output/'runs').mkdir(exist_ok=True)
            print('WARMUP Numba kernels; compilation excluded',flush=True); warm_kernels()
            metis=Metis(args.metis_header)
            for scale in args.scales:
                for family in args.families:
                    graphpath=args.graph_root/f'{family}-{scale}.npz'; graph=Graph.load(graphpath)
                    key=f'{family}-{scale}'; graphsha=sha(graphpath)
                    jobs=[(a,p,seed) for a in args.algorithms for p in args.ranks
                          for seed in (args.seeds if a.startswith('metis-') else [0])]
                    random.Random(20261007).shuffle(jobs)
                    for algorithm,parts,seed in jobs:
                        name=f'{key}__{algorithm}__p{parts}__seed{seed}'
                        recordpath=args.output/'runs'/f'{name}.json'
                        if args.resume and recordpath.exists(): continue
                        print(f'RUN {name} cells={graph.cells}',flush=True)
                        start=time.time()
                        try:
                            result=run_job(graph,algorithm,parts,seed,metis); result['status']='ok'
                        except Exception as error:
                            result=dict(status='failed',error=f'{type(error).__name__}: {error}',algorithm=algorithm,parts=parts,seed=seed)
                        result.update(graph=key,graphSHA256=graphsha,graphMetadata=graph.metadata,unixStarted=start,unixEnded=time.time())
                        recordpath.write_text(json.dumps(result,indent=2)+'\n')
                        if result['status']=='ok':
                            print(f"DONE partition={result['partitionSeconds']:.4f}s cut={result['cutFaces']} ghosts={result['totalGhostCells']} workImbalance={result['ownedFaceWorkImbalance']:.4f} components={max(result['componentsPerPart'])}",flush=True)
                        else: print(f"FAILED {result['error']}",flush=True)
                    del graph; gc.collect(); summarize(args.output)
        if args.analyze or args.run:
            summaries,ranking=summarize(args.output)
            if summaries: plot(args.output,summaries)
            print(json.dumps(ranking,indent=2),flush=True)
    return 0

if __name__=='__main__': raise SystemExit(main())
