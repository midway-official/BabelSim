"""Python partition prototypes and compiled graph kernels; no solver access."""
from __future__ import annotations

import ctypes as ct
import ctypes.util
from heapq import heapify, heappop, heappush
from pathlib import Path
import re
import numpy as np
from numba import njit

ALGORITHMS = ('graph-growth', 'rcb', 'rib', 'morton', 'morton-isotropic', 'metis-cut',
              'metis-cut-work', 'metis-volume-work', 'metis-volume-contig')


def quotas(cells, parts):
    return np.full(parts,cells//parts,dtype=np.int32)+(np.arange(parts)<cells%parts)


@njit(cache=True,nogil=True)
def graph_growth(offsets, neighbours, parts):
    """Faithful seed/round-robin model of the existing C++ graph partitioner."""
    cells=len(offsets)-1
    capacity=np.full(parts,cells//parts,np.int32)
    for part in range(cells%parts): capacity[part]+=1
    filled=np.zeros(parts,np.int32)
    nearest=np.full(cells,2147483647,np.int32)
    seeded=np.zeros(cells,np.uint8)
    seeds=np.empty(parts,np.int32)
    queue=np.empty(cells,np.int32)
    for part in range(parts):
        seed=-1
        for cell in range(cells):
            if not seeded[cell] and (seed==-1 or nearest[cell]>nearest[seed]): seed=cell
        seeds[part]=seed; seeded[seed]=1; nearest[seed]=0
        head,tail=0,1; queue[0]=seed
        while head<tail:
            cell=queue[head]; head+=1; distance=nearest[cell]+1
            for edge in range(offsets[cell],offsets[cell+1]):
                neighbour=neighbours[edge]
                if nearest[neighbour]>distance:
                    nearest[neighbour]=distance; queue[tail]=neighbour; tail+=1
    owners=np.full(cells,-1,np.int32)
    following=np.full(cells,-1,np.int32)
    heads=seeds.copy(); tails=seeds.copy()
    for part in range(parts): owners[seeds[part]]=part; filled[part]=1
    remaining=cells-parts
    while remaining:
        grew=False
        for part in range(parts):
            if filled[part]>=capacity[part]: continue
            while heads[part]!=-1:
                cell=heads[part]; selected=-1
                for edge in range(offsets[cell],offsets[cell+1]):
                    neighbour=neighbours[edge]
                    if owners[neighbour]==-1: selected=neighbour; break
                if selected!=-1:
                    owners[selected]=part; filled[part]+=1; remaining-=1
                    following[tails[part]]=selected; tails[part]=selected
                    grew=True; break
                heads[part]=following[cell]
                if heads[part]==-1: tails[part]=-1
        if grew: continue
        part=0
        while filled[part]>=capacity[part]: part+=1
        seed=0
        while owners[seed]!=-1: seed+=1
        owners[seed]=part; filled[part]+=1; remaining-=1
        if heads[part]==-1: heads[part]=seed
        else: following[tails[part]]=seed
        tails[part]=seed
    return owners


@njit(cache=True,nogil=True)
def transfer_cost(offsets, neighbours, owners, cell, target):
    cost=0; source=owners[cell]
    for edge in range(offsets[cell],offsets[cell+1]):
        part=owners[neighbours[edge]]
        cost+=int(part==source)-int(part==target)
    return cost


@njit(cache=True,nogil=True)
def removable(offsets,neighbours,owners,cell,visited,goals,queue,sequence,budget):
    """A bounded witness proves safe removal; full checking is used for final fallback moves."""
    source=owners[cell]; sequence[0]+=1; stamp=sequence[0]; remaining=0; seed=-1
    for edge in range(offsets[cell],offsets[cell+1]):
        neighbour=neighbours[edge]
        if owners[neighbour]==source:
            goals[neighbour]=stamp; seed=neighbour; remaining+=1
    if remaining<=1: return True
    visited[cell]=stamp; visited[seed]=stamp; remaining-=1
    queue[0]=seed; head,tail=0,1
    while head<tail and head<budget:
        current=queue[head]; head+=1
        for edge in range(offsets[current],offsets[current+1]):
            neighbour=neighbours[edge]
            if owners[neighbour]!=source or visited[neighbour]==stamp: continue
            visited[neighbour]=stamp
            if goals[neighbour]==stamp:
                remaining-=1
                if remaining==0: return True
            queue[tail]=neighbour; tail+=1
    return False


@njit(cache=True,nogil=True)
def balance_path(adjacency,counts,capacity):
    parts=len(capacity); parent=np.full(parts,-1,np.int32)
    queue=np.empty(parts,np.int32); head=tail=0
    for part in range(parts):
        if counts[part]>capacity[part]: parent[part]=part;queue[tail]=part;tail+=1
    deficit=-1
    while head<tail and deficit==-1:
        source=queue[head];head+=1
        for target in range(parts):
            if not adjacency[source,target] or parent[target]!=-1: continue
            parent[target]=source
            if counts[target]<capacity[target]: deficit=target;break
            queue[tail]=target;tail+=1
    path=np.empty(parts,np.int32);length=0
    if deficit!=-1:
        while parent[deficit]!=deficit:
            path[length]=deficit;length+=1;deficit=parent[deficit]
        path[length]=deficit;length+=1
        path[:length]=path[:length][::-1]
    return path,length


@njit(cache=True,nogil=True)
def push_candidates(heap,offsets,neighbours,owners,counts,capacity,versions,cell):
    source=owners[cell]
    if counts[source]<=capacity[source]:return
    targets=np.uint64(0)
    for edge in range(offsets[cell],offsets[cell+1]):
        target=owners[neighbours[edge]]
        bit=np.uint64(1)<<np.uint64(target)
        if target==source or counts[target]>=capacity[target] or targets&bit:continue
        targets|=bit
        cost=transfer_cost(offsets,neighbours,owners,cell,target)
        heappush(heap,(np.int64(cost),np.int64(cell),np.int64(target),versions[cell]))


@njit(cache=True,nogil=True)
def rebalance(offsets, neighbours, initial, capacity):
    owners=initial.copy();parts=len(capacity);cells=len(owners)
    counts=np.zeros(parts,np.int32)
    for owner in owners:counts[owner]+=1
    if np.all(counts==capacity):return owners,0,0
    if len(neighbours)==0:
        # No interfaces exist; linear-time quota repair avoids scanning an
        # ever-growing list of disconnected components for every transfer.
        target=0;transfers=0
        for cell in range(cells):
            source=owners[cell]
            if counts[source]<=capacity[source]:continue
            while target<parts and counts[target]>=capacity[target]:target+=1
            assert target<parts
            owners[cell]=target;counts[source]-=1;counts[target]+=1;transfers+=1
        return owners,transfers,0
    boundary=np.empty(cells,np.int32);recorded=np.zeros(cells,np.uint8);length=0
    visited=np.zeros(cells,np.int64);goals=np.zeros(cells,np.int64)
    queue=np.empty(cells,np.int32);sequence=np.zeros(1,np.int64)
    mobility=np.zeros(cells,np.uint8);mobility_epoch=np.full(cells,-1,np.int64)
    versions=np.zeros(cells,np.int64);state=np.int64(0)
    for cell in range(cells):
        for edge in range(offsets[cell],offsets[cell+1]):
            if owners[cell]!=owners[neighbours[edge]]:
                recorded[cell]=1;boundary[length]=cell;length+=1;break
    heap=[(np.int64(0),np.int64(0),np.int64(0),np.int64(0))];heap.pop()
    for at in range(length):push_candidates(heap,offsets,neighbours,owners,counts,capacity,versions,boundary[at])
    transfers=unsafe=0
    while np.any(counts!=capacity):
        bestcell=-1;besttarget=-1
        while heap:
            cost,cell,target,version=heappop(heap)
            source=owners[cell]
            if version!=versions[cell] or counts[source]<=capacity[source] or counts[target]>=capacity[target]:continue
            if mobility_epoch[cell]!=state:
                mobility[cell]=removable(offsets,neighbours,owners,cell,visited,goals,queue,sequence,64)
                mobility_epoch[cell]=state
            if not mobility[cell]:continue
            bestcell,besttarget=cell,target;break
        path=np.empty(parts,np.int32);pathlen=0;routed=bestcell==-1
        if bestcell!=-1:
            path[0]=besttarget;pathlen=1
        else:
            adjacency=np.zeros((parts,parts),np.uint8);safe_adjacency=adjacency.copy()
            for at in range(length):
                cell=boundary[at];source=owners[cell]
                if mobility_epoch[cell]!=state:
                    mobility[cell]=removable(offsets,neighbours,owners,cell,visited,goals,queue,sequence,64)
                    mobility_epoch[cell]=state
                for edge in range(offsets[cell],offsets[cell+1]):
                    target=owners[neighbours[edge]]
                    if source!=target:
                        adjacency[source,target]=1
                        if mobility[cell]:safe_adjacency[source,target]=1
            path,pathlen=balance_path(safe_adjacency,counts,capacity)
            if not pathlen:path,pathlen=balance_path(adjacency,counts,capacity)
            if not pathlen:
                target=0
                while counts[target]>=capacity[target]:target+=1
                bestcost=2147483647
                for attempt in range(2):
                    for cell in range(cells):
                        if counts[owners[cell]]<=capacity[owners[cell]]:continue
                        cost=transfer_cost(offsets,neighbours,owners,cell,target)
                        better=cost<bestcost or (cost==bestcost and (bestcell==-1 or cell<bestcell))
                        if not better:continue
                        if attempt==0 and mobility_epoch[cell]!=state:
                            mobility[cell]=removable(offsets,neighbours,owners,cell,visited,goals,queue,sequence,64)
                            mobility_epoch[cell]=state
                        if attempt==1 or mobility[cell]:bestcost,bestcell=cost,cell
                    if bestcell!=-1:break
                path[0]=target;pathlen=1
        begin=0 if bestcell!=-1 else 1
        for step in range(begin,pathlen):
            target=path[step]
            if bestcell==-1 or step>begin:
                source=path[step-1];bestcost=2147483647;bestcell=-1
                for attempt in range(2):
                    for at in range(length):
                        cell=boundary[at]
                        if owners[cell]!=source:continue
                        adjacent=False
                        for edge in range(offsets[cell],offsets[cell+1]):adjacent|=owners[neighbours[edge]]==target
                        if not adjacent:continue
                        cost=transfer_cost(offsets,neighbours,owners,cell,target)
                        better=cost<bestcost or (cost==bestcost and (bestcell==-1 or cell<bestcell))
                        if not better:continue
                        if attempt==0 and mobility_epoch[cell]!=state:
                            mobility[cell]=removable(offsets,neighbours,owners,cell,visited,goals,queue,sequence,64)
                            mobility_epoch[cell]=state
                        if attempt==1 or mobility[cell]:bestcost,bestcell=cost,cell
                    if bestcell!=-1:break
            assert bestcell!=-1
            if mobility_epoch[bestcell]!=state or not mobility[bestcell]:
                unsafe+=int(not removable(offsets,neighbours,owners,bestcell,visited,goals,queue,sequence,cells))
            counts[owners[bestcell]]-=1;counts[target]+=1;owners[bestcell]=target;transfers+=1;state+=1
            versions[bestcell]+=1
            if not recorded[bestcell]:recorded[bestcell]=1;boundary[length]=bestcell;length+=1
            if not routed:push_candidates(heap,offsets,neighbours,owners,counts,capacity,versions,bestcell)
            for edge in range(offsets[bestcell],offsets[bestcell+1]):
                cell=neighbours[edge];versions[cell]+=1
                if not recorded[cell]:recorded[cell]=1;boundary[length]=cell;length+=1
                if not routed:push_candidates(heap,offsets,neighbours,owners,counts,capacity,versions,cell)
        if routed:
            # A routing step can make a formerly balanced part a temporary
            # surplus. Rebuild once after the route, rather than per transfer.
            heap.clear()
            for at in range(length):push_candidates(heap,offsets,neighbours,owners,counts,capacity,versions,boundary[at])
    return owners,transfers,unsafe


@njit(cache=True,nogil=True)
def halo_masks(offsets,neighbours,owners,layers):
    masks=np.empty(len(owners),np.uint64)
    for cell in range(len(owners)): masks[cell]=np.uint64(1)<<np.uint64(owners[cell])
    first=masks.copy()
    for layer in range(layers):
        following=masks.copy()
        for cell in range(len(owners)):
            for edge in range(offsets[cell],offsets[cell+1]): following[cell]|=masks[neighbours[edge]]
        masks=following
        if layer==0: first=masks.copy()
    return masks,first


@njit(cache=True,nogil=True)
def components_per_part(offsets,neighbours,owners,parts):
    visited=np.zeros(len(owners),np.uint8); queue=np.empty(len(owners),np.int32)
    counts=np.zeros(parts,np.int32)
    for seed in range(len(owners)):
        if visited[seed]: continue
        part=owners[seed]; counts[part]+=1; visited[seed]=1; queue[0]=seed; head,tail=0,1
        while head<tail:
            cell=queue[head]; head+=1
            for edge in range(offsets[cell],offsets[cell+1]):
                neighbour=neighbours[edge]
                if not visited[neighbour] and owners[neighbour]==part:
                    visited[neighbour]=1; queue[tail]=neighbour; tail+=1
    return counts


def geometric(graph,parts,algorithm):
    capacity=quotas(graph.cells,parts); owners=np.empty(graph.cells,np.int32)
    if algorithm.startswith('morton'):
        xyz=graph.centres
        extent=np.ptp(xyz,axis=0); extent[extent==0]=1
        if algorithm=='morton-isotropic': extent[:]=extent.max()
        integer=np.minimum((xyz-xyz.min(axis=0))/extent*((1<<21)-1),(1<<21)-1).astype(np.uint64)
        codes=np.zeros(graph.cells,np.uint64)
        for bit in range(21):
            for axis in range(3): codes|=((integer[:,axis]>>np.uint64(bit))&np.uint64(1))<<np.uint64(3*bit+axis)
        order=np.argsort(codes,kind='stable'); offset=0
        for part,count in enumerate(capacity): owners[order[offset:offset+count]]=part; offset+=count
        return owners
    def split(indices,begin,end):
        if end-begin==1: owners[indices]=begin; return
        mid=(begin+end)//2; count=int(capacity[begin:mid].sum())
        points=graph.centres[indices]
        if algorithm=='rcb': projection=points[:,np.argmax(np.ptp(points,axis=0))]
        else:
            shifted=points-points.mean(axis=0)
            values,vectors=np.linalg.eigh(shifted.T@shifted)
            direction=vectors[:,-1]
            if direction[np.argmax(np.abs(direction))]<0: direction=-direction
            projection=shifted@direction
        order=np.lexsort((indices,projection)); ranked=indices[order]
        split(ranked[:count],begin,mid); split(ranked[count:],mid,end)
    split(np.arange(graph.cells,dtype=np.int32),0,parts)
    return owners


class Metis:
    """ctypes binding derived from the installed, matching public metis.h."""
    def __init__(self,header):
        self.header=Path(header); source=self.header.read_text()
        width=int(re.search(r'#define IDXTYPEWIDTH\s+(\d+)',source)[1])
        realwidth=int(re.search(r'#define REALTYPEWIDTH\s+(\d+)',source)[1])
        self.index,self.dtype=(ct.c_int32,np.int32) if width==32 else (ct.c_int64,np.int64)
        self.real,self.real_dtype=(ct.c_float,np.float32) if realwidth==32 else (ct.c_double,np.float64)
        self.constants={}
        for enum in re.findall(r'typedef\s+enum\s*\{(.*?)\}\s*\w+\s*;',source,re.S):
            enum=re.sub(r'/\*.*?\*/|//[^\n]*','',enum,flags=re.S)
            value=0
            for token in enum.split(','):
                token=token.strip()
                if not token: continue
                name,sep,explicit=token.partition('=')
                if sep: value=int(explicit.strip(),0)
                self.constants[name.strip()]=value; value+=1
        self.noptions=int(re.search(r'#define METIS_NOPTIONS\s+(\d+)',source)[1])
        self.lib=ct.CDLL(ctypes.util.find_library('metis'))
        ip=ct.POINTER(self.index); rp=ct.POINTER(self.real)
        self.lib.METIS_SetDefaultOptions.argtypes=[ip]
        self.lib.METIS_PartGraphKway.argtypes=[ip,ip,ip,ip,ip,ip,ip,ip,rp,rp,ip,ip,ip]
        self.lib.METIS_PartGraphKway.restype=ct.c_int

    def partition(self,graph,parts,algorithm,seed):
        weighted=algorithm!='metis-cut'
        ncon=2 if weighted else 1
        cells,con,nparts,objective=self.index(graph.cells),self.index(ncon),self.index(parts),self.index(0)
        offsets=graph.offsets.astype(self.dtype,copy=True); neighbours=graph.neighbours.astype(self.dtype,copy=True)
        weights=np.column_stack((np.ones(graph.cells,np.int32),graph.face_counts)).astype(self.dtype).ravel() if weighted else None
        capacity=quotas(graph.cells,parts)
        targets=np.repeat(capacity.astype(self.real_dtype)/graph.cells,ncon).astype(self.real_dtype)
        imbalance=np.array([1.001,1.05][:ncon],dtype=self.real_dtype)
        options=np.empty(self.noptions,dtype=self.dtype)
        ip=ct.POINTER(self.index); rp=ct.POINTER(self.real)
        self.lib.METIS_SetDefaultOptions(options.ctypes.data_as(ip))
        for name,value in {'SEED':seed,'NCUTS':2,'MINCONN':1,'NUMBERING':0,
                           'OBJTYPE':self.constants['METIS_OBJTYPE_VOL' if 'volume' in algorithm else 'METIS_OBJTYPE_CUT'],
                           'CONTIG':int('contig' in algorithm and graph.metadata['components']==1)}.items():
            options[self.constants['METIS_OPTION_'+name]]=value
        owners=np.empty(graph.cells,dtype=self.dtype)
        status=self.lib.METIS_PartGraphKway(ct.byref(cells),ct.byref(con),offsets.ctypes.data_as(ip),
            neighbours.ctypes.data_as(ip),None if weights is None else weights.ctypes.data_as(ip),None,None,
            ct.byref(nparts),targets.ctypes.data_as(rp),imbalance.ctypes.data_as(rp),options.ctypes.data_as(ip),
            ct.byref(objective),owners.ctypes.data_as(ip))
        if status!=1: raise RuntimeError(f'METIS returned {status}')
        return owners.astype(np.int32),int(objective.value)


def warm_kernels():
    offsets=np.array([0,1,3,4],np.int32); neighbours=np.array([1,0,2,1],np.int32)
    labels=graph_growth(offsets,neighbours,2)
    rebalance(offsets,neighbours,labels,np.array([2,1],np.int32))
    halo_masks(offsets,neighbours,labels,3); components_per_part(offsets,neighbours,labels,2)
