#!/usr/bin/env python3
"""Deterministic hybrid unstructured mesh: normal quad layers + constrained
Delaunay transition + Cartesian background, extruded by one spanwise cell.
Dependencies: numpy, triangle (see requirements-mesh.txt).
"""
from pathlib import Path
from collections import defaultdict
import json
import math
import numpy as np
import triangle
from mesh_io import write_volume_mesh, signed_area

ROOT=Path(__file__).resolve().parent
X_MIN,X_MAX,Y_MIN,Y_MAX=-8.,4.,-5.,5.
SPAN_THICKNESS=.02
FIRST_HEIGHT=3e-5
GROWTH=1.12
LAYERS=48
REYNOLDS=1e6

def graded(a,b,step=.03,ratio=1.055):
    n=math.ceil(math.log(1+abs(b-a)*(ratio-1)/step)/math.log(ratio))
    w=ratio**np.arange(n); w*=abs(b-a)/sum(w)
    return np.r_[a,a+np.sign(b-a)*np.cumsum(w)]

def axes():
    # Fine wake occupies y=[-1,1], x=[-6,1]. Interfaces have matching nodes.
    x=np.r_[graded(-6,X_MIN)[::-1][:-1],np.linspace(-6,-1,201)[:-1],
            np.linspace(-1,1,101)[:-1],graded(1,X_MAX)]
    y=np.r_[graded(-1,Y_MIN)[::-1][:-1],np.linspace(-1,1,201)[:-1],graded(1,Y_MAX)]
    return x,y

def wall_points():
    # Standard NACA0012 except a smooth cap confined to the final 0.15% chord.
    # Avoid a singular sharp TE and retain all-hex normal boundary layers.
    t=np.linspace(0,math.acos(1-2*.9985),361)
    x=.5*(1-np.cos(t))
    yt=lambda x:.6*(.2969*np.sqrt(x)-.126*x-.3516*x*x+.2843*x**3-.1015*x**4)
    upper=np.column_stack((x,yt(x)))
    a=upper[-1]; b=a+np.array([.0008,-.0008*.1403]); c=np.array([1.,.0007]); d=np.array([1.,0.])
    cap=[]
    for s in np.linspace(0,1,13)[1:]: cap.append((1-s)**3*a+3*(1-s)**2*s*b+3*(1-s)*s*s*c+s**3*d)
    upper=np.vstack((upper,cap))
    wall=np.vstack((upper,upper[-2:0:-1]*[1,-1]))
    # CCW solid polygon; LE at (+.2415,-.0647), TE at (-.7244,+.1941).
    if signed_area(wall.tolist())<0: wall=wall[::-1]
    theta=math.radians(165)
    rot=np.array([[math.cos(theta),-math.sin(theta)],[math.sin(theta),math.cos(theta)]])
    return (wall-[.25,0])@rot.T

def recombine(vertices, triangles):
    """Merge well-shaped adjacent triangles; keep difficult pairs as prisms."""
    adjacency=defaultdict(list)
    for i,t in enumerate(triangles):
        for a,b in zip(t,np.roll(t,-1)): adjacency[tuple(sorted((int(a),int(b))))].append(i)
    candidates=[]
    for pair in adjacency.values():
        if len(pair)!=2: continue
        i,j=pair; ring=list(set(triangles[i])|set(triangles[j]))
        if len(ring)!=4: continue
        p=vertices[ring]; center=p.mean(axis=0)
        ring=[ring[k] for k in np.argsort(np.arctan2(p[:,1]-center[1],p[:,0]-center[0]))]
        p=vertices[ring]; edges=np.roll(p,-1,axis=0)-p
        turns=edges[:,0]*np.roll(edges,-1,axis=0)[:,1]-edges[:,1]*np.roll(edges,-1,axis=0)[:,0]
        if min(turns)<=0: continue
        lengths=np.linalg.norm(edges,axis=1)
        angles=np.degrees(np.arccos(np.clip(np.sum(-np.roll(edges,1,axis=0)*edges,axis=1)/(np.roll(lengths,1)*lengths),-1,1)))
        score=max(abs(angles-90))
        if score<=35 and max(lengths)/min(lengths)<3: candidates.append((score,i,j,ring))
    used=set(); result=[]
    for _,i,j,ring in sorted(candidates):
        if i in used or j in used: continue
        used.update((i,j)); result.append(ring)
    result.extend(t for i,t in enumerate(triangles) if i not in used)
    return result

def build_2d_mesh():
    wall=wall_points(); n=len(wall)
    edges=np.roll(wall,-1,axis=0)-wall
    length=np.linalg.norm(edges,axis=1)
    normal=np.column_stack((edges[:,1],-edges[:,0]))/length[:,None]
    prev=np.roll(normal,1,axis=0)
    bisector=(normal+prev)/(1+np.sum(normal*prev,axis=1))[:,None]
    distances=np.r_[0,np.cumsum(FIRST_HEIGHT*GROWTH**np.arange(LAYERS))]
    rows=wall[None,:,:]+distances[:,None,None]*bisector[None,:,:]
    points=[]; lookup={}; cells=[]; boundary={}
    def vertex(p):
        key=tuple(np.round(p,12))
        if key not in lookup: lookup[key]=len(points); points.append(tuple(map(float,p)))
        return lookup[key]
    ids=[[vertex(p) for p in row] for row in rows]
    for j in range(LAYERS):
        for i in range(n):
            k=(i+1)%n
            cells.append((ids[j][k],ids[j][i],ids[j+1][i],ids[j+1][k]))
    for i in range(n): boundary[tuple(sorted((ids[0][i],ids[0][(i+1)%n])))]='airfoil'
    bl_end=len(cells)
    x,y=axes()
    xi=x[(x>=-1-1e-9)&(x<=1+1e-9)]; yi=y[(y>=-1-1e-9)&(y<=1+1e-9)]
    outer=np.array([(v,-1) for v in xi[:-1]]+[(1,v) for v in yi[:-1]]+
                   [(v,1) for v in xi[:0:-1]]+[(-1,v) for v in yi[:0:-1]])
    inner=rows[-1]; m=len(outer)
    vertices=np.vstack((outer,inner))
    seg=[(i,(i+1)%m) for i in range(m)]+[(m+i,m+(i+1)%n) for i in range(n)]
    mesh=triangle.triangulate({'vertices':vertices,'segments':np.array(seg),'holes':np.array([[0.,0.]])},'pYq28a0.0006Q')
    mapping=[vertex(p) for p in mesh['vertices']]
    for tri in recombine(mesh['vertices'],mesh['triangles']):
        c=tuple(mapping[i] for i in tri)
        if signed_area([points[i] for i in c])<0: c=c[::-1]
        cells.append(c)
    transition_end=len(cells)
    for j in range(len(y)-1):
        for i in range(len(x)-1):
            xm=.5*(x[i]+x[i+1]); ym=.5*(y[j]+y[j+1])
            if -1<xm<1 and -1<ym<1: continue
            cells.append(tuple(vertex(p) for p in [(x[i],y[j]),(x[i+1],y[j]),(x[i+1],y[j+1]),(x[i],y[j+1])]))
    owners=defaultdict(list)
    for ci,c in enumerate(cells):
        assert signed_area([points[v] for v in c])>0, ('non-positive area',ci)
        for a,b in zip(c,c[1:]+c[:1]): owners[tuple(sorted((a,b)))].append(ci)
    for e,own in owners.items():
        assert len(own)<=2,('non-manifold',e)
        if len(own)==2:
            assert e not in boundary
            continue
        if e in boundary: continue
        a,b=np.array([points[v] for v in e])
        if abs(a[0]-X_MIN)<1e-9 and abs(b[0]-X_MIN)<1e-9: patch='outlet'
        elif abs(a[0]-X_MAX)<1e-9 and abs(b[0]-X_MAX)<1e-9: patch='inlet'
        elif abs(abs(a[1])-5)<1e-9 and abs(a[1]-b[1])<1e-9: patch='farfield'
        else: raise ValueError(('unmatched edge',e,a,b))
        boundary[e]=patch
    return points,cells,boundary,{'region_cell_ranges':{'wall_boundary_layer':[0,bl_end],
        'unstructured_transition':[bl_end,transition_end],'cartesian_background':[transition_end,len(cells)]},
        'wall_surface_faces':n,'layers':LAYERS,'first_layer_height':FIRST_HEIGHT,'growth':GROWTH,
        'boundary_layer_thickness':float(distances[-1]),'domain':[X_MIN,X_MAX,Y_MIN,Y_MAX],
        'angle_of_attack_deg':15,'freestream_velocity':[-1,0,0],'reynolds_number':REYNOLDS,
        'geometry_note':'NACA0012 with smooth trailing cap over last 0.0015 chord; chord=1',
        'first_cell_y_plus_flat_plate_estimate':math.sqrt(.5*.026/REYNOLDS**(1/7))*FIRST_HEIGHT*.5*REYNOLDS}

def main():
    points,cells,boundary,metadata=build_2d_mesh()
    _,_,quality=write_volume_mesh(ROOT/'mesh/naca0012.mesh',points,cells,boundary,metadata['region_cell_ranges'])
    report=metadata|quality
    report['hexahedral_cells']=len(cells)-quality['non_hexahedral_cells']
    report['hexahedral_fraction']=report['hexahedral_cells']/len(cells)
    report['quality_gate']={'cell_count_100k_200k':100000<=len(cells)<=200000,
        'positive_volumes':quality['negative_or_zero_volume_cells']==0,
        'convex_cells':quality['concave_cells']==0,
        'max_nonorthogonality_below_70':quality['internal_nonorthogonality_deg_max']<70,
        'min_corner_above_15':quality['minimum_cell_corner_angle_deg']>15}
    (ROOT/'mesh/mesh_quality.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if not all(report['quality_gate'].values()): raise RuntimeError('mesh quality gates failed')
if __name__=='__main__': main()
