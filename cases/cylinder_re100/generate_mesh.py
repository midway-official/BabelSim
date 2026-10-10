#!/usr/bin/env python3
"""Generate a conformal triangular-prism cylinder mesh; Gmsh handles outer CDT.
Python constructs all annular layers, merges the interface, writes explicit
polyhedral faces, audits topology/geometry, and creates a solenoidal initial field.
Requires numpy, matplotlib, and the gmsh executable (MSH 2 ASCII output).
"""
from pathlib import Path
import argparse, collections, json, math, subprocess
import numpy as np

CASE = Path(__file__).resolve().parent
NTHETA, NLAYERS, FIRST, GROWTH = 160, 12, .015, 1.10
R, THICKNESS = .5, .1
XMIN, XMAX, YMIN, YMAX = -10., 25., -10., 10.
PATCHES = [('cylinder','wall'),('inlet','inlet'),('outlet','outlet'),
           ('farfield','symmetry'),('front','symmetry'),('back','symmetry')]

def cross(a,b): return a[...,0]*b[...,1]-a[...,1]*b[...,0]

def generate(wake_size):
    radii = np.r_[R, R+np.cumsum(FIRST*GROWTH**np.arange(NLAYERS))]
    angles = np.arange(NTHETA)*2*np.pi/NTHETA
    points = [(r*np.cos(t),r*np.sin(t)) for r in radii for t in angles]
    triangles=[]
    for j in range(NLAYERS):
        for i in range(NTHETA):
            a=j*NTHETA+i; b=j*NTHETA+(i+1)%NTHETA
            c=b+NTHETA; d=a+NTHETA
            triangles.extend([(a,d,c),(a,c,b)] if (i+j)%2 == 0 else [(a,d,b),(b,d,c)])
    layer_cells=len(triangles)
    outer=points[-NTHETA:]
    geo=['Mesh.MshFileVersion = 2.2;', 'Mesh.Algorithm = 6;', 'Mesh.Smoothing = 20;',
         'General.NumThreads = 1;', 'Mesh.MeshSizeFromPoints = 0;',
         'Mesh.MeshSizeFromCurvature = 0;', 'Mesh.MeshSizeExtendFromBoundary = 0;']
    for i,(x,y) in enumerate(outer,1): geo.append(f'Point({i}) = {{{x:.16g},{y:.16g},0,1}};')
    for i in range(1,NTHETA+1):geo.append(f'Line({i}) = {{{i},{i%NTHETA+1}}};')
    geo.append(f'Transfinite Curve {{1:{NTHETA}}} = 2;')
    for i,(x,y) in enumerate([(XMIN,YMIN),(XMAX,YMIN),(XMAX,YMAX),(XMIN,YMAX)],NTHETA+1):
        geo.append(f'Point({i}) = {{{x},{y},0,1}};')
    for k in range(4):geo.append(f'Line({NTHETA+1+k}) = {{{NTHETA+1+k},{NTHETA+1+(k+1)%4}}};')
    geo += [f'Curve Loop(1) = {{1:{NTHETA}}};', f'Curve Loop(2) = {{{NTHETA+1}:{NTHETA+4}}};',
            'Plane Surface(1) = {2,1};', 'Field[1] = MathEval;',
            f'Field[1].F = "0.031 + 0.19*Max(0,Sqrt(x*x+y*y)-{radii[-1]:.16g})";',
            'Field[2] = Box;', f'Field[2].VIn = {wake_size};', 'Field[2].VOut = 1.2;',
            'Field[2].XMin = 0;', 'Field[2].XMax = 15;',
            'Field[2].YMin = -2;', 'Field[2].YMax = 2;',
            'Field[2].Thickness = 2.5;', 'Field[3] = Min;',
            'Field[3].FieldsList = {1,2};', 'Background Field = 3;',
            'Mesh.MeshSizeMax = 1.2;', 'Physical Surface(1) = {1};']
    meshdir=CASE/'mesh';meshdir.mkdir(exist_ok=True)
    (meshdir/'outer.geo').write_text('\n'.join(geo)+'\n')
    with (meshdir/'gmsh.log').open('w') as log:
        subprocess.run(['gmsh',str(meshdir/'outer.geo'),'-2','-format','msh2','-o',str(meshdir/'outer.msh')],stdout=log,stderr=subprocess.STDOUT,check=True)
    lines=(meshdir/'outer.msh').read_text().splitlines()
    lookup={tuple(np.round(p,11)):i for i,p in enumerate(points)}; mapping={}
    start=lines.index('$Nodes'); n=int(lines[start+1])
    for line in lines[start+2:start+2+n]:
        tag,x,y,_=line.split();p=(float(x),float(y));key=tuple(np.round(p,11))
        if key not in lookup:lookup[key]=len(points);points.append(p)
        mapping[int(tag)]=lookup[key]
    start=lines.index('$Elements');n=int(lines[start+1])
    for line in lines[start+2:start+2+n]:
        s=list(map(int,line.split()))
        if s[1]==2:triangles.append(tuple(mapping[k] for k in s[3+s[2]:]))
    xy=np.asarray(points);tri=np.asarray(triangles,dtype=int)
    a=xy[tri];signed=cross(a[:,1]-a[:,0],a[:,2]-a[:,0]);bad=signed<0
    tri[bad]=tri[bad][:,[0,2,1]]
    return xy,tri,layer_cells,radii

def audit_and_write(xy,tri,layer_cells,radii,wake_size):
    meshdir=CASE/'mesh'; n=len(xy);nt=len(tri)
    p=xy[tri]; centres=p.mean(axis=1);areas=cross(p[:,1]-p[:,0],p[:,2]-p[:,0])/2
    edge_lengths=np.linalg.norm(np.roll(p,-1,axis=1)-p,axis=2)
    corner_angles=[]
    for i in range(3):
        u=p[:,(i+1)%3]-p[:,i];v=p[:,(i+2)%3]-p[:,i]
        corner_angles.append(np.degrees(np.arccos(np.clip(np.sum(u*v,axis=1)/np.linalg.norm(u,axis=1)/np.linalg.norm(v,axis=1),-1,1))))
    edge_map={}
    for c,t in enumerate(tri):
        for i in range(3):
            a,b=int(t[i]),int(t[(i+1)%3]);key=tuple(sorted((a,b)))
            if key in edge_map:
                assert edge_map[key][3]==-1,'nonmanifold edge'
                edge_map[key][3]=c
            else:edge_map[key]=[a,b,c,-1]
    records=[];patch_counts=collections.Counter();nonorth=[];skew=[];boundary_nonorth=[]
    edge_cells=[];wall=[];interface=0
    for a,b,o,nb in edge_map.values():
        centre=(xy[a]+xy[b])/2;e=xy[b]-xy[a];normal=np.array([e[1],-e[0]])/np.linalg.norm(e)
        patch=-1
        if nb<0:
            if a<NTHETA and b<NTHETA:patch=0;wall.append([a,b,o])
            elif abs(centre[0]-XMIN)<1e-9:patch=1
            elif abs(centre[0]-XMAX)<1e-9:patch=2
            elif abs(abs(centre[1])-YMAX)<1e-9:patch=3
            else:raise ValueError(f'unassigned boundary edge {a,b,centre}')
            patch_counts[PATCHES[patch][0]]+=1
            d=centre-centres[o]
            boundary_nonorth.append(float(np.degrees(np.arccos(np.clip(np.dot(d,normal)/np.linalg.norm(d),-1,1)))))
        else:
            d=centres[nb]-centres[o];cos=np.dot(d,normal)/np.linalg.norm(d)
            nonorth.append(float(np.degrees(np.arccos(np.clip(cos,-1,1)))))
            f=np.dot(centre-centres[o],normal)/np.dot(d,normal)
            skew.append(float(np.linalg.norm(centre-centres[o]-f*d)/np.linalg.norm(d)))
            edge_cells.append([o,nb,a,b])
            interface+=int((o<layer_cells)!=(nb<layer_cells))
        records.append(([a,b,b+n,a+n],o,nb,patch))
    for c,(a,b,d) in enumerate(tri):
        records.append(([int(d),int(b),int(a)],c,-1,4))
        records.append(([int(a+n),int(b+n),int(d+n)],c,-1,5))
    expected_area=(XMAX-XMIN)*(YMAX-YMIN)-.5*NTHETA*R*R*np.sin(2*np.pi/NTHETA)
    quality=dict(cells=nt,triangular_prisms=nt,other_cell_types=0,vertices=2*n,faces=len(records),
        annular_layer_cells=layer_cells,outer_triangular_cells=nt-layer_cells,annular_layers=NLAYERS,
        wall_segments=NTHETA,first_layer_height_D=FIRST,radial_growth=GROWTH,layer_outer_radius_D=float(radii[-1]),
        wake_target_size_D=wake_size,domain_D=[XMIN,XMAX,YMIN,YMAX],span_thickness_D=THICKNESS,
        min_volume=float(areas.min()*THICKNESS),max_volume=float(areas.max()*THICKNESS),
        min_triangle_angle_deg=float(np.min(corner_angles)),max_edge_aspect_ratio=float((edge_lengths.max(axis=1)/edge_lengths.min(axis=1)).max()),
        max_nonorthogonality_deg=max(nonorth),p95_nonorthogonality_deg=float(np.percentile(nonorth,95)),
        max_normalized_skewness=max(skew),max_boundary_nonorthogonality_deg=max(boundary_nonorth),
        max_circle_sagitta_D=float(R*(1-np.cos(np.pi/NTHETA))),
        area_relative_error=float(abs(areas.sum()-expected_area)/expected_area),
        interface_edges=interface,boundary_face_counts={**patch_counts,'front':nt,'back':nt})
    gates=dict(cell_count=18000<=nt<=23000,positive_volume=bool(np.all(areas>0)),
        minimum_angle=quality['min_triangle_angle_deg']>=22,
        nonorthogonality=quality['max_nonorthogonality_deg']<65,
        skewness=quality['max_normalized_skewness']<.5,
        area=quality['area_relative_error']<1e-11,conformal_interface=interface==NTHETA,
        unique_triangles=len(set(map(tuple,np.sort(tri,axis=1))))==nt,
        wall_complete=len(wall)==NTHETA)
    quality['gates']=gates;quality['passed']=all(gates.values())
    (meshdir/'mesh_quality.json').write_text(json.dumps(quality,indent=2)+'\n')
    np.savez_compressed(meshdir/'layout.npz',xy=xy,tri=tri,centres=centres,areas=areas,edge_cells=edge_cells,wall=wall)
    if not quality['passed']:
        print(json.dumps(quality,indent=2));raise RuntimeError('mesh quality gate failed')
    with (meshdir/'cylinder_re100.mesh').open('w') as out:
        out.write(f'BABELSIM_MESH 3\nvertices\n{2*n}\n')
        for z in (0.,THICKNESS):
            for x,y in xy:out.write(f'{x:.16g} {y:.16g} {z:.16g}\n')
        out.write(f'faces\n{len(records)}\n')
        for ring,o,nb,patch in records:out.write(f'face {len(ring)} {o} {nb} {patch} '+ ' '.join(map(str,ring))+'\n')
        out.write(f'patches\n{len(PATCHES)}\n')
        for name,kind in PATCHES:out.write(f'patch {name} {kind}\n')
        out.write('end\n')
    # Streamfunction psi=(r-R^2/r) sin(theta) f(r), with f=1-exp(-(r-R)/delta).
    # Both velocity components vanish on the wall; divergence is zero analytically.
    x,y=centres.T;r=np.hypot(x,y);ct=x/r;st=y/r;delta=.08
    q=1-np.exp(-(r-R)/delta);dq=np.exp(-(r-R)/delta)/delta
    ur=(1-R*R/r**2)*ct*q
    ut=-st*((1+R*R/r**2)*q+(r-R*R/r)*dq)
    u=ur*ct-ut*st;v=ur*st+ut*ct
    # Compact streamfunction perturbation entirely downstream, no continuous forcing.
    s=(x-1.5)/.8;t=y/.8;mask=(abs(s)<1)&(abs(t)<1)
    u[mask]+=-6*.002*t[mask]*(1-t[mask]**2)**2*(1-s[mask]**2)**3/.8
    v[mask]+=6*.002*s[mask]*(1-s[mask]**2)**2*(1-t[mask]**2)**3/.8
    data=CASE/'fields/initial';data.mkdir(parents=True,exist_ok=True)
    np.savetxt(data/'U.initial.dat',np.column_stack([np.arange(nt),u,v,np.zeros(nt)]),fmt=['%d','%.15g','%.15g','%.15g'],header='globalCellId Ux Uy Uz')
    np.savetxt(data/'p.initial.dat',np.column_stack([np.arange(nt),.5*(1-u*u-v*v)]),fmt=['%d','%.15g'],header='globalCellId p')
    return quality

def plot(xy,tri):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(17,5),constrained_layout=True)
    for ax,bounds,title in zip(axes,[(-10,25,-10,10),(-1,3,-2,2),(-.65,.65,-.65,.65)],['Domain','Near wake','Wall layers']):
        ax.triplot(xy[:,0],xy[:,1],tri,lw=.22,color='#27465a')
        ax.add_patch(plt.Circle((0,0),R,color='#bbbbbb'))
        ax.set(xlim=bounds[:2],ylim=bounds[2:],aspect='equal',xlabel='x/D',ylabel='y/D',title=title)
    fig.suptitle(f'Re=100 cylinder: {len(tri):,} triangular prisms, one spanwise layer')
    fig.savefig(CASE/'mesh/mesh.png',dpi=180);plt.close(fig)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--wake-size',type=float,default=.11);args=parser.parse_args()
    xy,tri,nlayer,radii=generate(args.wake_size)
    quality=audit_and_write(xy,tri,nlayer,radii,args.wake_size)
    plot(xy,tri);print(json.dumps(quality,indent=2))
