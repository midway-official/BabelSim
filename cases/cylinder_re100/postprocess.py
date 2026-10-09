#!/usr/bin/env python3
"""Read owned MPI CSV fields; plot signed vorticity, probes, forces and spectra.
Vorticity is a weighted least-squares cell gradient, with no-slip wall samples.
Force integration is a first-layer estimate and does not replace grid convergence.
"""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
CASE=Path(__file__).resolve().parent


def field(path,name,n):
    files=sorted(path.glob(f'rank-*/{name}.csv'))
    if len(files)!=2:raise ValueError(f'Expected two ranks at {path}')
    a=np.concatenate([np.loadtxt(p,delimiter=',',skiprows=1,ndmin=2) for p in files])
    a=a[np.argsort(a[:,0])]
    if len(a)!=n or not np.array_equal(a[:,0],np.arange(n)) or not np.isfinite(a).all():
        raise ValueError(f'Incomplete/nonfinite {name} at {path}')
    return a[:,4:]

class Geometry:
    def __init__(self):
        data=np.load(CASE/'mesh/layout.npz')
        self.xy=data['xy'];self.tri=data['tri'];self.c=data['centres'];self.area=data['areas'];self.n=len(self.c)
        edges=data['edge_cells'];self.o=edges[:,0];self.nb=edges[:,1]
        self.sf=np.column_stack([(self.xy[edges[:,3]]-self.xy[edges[:,2]])[:,1],-(self.xy[edges[:,3]]-self.xy[edges[:,2]])[:,0]])
        wall=data['wall'];self.wo=wall[:,2];self.wc=(self.xy[wall[:,0]]+self.xy[wall[:,1]])/2
        d=self.xy[wall[:,1]]-self.xy[wall[:,0]];self.wlen=np.linalg.norm(d,axis=1)
        self.wn=np.column_stack([d[:,1],-d[:,0]])/self.wlen[:,None] # fluid outward (into solid)
        self.wdist=np.sum((self.wc-self.c[self.wo])*self.wn,axis=1)
        self.i=np.r_[self.o,self.nb,self.wo]
        self.j=np.r_[self.nb,self.o]
        delta=np.vstack([self.c[self.nb]-self.c[self.o],self.c[self.o]-self.c[self.nb],self.wc-self.c[self.wo]])
        self.w=1/np.sum(delta*delta,axis=1)
        self.delta=delta
        self.matrix=np.zeros((self.n,2,2))
        np.add.at(self.matrix,self.i,self.w[:,None,None]*delta[:,:,None]*delta[:,None,:])
        assert np.all(np.linalg.det(self.matrix)>0)

    def gradient(self,u,wall_values=None):
        if wall_values is None:wall_values=np.zeros((len(self.wo),u.shape[1]))
        diff=np.vstack([u[self.j]-u[self.i[:len(self.j)]],wall_values-u[self.wo]])
        rhs=np.zeros((self.n,2,u.shape[1]))
        np.add.at(rhs,self.i,self.w[:,None,None]*self.delta[:,:,None]*diff[:,None,:])
        return np.linalg.solve(self.matrix,rhs)

    def diagnostics(self,u,p,dt):
        grad=self.gradient(u[:,:2]);omega=grad[:,0,1]-grad[:,1,0]
        flux=abs(np.sum(.5*(u[self.o,:2]+u[self.nb,:2])*self.sf,axis=1))
        accum=np.zeros(self.n);np.add.at(accum,self.o,flux);np.add.at(accum,self.nb,flux)
        # Interior flux proxy; boundary contributions are zero at the wall and far from max Co.
        co=.5*dt*accum/self.area
        vel=u[self.wo,:2];tangent=vel-np.sum(vel*self.wn,axis=1)[:,None]*self.wn
        traction=p[self.wo]*self.wn+.01*tangent/self.wdist[:,None]
        force=np.sum(traction*self.wlen[:,None],axis=0)/.5 # span cancels
        return omega,co,force

    def check_affine_gradient(self):
        a=np.array([[1.2,-.4],[.7,2.1]])
        u=self.c@a+np.array([.3,-.1]);wu=self.wc@a+np.array([.3,-.1])
        err=float(np.max(abs(self.gradient(u,wu)-a)))
        assert err<1e-10
        return err

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--results',type=Path,default=CASE/'results/run')
    parser.add_argument('--out',type=Path,default=CASE/'validation');parser.add_argument('--dt',type=float,default=.001)
    parser.add_argument('--gif',action='store_true');args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    geo=Geometry();graderr=geo.check_affine_gradient();paths=[]
    for p in args.results.iterdir():
        try:t=float(p.name)
        except ValueError:continue
        if not p.is_dir():continue
        try:u=field(p,'U',geo.n);pressure=field(p,'p',geo.n)
        except (ValueError,OSError):continue
        paths.append((t,p))
    paths.sort()
    if not paths:raise RuntimeError('No complete finite two-rank snapshots')
    targets=np.array([[2.,.5],[4.,.5],[6.,0.],[10.,.5]])
    probes=np.argmin(np.sum((geo.c[:,None,:]-targets[None,:,:])**2,axis=2),axis=0)
    rows=[];last_omega=None;last_u=None;last_p=None;last_co=None
    for t,path in paths:
        u=field(path,'U',geo.n);p=field(path,'p',geo.n)
        omega,co,force=geo.diagnostics(u,p,args.dt)
        rows.append([t,*u[probes,1],*force,float(np.max(np.linalg.norm(u,axis=1))),float(np.max(abs(p))),float(np.max(co))])
        last_omega,last_u,last_p,last_co=omega,u,p,co
    rows=np.asarray(rows)
    np.savetxt(args.out/'probes_forces.csv',rows,delimiter=',',header='time,Uy_2_0.5,Uy_4_0.5,Uy_6_0,Uy_10_0.5,Cd_estimate,Cl_estimate,max_speed,max_abs_p,max_Co_interior',comments='')
    fig,axes=plt.subplots(2,1,figsize=(13,7),constrained_layout=True)
    for ax,values,limit,title,cmap in [(axes[0],last_omega,2.,r'Signed vorticity $\omega_z D/U_\infty$','RdBu_r'),
        (axes[1],np.linalg.norm(last_u,axis=1),None,r'Speed $|U|/U_\infty$','viridis')]:
        im=ax.tripcolor(geo.xy[:,0],geo.xy[:,1],geo.tri,facecolors=values,shading='flat',cmap=cmap,
            **({'vmin':-limit,'vmax':limit} if limit else {'vmin':0,'vmax':1.6}))
        ax.add_patch(plt.Circle((0,0),.5,color='.6',zorder=3));ax.set(xlim=(-2,16),ylim=(-3,3),aspect='equal',xlabel='x/D',ylabel='y/D',title=title)
        fig.colorbar(im,ax=ax,extend='both')
    fig.suptitle(f'Re=100 | PISO/BDF2, laminar | t U/D={paths[-1][0]:g}')
    fig.savefig(args.out/'flow_latest.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,9),constrained_layout=True)
    for k,target in enumerate(targets):axes[0].plot(rows[:,0],rows[:,1+k],label=f'{tuple(target)}')
    axes[0].set(ylabel=r'$U_y/U_\infty$',title='Wake probes (actual nearest-cell locations in diagnostics.json)');axes[0].legend(ncol=4)
    axes[1].plot(rows[:,0],rows[:,5],label='Cd estimate');axes[1].plot(rows[:,0],rows[:,6],label='Cl estimate')
    axes[1].legend();axes[1].set(ylabel='Force coefficient',title='Wall pressure + first-layer tangential viscous traction')
    spectrum={'status':'insufficient developed-flow duration'}
    sample=rows[rows[:,0]>=max(50.,.5*rows[-1,0])]
    if len(sample)>=32 and sample[-1,0]-sample[0,0]>=20:
        ts=sample[:,0];signal=sample[:,2];step=float(np.median(np.diff(ts)))
        if not np.allclose(np.diff(ts),step,atol=1e-8):raise ValueError('Nonuniform PSD samples')
        y=signal-np.polyval(np.polyfit(ts,signal,1),ts);window=np.hanning(len(y));power=abs(np.fft.rfft(y*window))**2
        freq=np.fft.rfftfreq(len(y),step);valid=(freq>.05)&(freq<.4)
        peak=int(np.flatnonzero(valid)[np.argmax(power[valid])]);axes[2].plot(freq,power)
        axes[2].set(xlim=(0,.5),xlabel=r'$fD/U_\infty=St$',ylabel='Windowed spectral power')
        spectrum={'status':'finite-window estimate; inspect periodicity and refine dt/mesh for validation',
            'St_peak':float(freq[peak]),'frequency_bin_width':float(freq[1]-freq[0]),'window':[float(ts[0]),float(ts[-1])],
            'probe_Uy_rms':float(np.std(signal)),'cycles_in_window':float(freq[peak]*(ts[-1]-ts[0]))}
    else:axes[2].text(.5,.5,'More developed-flow samples needed for a shedding spectrum',ha='center',transform=axes[2].transAxes)
    for ax in axes:ax.grid(alpha=.2);ax.set_xlabel(r'$tU_\infty/D$' if ax!=axes[2] else ax.get_xlabel())
    fig.savefig(args.out/'probes_spectrum.png',dpi=160);plt.close(fig)
    report={'cells':geo.n,'ranks_per_snapshot':2,'complete_finite_snapshots':len(paths),'first_time':paths[0][0],'last_time':paths[-1][0],
        'affine_gradient_max_error':graderr,'probe_cell_ids':probes.tolist(),'actual_probe_xy':geo.c[probes].tolist(),
        'latest_U_min':last_u.min(axis=0).tolist(),'latest_U_max':last_u.max(axis=0).tolist(),
        'latest_pressure_range':[float(last_p.min()),float(last_p.max())],'max_speed_over_saved_times':float(rows[:,7].max()),
        'max_abs_pressure_over_saved_times':float(rows[:,8].max()),'max_Co_interior_over_saved_times':float(rows[:,9].max()),
        'latest_omega_range':[float(last_omega.min()),float(last_omega.max())],'spectrum':spectrum,
        'physical_validation':'not established by this script; requires developed periodicity and sensitivity checks'}
    (args.out/'diagnostics.json').write_text(json.dumps(report,indent=2)+'\n')
    if args.gif:
        chosen=paths[max(0,len(paths)-100):][::2]
        fig,ax=plt.subplots(figsize=(11,4),constrained_layout=True)
        artist=ax.tripcolor(geo.xy[:,0],geo.xy[:,1],geo.tri,facecolors=last_omega,cmap='RdBu_r',vmin=-2,vmax=2)
        ax.add_patch(plt.Circle((0,0),.5,color='.6',zorder=3));ax.set(xlim=(-2,16),ylim=(-3,3),aspect='equal',xlabel='x/D',ylabel='y/D')
        fig.colorbar(artist,ax=ax,label=r'$\omega_z D/U_\infty$',extend='both')
        def update(k):
            t,path=chosen[k];u=field(path,'U',geo.n);g=geo.gradient(u[:,:2]);artist.set_array(g[:,0,1]-g[:,1,0]);ax.set_title(f'Re=100 cylinder, signed vorticity, tU/D={t:g}');return artist,
        FuncAnimation(fig,update,frames=len(chosen),blit=False).save(args.out/'vorticity.gif',writer=PillowWriter(fps=8),dpi=110);plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
