#!/usr/bin/env python3
"""Build the measured startup-instability comparison figure and compact evidence."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from postprocess import Geometry,field,CASE

def main():
    geo=Geometry(); reports=[]
    for f in ['stability_study','stability_followup']:reports+=json.loads((CASE/f'validation/{f}.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(13,8),constrained_layout=True)
    labels={'bdf2_dt005':'BDF2 dt=.005','bdf2_dt005_one_corrector':'dt=.005, 1 corrector',
        'euler_dt005':'Euler dt=.005','bdf2_dt005_limited':'dt=.005, limited diffusion',
        'bdf2_dt005_upwind':'dt=.005, upwind','bdf2_dt002':'BDF2 dt=.002',
        'bdf2_dt005_greenGauss':'dt=.005, Green-Gauss','bdf2_dt005_orthogonal':'dt=.005, orthogonal only',
        'bdf2_dt005_alpha07':'dt=.005, relaxation=.7'}
    for report in reports:
        samples=[x for x in report['samples'] if x['max_abs_p']<1e5]
        axes[0,0].semilogy([x['t'] for x in samples],[x['max_abs_p'] for x in samples],label=labels[report['name']],lw=1)
    stable=CASE/'results/pilot_dt0001/results/run';samples=[]
    for p in stable.iterdir():
        try:t=float(p.name)
        except ValueError:continue
        pressure=field(p,'p',geo.n);samples.append((t,float(abs(pressure).max())))
    samples.sort();a=np.asarray(samples)
    axes[0,0].semilogy(a[:,0],a[:,1],'k.-',label='BDF2 dt=.001',lw=1.8)
    axes[0,0].set(xlim=(0,.4),ylim=(.5,1e5),xlabel='t U/D',ylabel='max |p| / (rho U^2)',title='Same mesh: startup stability');axes[0,0].legend(fontsize=7,ncol=2)
    j=int(np.argmin(np.sum((geo.c-np.array([0.,.505]))**2,axis=1)))
    oscillation=[]
    for t in [.025,.05,.075,.1,.125,.15]:
        p=field(CASE/f'results/stability/bdf2_dt005/results/run/{t:g}','p',geo.n)
        oscillation.append([t,float(p[j,0])])
    a=np.asarray(oscillation);axes[0,1].plot(a[:,0],a[:,1],'o-',color='crimson')
    axes[0,1].axhline(0,color='gray',lw=.5);axes[0,1].set(xlabel='t U/D',ylabel='p / (rho U^2)',title=f'Near-wall pressure at ({geo.c[j,0]:.4f}, {geo.c[j,1]:.4f})')
    for ax,folder,title,limit in [(axes[1,0],CASE/'results/stability/bdf2_dt005/results/run/0.1','dt=.005, t=.1',250),
        (axes[1,1],stable/'0.1','dt=.001, t=.1',1.6)]:
        p=field(folder,'p',geo.n)[:,0]
        im=ax.tripcolor(geo.xy[:,0],geo.xy[:,1],geo.tri,facecolors=p,cmap='RdBu_r',vmin=-limit,vmax=limit)
        ax.add_patch(plt.Circle((0,0),.5,color='.65',zorder=3));ax.set(xlim=(-.9,.9),ylim=(-.9,.9),aspect='equal',xlabel='x/D',ylabel='y/D',title=title+' (different color scales)');fig.colorbar(im,ax=ax)
    for ax in axes[0]:ax.grid(alpha=.2)
    fig.savefig(CASE/'validation/stability_comparison.png',dpi=170);plt.close(fig)
    ncorrect=[]
    for t in ['0.025','0.05','0.1','0.15']:
        base=CASE/'results/stability';a=field(base/f'bdf2_dt005/results/run/{t}','U',geo.n);b=field(base/f'bdf2_dt005_one_corrector/results/run/{t}','U',geo.n)
        ncorrect.append({'t':float(t),'relative_L2_U_difference':float(np.linalg.norm(a-b)/np.linalg.norm(a))})
    (CASE/'validation/mechanism_evidence.json').write_text(json.dumps({'near_wall_cell':j,'near_wall_xy':geo.c[j].tolist(),
        'near_wall_pressure_series':oscillation,'one_vs_three_correctors':ncorrect},indent=2)+'\n')
if __name__=='__main__':main()
