#!/usr/bin/env python3
"""Check all saved owned-cell fields; solver exit alone is not acceptance."""
import argparse
import json
from pathlib import Path
import re
import numpy as np

def summarize(run):
    report=json.loads((run/'run_summary.json').read_text())
    text=(run/'solver.log').read_text()
    lines=[line for line in text.splitlines() if line.startswith('PISO ')]
    mesh_report=json.loads((run/'mesh'/'mesh_quality.json').read_text())
    cell_count=mesh_report['cell_count']
    steps={}
    folders=sorted((run/'results').iterdir()) if (run/'results').exists() else []
    for folder in folders:
        if not folder.is_dir() or folder.name=='final': continue
        try: float(folder.name)
        except ValueError: continue
        stats={}
        for name in ('U','p','k','omega','mut'):
            files=sorted(folder.glob(f'rank-*/{name}.csv'))
            arrays=[np.loadtxt(f,delimiter=',',skiprows=1,ndmin=2) for f in files]
            if not arrays: raise RuntimeError(f'missing {name} at {folder}')
            data=np.vstack(arrays)
            value=np.linalg.norm(data[:,4:7],axis=1) if name=='U' else data[:,4]
            stats[name]={'min':float(value.min()),'max':float(value.max()),
                         'finite':bool(np.isfinite(data).all()),'cells':len(value),
                         'unique_cells':len(np.unique(data[:,0]))}
        steps[folder.name]=stats
    report['saved_steps']=steps
    report['completed_step_records']=len(lines)
    report['inexact_linear_records']=sum('linear=inexact' in s for s in lines)
    report['max_reported_mass_error']=max([float(re.search(r'mass=([^ ]+)',s)[1]) for s in lines],default=None)
    report['convection']={}
    for line in (run/'numerics/methods.bs').read_text().splitlines():
        if '.convection ' in line:
            key,value=line.split(); report['convection'][key]=value
    report['acceptance']={
        'exit_zero':report['returncode']==0,
        'all_requested_steps':len(lines)==report['steps_requested'] and len(steps)==report['steps_requested'],
        'all_linear_records_ok':report['inexact_linear_records']==0,
        'no_startup_turbulence_floor_clipping':bool(steps) and all(s['k']['min']>1.01e-12 and s['omega']['min']>1.01e-10 for s in steps.values()),
        'finite_fields':bool(steps) and all(s[f]['finite'] for s in steps.values() for f in s),
        'complete_unique_cell_fields':bool(steps) and all(s[f]['cells']==s[f]['unique_cells']==cell_count for s in steps.values() for f in s),
        # Generous sanity limits, not a physical validation or accuracy criterion.
        'bounded_startup_fields':bool(steps) and all(s['U']['max']<10 and s['k']['max']<10 and s['omega']['max']<1e8 and s['mut']['max']<.1 and s['k']['min']>0 and s['omega']['min']>0 for s in steps.values())}
    report['smoke_passed']=all(report['acceptance'].values())
    (run/'field_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path);a=p.parse_args()
    r=summarize(a.run)
    print(json.dumps({k:v for k,v in r.items() if k!='saved_steps'},indent=2))
    raise SystemExit(0 if r['smoke_passed'] else 1)
