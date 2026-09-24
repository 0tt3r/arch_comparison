"""Invert saved week-cutoff costs into maximum tested lattice size versus Q.

The staircase is a discrete reachability envelope. It does not interpolate
untested lattice sizes or certify the common four-step Trotter circuit.
"""
import csv
import json
from math import ceil, exp, log10
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from model import accuracy, finish, workload
from project_paths import ROOT, project_file

COLORS={'nisq':'#168C9E','star':'#D9872C','ft':'#6750A4','hybrid':'#56815B'}
LABELS={'nisq':'Physical · gatewise PEC bound','star':'STAR · parallel retries',
        'ft':'FTQC · cultivated T states','hybrid':'Hybrid · fixed T + analog'}


def maximum_tested_side(points,qubits):
    """Largest tested L with required area <= Q; NaN means none is feasible.

    Costs need not be monotone in L. This operation takes the upper envelope,
    rather than pretending that swapping plot axes defines a single-valued curve.
    """
    if qubits<=0:raise ValueError('Q must be positive')
    available=log10(qubits)
    feasible=[p['L'] for p in points if p['log10_Q_min']<=available+1e-12]
    return max(feasible) if feasible else float('nan')


def read(relative):
    return json.loads((ROOT/relative).read_text())


def costs():
    config=json.loads(project_file('assumptions.json').read_text())
    cult=read('results/refinement_v3/cultivation.json')
    star=read('results/refinement_v3/star_parallel.json')
    rows=[]
    for eps in (.001,.01):
        scenario={**config,'relative_error':eps/config['doublon_reference']}
        for L in (2,4,6,8):
            native=read(f'results/refinement_v3/parallel_L{L}_observable/result.json')
            w=workload(L,1.,4,'plaquette',config['u_over_j'])
            row=finish(w,'nisq',native['observable_exposure']['chi'],
                       native['schedule']['physical_copy_qubits'],native['schedule']['seconds_per_attempt'],scenario)
            rows.append({**row,'absolute_error':eps,'reference':'compiled native, nominal idle'})
            best=min((r for r in star if r['L']==L and r['absolute_error']==eps and
                      r['logical_bias']=='worst_case'),key=lambda r:r['log10_Q_min'])
            rows.append({**best,'arm':'star','reference':'best of 1, 2, 4, 8 reference lanes'})
            for arm in ('ft','hybrid'):
                row=next(r for r in cult if r['L']==L and r['arm']==arm and
                         r['absolute_error']==eps and r['source']=='cultivation packed' and r['handoff_error']==0)
                rows.append({**row,'reference':'packed cultivation; zero additional handoff error'})
    return rows


def nisq_sampling_audit():
    """Decompose the L=4 headline; noise ablations are diagnostic scenarios."""
    config=json.loads(project_file('assumptions.json').read_text())
    native=read('results/refinement_v3/parallel_L4_observable/result.json')
    schedule=native['schedule'];exposure=native['observable_exposure']
    eps,n0=accuracy(config);q=schedule['physical_copy_qubits'];t=schedule['seconds_per_attempt']
    assert abs(sum(exposure['exposure'].values())-exposure['chi'])<1e-12
    scenarios=[]
    for name,chi,scope in [
        ('reference',exposure['chi'],'Current sufficient PEC shot budget'),
        ('zero_idle',exposure['chi']-exposure['exposure']['idle'],'Noise sensitivity; same gates and duration'),
        ('two_qubit_noise_only',exposure['exposure']['two_qubit'],'Diagnostic: all other noise channels removed'),
        ('no_mitigation_multiplier',0.,'Shot throughput only; hardware bias not corrected'),
    ]:
        shots=n0*exp(chi)
        copies=max(1,ceil(shots*t/config['week_seconds']))
        row=finish(workload(4,1.,4,'plaquette',config['u_over_j']),'nisq',chi,q,t,config)
        assert copies*q==round(10**row['log10_Q_min'])
        scenarios.append(dict(scenario=name,scope=scope,chi=chi,gamma2=exp(chi),
                              shot_budget=ceil(shots),single_copy_seconds=shots*t,
                              copies=copies,physical_qubit_budget=copies*q))
    return dict(L=4,tau=1.,r=4,absolute_error=eps,statistical_error=eps*config['stat_fraction'],
                confidence=1-config['delta'],ideal_shot_bound=n0,circuit_qubits=schedule['width'],
                qubits_per_copy=q,seconds_per_attempt=t,exposure=exposure['exposure'],scenarios=scenarios,
                interpretation='Sufficient budgets for a specified PEC implementation, not necessary NISQ resources')


def main():
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    data=ROOT/'results/presentation_v1';data.mkdir(exist_ok=True)
    rows=costs()
    (data/'nisq_sampling_audit.json').write_text(json.dumps(nisq_sampling_audit(),indent=2)+'\n')
    (data/'capacity_thresholds.json').write_text(json.dumps(rows,indent=2)+'\n')
    fields=['absolute_error','arm','L','N','log10_Q_min','reference']
    with (data/'capacity_thresholds.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(12.6,5.4),sharey=True)
    xmin,xmax=8.,1e9
    exported=[]
    for ax,eps,title in zip(axes,(.001,.01),('1% relative at D* = 0.1','0.01 absolute sensitivity')):
        ax.axhspan(7.75,8.25,color='#f1f1f1',zorder=0)
        for arm in ('nisq','star','ft','hybrid'):
            points=[r for r in rows if r['absolute_error']==eps and r['arm']==arm]
            thresholds=[10**r['log10_Q_min'] for r in points if log10(xmin)<=r['log10_Q_min']<=log10(xmax)]
            grid=sorted(set([xmin,xmax]+thresholds))
            values=[maximum_tested_side(points,q) for q in grid]
            ax.step(grid,values,where='post',label=LABELS[arm],color=COLORS[arm],
                    lw=2.5 if arm!='hybrid' else 2.,ls='-' if arm!='hybrid' else '--',zorder=3)
            # Mark only actual upper-envelope transitions, without hiding
            # equal-threshold jumps or connecting the unachievable region.
            prev=float('nan')
            for q,value in zip(grid,values):
                if np.isfinite(value) and (not np.isfinite(prev) or value!=prev):
                    ax.plot(q,value,'o',ms=4,color=COLORS[arm],zorder=4)
                exported.append({'absolute_error':eps,'arm':arm,'Q':q,'L_max_tested':value if np.isfinite(value) else None})
                prev=value
        ax.set_xscale('log');ax.set_xlim(xmin,xmax);ax.set_ylim(1.45,8.5)
        ax.set_xticks([10**k for k in (1,3,5,7,9)])
        ax.set_yticks([2,4,6,8],['2 × 2','4 × 4','6 × 6','8 × 8'])
        ax.grid(axis='x',alpha=.14);ax.grid(axis='y',alpha=.23)
        ax.set_xlabel('Available physical qubits, Q',labelpad=10)
        ax.set_title(title,loc='left',fontsize=13,pad=12)
        ax.text(.98,.965,'8 × 8: largest size tested',transform=ax.transAxes,ha='right',va='top',fontsize=8.5,color='#666666')
    axes[0].set_ylabel('Largest tested lattice within one week\n(reference implementations)',labelpad=12)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.072),ncol=2,frameon=False,fontsize=9.5)
    fig.suptitle('Hubbard capacity of the reference implementations',fontsize=17,x=.085,ha='left',y=.985)
    fig.text(.085,.91,'Conditional circuit costs: τ = 1, r = 4, U/J = 4, p = 10⁻³',fontsize=10,color='#555555')
    fig.text(.5,.024,'The physical arm uses a sufficient gatewise PEC budget, not an optimized hardware limit. Steps show tested sizes only.\n'
             'Four-step accuracy is unvalidated at L ≥ 4 and fails the smaller-patch checks.',
             ha='center',fontsize=8)
    fig.subplots_adjust(left=.085,right=.985,bottom=.27,top=.79,wspace=.10)
    for suffix in ('pdf','svg','png'):
        path=out/f'capacity_vs_qubits.{suffix}'
        fig.savefig(path,dpi=210)
        if suffix=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)
    (data/'capacity_frontiers.json').write_text(json.dumps(exported,indent=2,allow_nan=False)+'\n')
    print('Wrote capacity figure and 32 underlying cost thresholds; no interpolation beyond tested sizes.')


if __name__=='__main__':main()
