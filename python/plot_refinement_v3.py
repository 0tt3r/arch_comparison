"""Figures for parallel compilation, calibration, cultivation, and logical RUS."""

from project_paths import ROOT, project_file
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=ROOT/'results/refinement_v3'


def read(path):
    return json.loads(Path(path).read_text())


def save(fig,name,caption):
    fig.text(.5,.01,caption,ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.95))
    for ext in ('png','pdf'):fig.savefig(OUT/f'{name}.{ext}',dpi=180)
    plt.close(fig)


def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    records=[read(p) for p in sorted(OUT.glob('parallel*/result.json'))]
    chosen=sorted([r for r in records if r['steps']==4 and r['observable_pruned']],key=lambda r:r['L'])
    old=[read(ROOT/f"results/nisq_v2/plaquette_L{r['L']}_r4_fused_cone/result.json") for r in chosen]
    fig,axes=plt.subplots(1,3,figsize=(13,4.5))
    for rows,label,color,style in [(old,'v2 greedy + cone','#999999',':o'),(chosen,'v3 parallel + observable','#0072B2','-o')]:
        for ax,field in zip(axes[:2],('two_qubit_gates','depth_layers')):
            ax.plot([r['L'] for r in rows],[r['schedule'][field] for r in rows],style,color=color,label=label)
    axes[0].set_ylabel('Native entanglers');axes[0].set_yscale('log')
    axes[1].set_ylabel('10 ns primitive layers');axes[0].legend(fontsize=8)
    for ratio,color,label in [('1.0','#0072B2','nominal idle p / 100 ns'),('0.0','#999999','zero idle, same circuit')]:
        axes[2].plot([r['L'] for r in chosen],[r['costs_by_idle_ratio'][ratio]['log10_Q_min'] for r in chosen],
                     '-o',color=color,label=label)
    axes[2].set_ylabel('log₁₀ Q within one week');axes[2].legend(fontsize=8)
    for ax in axes:ax.set_xlabel('Lattice side L');ax.set_xticks([2,4,6,8]);ax.grid(alpha=.2)
    fig.suptitle('Parallel fermionic regrouping reduces both gate count and idle exposure')
    save(fig,'parallel_compilation','Plaquette τ = 1, r = 4; ε = 0.001. These costs remain conditional: r = 4 fails the small-patch Trotter checks.')

    data=read(OUT/'small_dynamics.json')
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    colors={'square2':'#0072B2','square3':'#009E73','cross5':'#D55E00'}
    labels={'square2':'2 × 2','square3':'3 × 3','cross5':'five-site cross'}
    for patch in colors:
        rows=[r for r in data['curves'] if r['patch']==patch]
        axes[0,0].plot([r['tau'] for r in rows],[r['doublon'] for r in rows],color=colors[patch],label=labels[patch])
        rows=[r for r in data['trotter'] if r['patch']==patch and r['tau']==1.]
        axes[0,1].loglog([r['steps'] for r in rows],[r['absolute_trotter_error'] for r in rows],'-o',color=colors[patch],label=labels[patch])
    axes[0,0].set_xlabel('Evolution time τ');axes[0,0].set_ylabel('Exact-patch doublon D');axes[0,0].legend(fontsize=8)
    axes[0,1].axhline(.000125,color='#555555',ls=':',label='strict split Trotter budget')
    axes[0,1].axhline(.00125,color='#999999',ls='--',label='loose split Trotter budget')
    axes[0,1].set_xlabel('Plaquette steps r');axes[0,1].set_ylabel('Local Trotter error at τ = 1');axes[0,1].legend(fontsize=7)
    rows=data['truncation']
    axes[1,0].semilogy([r['tau'] for r in rows],[r['finite_reference_difference'] for r in rows],'-o',color='#D55E00')
    axes[1,0].axhline(.000125,color='#555555',ls=':',label='strict split truncation budget')
    axes[1,0].set_xlabel('Evolution time τ');axes[1,0].set_ylabel('|D(cross5) − D(3 × 3)|');axes[1,0].legend(fontsize=8)
    for L,color in [(2,'#0072B2'),(4,'#D55E00')]:
        rows=sorted([r for r in records if r['L']==L and r['observable_pruned']],key=lambda r:r['steps'])
        axes[1,1].plot([r['steps'] for r in rows],[r['costs_by_idle_ratio']['1.0']['log10_Q_min'] for r in rows],'-o',color=color,label=f'L = {L}')
    axes[1,1].set_xscale('log',base=2);axes[1,1].set_xlabel('Steps r at τ = 1')
    axes[1,1].set_ylabel('log₁₀ Q, compiled local PEC');axes[1,1].legend(fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle('Accuracy and truncation now constrain the resource comparison')
    save(fig,'accuracy_calibration','Exact dynamics only on patches ≤ 9 sites. Finite-reference differences do not bound infinite-volume error; L = 4 accuracy is uncalibrated.')

    cult=read(OUT/'cultivation.json');star=read(OUT/'star_parallel.json')
    fig,axes=plt.subplots(1,2,figsize=(11,4.8))
    for ax,eps in zip(axes,(.001,.01)):
        for arm,color,label in [('ft','#202020','FT + cultivation'),('hybrid','#CC79A7','Hybrid + cultivation')]:
            selected=sorted([r for r in cult if r['absolute_error']==eps and r['arm']==arm
                             and r['source']=='cultivation packed' and r['handoff_error']==0],key=lambda r:r['L'])
            ax.plot([r['L'] for r in selected],[r['log10_Q_min'] for r in selected],'-o',color=color,label=label)
            # Range of declared supply and handoff scenarios, not a confidence interval.
            upper=[max(r['log10_Q_min'] for r in cult if r['absolute_error']==eps and r['arm']==arm
                       and r['L']==s['L'] and r['source'].startswith('cultivation')) for s in selected]
            ax.fill_between([r['L'] for r in selected],[r['log10_Q_min'] for r in selected],upper,color=color,alpha=.12)
        serial=sorted([r for r in cult if r['source']=='STAR reference' and r['absolute_error']==eps],key=lambda r:r['L'])
        ax.plot([r['L'] for r in serial],[r['log10_Q_min'] for r in serial],':o',color='#009E73',label='Serial STAR reference')
        optimized=[min((r for r in star if r['L']==L and r['absolute_error']==eps and r['logical_bias']=='worst_case'),
                       key=lambda r:r['log10_Q_min']) for L in (2,4,6,8)]
        ax.plot([r['L'] for r in optimized],[r['log10_Q_min'] for r in optimized],'-o',color='#009E73',label='Parallel STAR reference')
        ax.set_title(f'Absolute target ε = {eps:g}');ax.set_xlabel('Lattice side L');ax.set_ylabel('log₁₀ Q within one week')
        ax.set_xticks([2,4,6,8]);ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Cultivation helps both encoded arms; retry parallelism has an area cost')
    save(fig,'cultivation_parallel','τ = 1, r = 4 remains conditional. Shading spans declared supply/handoff scenarios. Encoded timing, layout and noise remain reference models.')

    noisy=read(OUT/'noisy_small.json');assert len(noisy)==3
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for name,color,label in [('unselected','#999999','No selection'),('spin_parity','#0072B2','Spin parity'),
                              ('spin_particle_number','#D55E00','Exact spin particle numbers')]:
        axes[0].plot([r['r'] for r in noisy],[r[name]['acceptance'] for r in noisy],'-o',color=color,label=label)
        axes[1].plot([r['r'] for r in noisy],[abs(r[name]['bias_from_noiseless_circuit']) for r in noisy],'-o',color=color,label=label)
    axes[1].axhline(.001,color='#555555',ls=':',label='total reference ε = 0.001')
    axes[1].set_yscale('log')
    for ax in axes:ax.set_xscale('log',base=2);ax.set_xlabel('Steps r');ax.grid(alpha=.2)
    axes[0].set_ylabel('Acceptance');axes[1].set_ylabel('|D(noisy, selected) − D(noiseless circuit)|')
    axes[0].legend(fontsize=8);axes[1].legend(fontsize=8)
    fig.suptitle('Postselection: actual bias on the eight-qubit circuit')
    save(fig,'postselection_bias','L = 2, τ = 1; native gates, idle p / 100 ns, preparation and readout noise. Bias here excludes the separately measured Trotter error.')

    rows=[]
    for r in records:
        rows.append({'L':r['L'],'r':r['steps'],'pruned':r['observable_pruned'],
            **{k:r['schedule'][k] for k in ('two_qubit_gates','one_qubit_gates','depth_layers','idle_qubit_layers','seconds_per_attempt')},
            **{'chi_'+k:v for k,v in r['observable_exposure']['exposure'].items()},
            **{'log10_Q_idle_'+k:v['log10_Q_min'] for k,v in r['costs_by_idle_ratio'].items()}})
    with (OUT/'parallel_summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print('Wrote four figure pairs and compiled summary',flush=True)


if __name__=='__main__':main()
