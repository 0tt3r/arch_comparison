"""Compare actual compiled circuits with v1 proxies and named sensitivities."""
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from compiled_model import compiled_physical, compiled_star, encoded
from model import estimate, workload

ROOT=Path(__file__).resolve().parent


def main():
    c=json.loads((ROOT/'assumptions.json').read_text())
    out=ROOT/'results/cudaq_v1'
    records=[json.loads(p.read_text()) for p in sorted(out.glob('*_L*_tau*_r*/result.json'))]
    rows=[]
    for record in records:
        w=workload(record['L'],record['tau'],record['steps'],record['algorithm'])
        cases={'physical_v1_proxy':estimate(w,c,'nisq'),
               'physical_compiled_layer':compiled_physical(record,c,'per_layer'),
               'physical_compiled_time':compiled_physical(record,c,'per_100ns'),
               'physical_fswap_time':compiled_physical(record,c,'per_100ns','fswap'),
               'FT_v1_fit':estimate(w,c,'ft'),
               'FT_compiled':encoded(record,c),
               'FT_one_clock_sensitivity':encoded(record,c,t_clocks=1),
               'STAR_compiled':compiled_star(record,c),
               'hybrid':encoded(record,c,'hybrid'),
               'STAR_commuting_faults':compiled_star(record,c,'commuting_pauli'),
               'hybrid_commuting_faults':encoded(record,c,'hybrid','commuting_pauli')}
        rows.extend({'scenario':name,**row} for name,row in cases.items())
    fields=sorted(set().union(*(row.keys() for row in rows)))
    with (out/'resource_comparison.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    (out/'analysis_runtime.json').write_text(json.dumps({
        'python':sys.version,'executable':sys.executable,
        'packages':{name:importlib.metadata.version(name) for name in ('numpy','scipy','matplotlib')},
        'configuration':c,
        'sources':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                   for name in ('model.py','compiled_model.py','plot_compiled.py')}
    },indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    physical=[('physical_v1_proxy','v1 gate/depth proxy','#999999',':'),
              ('physical_compiled_layer','Compiled JW; idle p / layer','#D55E00','--'),
              ('physical_compiled_time','Compiled JW; idle p / 100 ns','#0072B2','-')]
    physical.append(('physical_fswap_time','Greedy fSWAP; idle p / 100 ns','#CC79A7','-'))
    protected=[('FT_v1_fit','FT: analytical T fit','#999999',':'),
               ('FT_compiled','FT: CUDA-Q T count','#202020','-'),
               ('STAR_compiled','STAR: all angles injected','#009E73','--'),
               ('hybrid','Hybrid: fixed angles distilled','#CC79A7','-')]
    for i,alg in enumerate(('plaquette','termwise')):
        for j,series in enumerate((physical,protected)):
            ax=axes[i,j]
            for name,label,color,style in series:
                chosen=sorted([r for r in rows if r['algorithm']==alg and r['scenario']==name],key=lambda r:r['L'])
                ax.plot([r['L'] for r in chosen],[r['log10_Q_min'] for r in chosen],
                        label=label,color=color,ls=style,marker='o',lw=1.7)
                if j==0:
                    clipped=[r for r in chosen if r['log10_Q_min']>100]
                    ax.scatter([r['L'] for r in clipped],[99]*len(clipped),
                               marker='^',color=color,s=24,clip_on=True)
            ax.set_title(('Plaquette' if i==0 else 'Termwise Strang') + (' — physical PEC' if j==0 else ' — encoded references'))
            ax.set_xlabel('Lattice side L')
            ax.set_ylabel('log₁₀ Q required within one week')
            ax.set_xticks([2,4,6,8])
            ax.grid(alpha=.2)
            if j==0:
                ax.set_ylim(0,100)
                ax.text(.02,.96,'▲ exceeds plotted range; full values in CSV',
                        transform=ax.transAxes,fontsize=8,va='top')
            else:
                ax.set_ylim(3,10)
            if i==0: ax.legend(fontsize=8)
    fig.suptitle('CUDA-Q circuit counts: routing correction and a first hybrid test',fontsize=13)
    fig.text(.5,.012,'τ = 1, r = 4; 1% relative at D* = 0.1 (ε = 0.001); p = 10⁻³; full-lattice PEC\n'
             'Physical path is executable NN routing, not optimized. Encoded timing/factories remain reference models.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.065,1,.95))
    for ext in ('png','pdf'):
        fig.savefig(out/f'compiled_comparison.{ext}',dpi=180)
    plt.close(fig)
    # Counts and depth remain more informative than the huge whole-circuit
    # PEC extrapolations once the latter leave any plausible resource range.
    fig,axes=plt.subplots(1,2,figsize=(10,4.5))
    for alg,color in [('plaquette','#0072B2'),('termwise','#D55E00')]:
        selected=sorted([r for r in records if r['algorithm']==alg],key=lambda r:r['L'])
        for route,key,style in [('JW','native_schedule','-'),('fSWAP','fswap_schedule','--')]:
            for ax,field,scale in [(axes[0],'two_qubit_gates',1),
                                   (axes[1],'seconds_per_attempt',1e6)]:
                ax.plot([r['L'] for r in selected],[r[key][field]*scale for r in selected],
                        color=color,ls=style,marker='o',label=f'{alg}, {route}')
    for ax in axes:
        ax.set_yscale('log');ax.set_xticks([2,4,6,8]);ax.set_xlabel('Lattice side L');ax.grid(alpha=.2)
    axes[0].set_ylabel('Nearest-neighbor CNOT count');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('ASAP attempt duration (µs)')
    fig.suptitle('Same fermionic blocks, two explicit NN implementations; τ = 1, r = 4')
    fig.tight_layout()
    for ext in ('png','pdf'):
        fig.savefig(out/f'compilation_diagnostics.{ext}',dpi=180)
    plt.close(fig)
    print('Wrote',out/'compiled_comparison.png')
    for row in rows:
        if row['L']==4 and row['algorithm']=='plaquette':
            print(row['scenario'],round(row['log10_Q_min'],3),'d=',row.get('distance'))


if __name__ == '__main__':
    main()
