"""Plot compiled NISQ refinements, exact parity acceptance, and injection sensitivity."""
import csv
import importlib.metadata
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from compiled_model import compiled_physical

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/nisq_v2'


def save(fig,name):
    for ext in ('png','pdf'):fig.savefig(OUT/f'{name}.{ext}',dpi=180)
    plt.close(fig)


def main():
    records=[json.loads(p.read_text()) for p in sorted(OUT.glob('*_L*_r*_*/result.json'))]
    boundary={r['case']:r for r in json.loads((OUT/'boundary_timing.json').read_text())}
    base=json.loads((ROOT/'assumptions.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11,4.6))
    colors={'plaquette':'#0072B2','termwise':'#D55E00'}
    for alg in colors:
        chosen=sorted([r for r in records if r['algorithm']==alg and r['steps']==4
                       and r['variant']=='fused_cone'],key=lambda r:r['L'])
        for ax,field,scale in [(axes[0],'two_qubit_gates',1),(axes[1],'seconds_per_attempt',1e6)]:
            ax.plot([r['L'] for r in chosen],[r['schedule'][field]*scale for r in chosen],
                    '-o',color=colors[alg],label=f'{alg}: fused + cone')
            old=[json.loads((ROOT/f"results/cudaq_v1/{alg}_L{r['L']}_tau1_r4/result.json").read_text()) for r in chosen]
            ax.plot([r['L'] for r in old],[r['native_schedule'][field]*scale for r in old],
                    ':o',color=colors[alg],label=f'{alg}: v1 JW')
    for ax in axes:ax.set_yscale('log');ax.set_xlabel('Lattice side L');ax.set_xticks([2,4,6,8]);ax.grid(alpha=.2)
    axes[0].set_ylabel('Two-qubit gates');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('Attempt duration (µs)')
    fig.suptitle('Native compilation: same ordered logical blocks, τ = 1, r = 4')
    fig.text(.5,.01,'v1 counts native CNOTs; v2 counts native RZZ gates and explicit one-qubit pulses. Both entanglers cost 10 ns.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.055,1,.95));save(fig,'native_compilation')

    plaq=sorted([r for r in records if r['algorithm']=='plaquette' and r['steps']==4
                 and r['variant']=='fused_cone'],key=lambda r:r['L'])
    fig,axes=plt.subplots(1,2,figsize=(11,4.6))
    for ratio,color in [('0.0','#999999'),('0.1','#009E73'),('1.0','#0072B2')]:
        axes[0].plot([r['L'] for r in plaq],[r['costs_by_idle_ratio'][ratio]['log10_Q_min'] for r in plaq],
                     '-o',color=color,label=f'idle error = {ratio}p per 100 ns')
    axes[0].set_xlabel('Lattice side L');axes[0].set_ylabel('log₁₀ Q required within one week')
    axes[0].set_xticks([2,4,6,8]);axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    axes[0].set_title('Plaquette: specified local PEC estimator')
    current=next(r for r in plaq if r['L']==4)
    old=json.loads((ROOT/'results/cudaq_v1/plaquette_L4_tau1_r4/result.json').read_text())
    rows=[compiled_physical(old,base,routing=route) for route in ('jw','fswap')]
    labels=['v1 JW','v1 fSWAP','v2 fused + cone'];bottom=[0.,0.,0.]
    for key,color in [('two_qubit','#0072B2'),('one_qubit','#E69F00'),('idle','#D55E00'),('spam','#999999')]:
        heights=[r['chi_'+key] for r in rows]+[current['observable_exposure']['exposure'][key]]
        axes[1].bar(labels,heights,bottom=bottom,color=color,label=key.replace('_',' '))
        bottom=[x+y for x,y in zip(bottom,heights)]
    axes[1].set_ylabel('χ = log variance multiplier');axes[1].set_title('L = 4; idle p per 100 ns')
    axes[1].legend(fontsize=8)
    fig.suptitle('NISQ costs depend on compilation and idle exposure')
    fig.text(.5,.01,'1% relative at D* = 0.1 (ε = 0.001); 95% statistical confidence. Trotter bias uncalibrated; no postselection credit.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.055,1,.95));save(fig,'nisq_resources')

    cone=json.loads((OUT/'cone_sweep.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for steps,color in [(1,'#0072B2'),(2,'#009E73'),(4,'#D55E00')]:
        chosen=[r for r in cone if r['algorithm']=='plaquette' and r['r']==steps]
        axes[0,0].plot([r['L'] for r in chosen],[r['cone_modes']/2 for r in chosen],'-o',color=color,label=f'plaquette r = {steps}')
    chosen=[r for r in cone if r['algorithm']=='termwise' and r['r']==1]
    axes[0,0].plot([r['L'] for r in chosen],[r['cone_modes']/2 for r in chosen],':',color='#555555',label='termwise r = 1 (whole array)')
    axes[0,0].set_ylabel('Sites in sufficient fermionic cone');axes[0,0].set_title('Exact circuit dependency, no spatial tolerance')
    axes[0,0].legend(fontsize=8)
    one=sorted([r for r in records if r['steps']==1 and r['variant']=='fused_cone'],key=lambda r:r['L'])
    axes[1,0].plot([r['L'] for r in one],[r['schedule']['two_qubit_gates'] for r in one],'-o',color='#0072B2')
    axes[1,0].set_ylabel('Compiled RZZ gates');axes[1,0].set_title('Plaquette r = 1: compilation saturates')
    axes[0,1].plot([r['L'] for r in plaq],[r['postselection']['acceptance'] for r in plaq],'-o',color='#0072B2',label='common preparation/readout')
    axes[0,1].axhline(.25,color='#999999',ls=':',label='fully randomized parity limit')
    axes[0,1].set_ylim(.2,1);axes[0,1].set_ylabel('Acceptance: both spin parities');axes[0,1].set_title('Plaquette r = 4: postselection')
    axes[0,1].legend(fontsize=8)
    axes[1,1].plot([r['L'] for r in plaq],[r['postselection']['conditional_bias_upper_bound'] for r in plaq],'-o',color='#D55E00')
    axes[1,1].axhline(.001,color='#555555',ls=':',label='total target ε = 0.001')
    axes[1,1].set_yscale('log');axes[1,1].set_ylabel('Sufficient worst-case bias bound')
    axes[1,1].set_title('Acceptance does not certify precision');axes[1,1].legend(fontsize=8)
    for ax in axes.flat:ax.set_xlabel('Lattice side L');ax.grid(alpha=.2)
    fig.suptitle('Locality and postselection constrain different parts of the calculation')
    fig.text(.5,.01,'Bounds are sufficient, not measured biases. Cone restriction applies to all architectures; r is not Trotter-calibrated.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.045,1,.95));save(fig,'locality_postselection')

    factory=json.loads((OUT/'factory_sensitivity.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(11,4.7))
    series=[('STAR',None,'STAR','#009E73','-'),('hybrid',1.,'Hybrid: p_in = p','#CC79A7',':'),
            ('hybrid',46/15,'Hybrid: p_in = 46p/15','#CC79A7','-'),
            ('ft',1.,'FT: p_in = p','#202020',':'),('ft',46/15,'FT: p_in = 46p/15','#202020','-')]
    for ax,eps in zip(axes,(.001,.01)):
        for arm,pin,label,color,style in series:
            chosen=sorted([r for r in factory if r['absolute_error']==eps and r['scenario']==arm
                           and (pin is None or r.get('injection_factor')==pin)],key=lambda r:r['L'])
            ax.plot([r['L'] for r in chosen],[r['log10_Q_min'] for r in chosen],style+'o',color=color,label=label)
        ax.set_title(f'Absolute target ε = {eps:g}');ax.set_xlabel('Lattice side L')
        ax.set_ylabel('log₁₀ Q required within one week');ax.set_xticks([2,4,6,8]);ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('The hybrid gain depends on injected-state quality and precision')
    fig.text(.5,.01,'Plaquette τ = 1, r = 4. Same compiled FT words in both panels; encoded timings are reference models.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.055,1,.95));save(fig,'factory_reassessment')

    rows=[]
    for r in records:
        case=f"{r['algorithm']}_L{r['L']}_r{r['steps']}_{r['variant']}"
        row={k:r[k] for k in ('L','tau','steps','algorithm','variant','modes_in_fermionic_cone')}
        row.update({k:r['schedule'][k] for k in ('two_qubit_gates','one_qubit_gates','depth_layers','idle_qubit_layers','seconds_per_attempt')})
        row.update({'chi':r['observable_exposure']['chi'],'full_circuit_chi':r['full_circuit_chi'],
                    'parity_acceptance':r['postselection']['acceptance'],
                    'postselection_bias_upper_bound':r['postselection']['conditional_bias_upper_bound'],
                    'jit_chi':boundary[case]['chi'],'jit_parity_acceptance':boundary[case]['postselection_acceptance']})
        row.update({f'log10_Q_idle_{ratio}':data['log10_Q_min'] for ratio,data in r['costs_by_idle_ratio'].items()})
        rows.append(row)
    with (OUT/'resource_summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (OUT/'plot_runtime.json').write_text(json.dumps({'python':sys.version,'executable':sys.executable,
        'packages':{name:importlib.metadata.version(name) for name in ('numpy','matplotlib')}},indent=2)+'\n')
    print('Wrote four figure pairs and resource_summary.csv')


if __name__=='__main__':main()
