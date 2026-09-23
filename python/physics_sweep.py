"""Frozen-baseline p/precision sweep and bounded algorithm alternatives.

Run normally to calculate; --plot rebuilds presentation from saved numbers;
--check audits ledgers, baseline replay, and source assumptions.
"""
import argparse
import csv
import hashlib
import json
from math import ceil, exp, isclose, log, log10

from budget_optimization import algorithm_for, coefficient, optimize
from discussion_model import assumptions, capacity
from dynamics_followup_analysis import extended_coefficient
from project_paths import ROOT
from resource_extensions import source_error

OUT = ROOT/'results/physics_sweep_v1'
ARMS = ('nisq', 'star', 'hybrid', 'ft')


def read(path):
    return json.loads(path.read_text())


def settings():
    return read(ROOT/'configs/physics_sweep.json')


def write(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT/name
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def csv_write(name, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def configuration(p, scenario='held'):
    config = settings()
    choices = config['source_scenarios'][scenario]
    c = {**assumptions(), 'p':p, 'factory_choices':config['factory_choices'], **choices}
    if choices['floor_policy'] == 'linear':
        c['star_floor'] *= p/c['nisq_anchor_p']
    return c


def planning_coefficient(L, tau, algorithm):
    if L == 2:
        from budget_optimization import settings as budgets
        return coefficient(tau, algorithm, {**budgets(), 'coefficient_patch':'2x2'})
    return extended_coefficient(tau, algorithm, read(ROOT/'results/dynamics_followup_v1/measurements.json'))


def solve(L, tau, arm, epsilon, c, **labels):
    algorithm = c.get('workload_algorithm', algorithm_for(arm))
    C = planning_coefficient(L, tau, algorithm)
    row = optimize(L, tau, arm, c, epsilon, C_override=C)
    row.update(p=c['p'], **labels)
    row['model_options'] = {k:c[k] for k in (
        'source_p_policy','source_cycle_multiplier','star_floor','star_law','local_encoding',
        'local_gates_per_hopping_pauli','local_routing_multiplier','local_layout_factor',
        'parallel_lanes','rotation_synthesis','hwp_batch_size','workload_algorithm') if k in c}
    row['coefficient_patch'] = '2x2' if L == 2 else 'envelope of 3x3, 4x3, 4x4'
    if arm in ('ft','hybrid'):
        row['input_source_error'] = source_error(c)
    return row


def run_sweep():
    config = settings(); rows = []
    for p in config['physical_errors']:
        for epsilon in config['precisions']:
            for L in config['sizes']:
                for arm in ARMS:
                    scenarios = ('held','scaled') if arm == 'star' else config['source_scenarios']
                    if arm == 'nisq':
                        for model in ('effective','pec_two_qubit','compiled_pec'):
                            c = {**configuration(p), 'nisq_model':model}
                            rows.append(solve(L,config['tau'],arm,epsilon,c,scenario=model))
                    else:
                        for scenario in scenarios:
                            rows.append(solve(L,config['tau'],arm,epsilon,configuration(p,scenario),scenario=scenario))
                        if arm == 'star':
                            c = {**configuration(p), 'star_law':'original_count'}
                            rows.append(solve(L,config['tau'],arm,epsilon,c,scenario='original_count'))
            write('costs.json',rows)
            print('physics sweep',p,epsilon,len(rows),'rows',flush=True)
    baseline = read(ROOT/'results/dynamics_followup_v1/costs.json')
    frozen = [v for v in baseline if v['patch_policy']=='added_patches' and v['tau']==config['tau']]
    write('frozen_baseline.json',dict(source='results/dynamics_followup_v1/costs.json',
          source_sha256=hashlib.sha256((ROOT/'results/dynamics_followup_v1/costs.json').read_bytes()).hexdigest(), rows=frozen))
    write('assumptions.json',dict(hardware=assumptions(),sweep=config,
        scope='Conditional source and injection scenarios; no calibration of physical gadgets at lower p'))
    return rows


def run_algorithms():
    config = settings(); rows = []
    for epsilon in config['precisions']:
        for L in config['algorithm_sizes']:
            for model in ('effective','pec_two_qubit'):
                c = {**configuration(.001),'nisq_model':model}
                rows.append(solve(L,1,'nisq',epsilon,c,variant='JW reference',scenario=model))
                for name,changes in [
                    ('compact parity gadgets',{'local_gates_per_hopping_pauli':4}),
                    ('compact fused hypothesis',{'local_gates_per_hopping_pauli':2}),
                    ('compact routing x2',{'local_gates_per_hopping_pauli':4,'local_routing_multiplier':2,'local_layout_factor':2})]:
                    rows.append(solve(L,1,'nisq',epsilon,{**c,'local_encoding':True,'workload_algorithm':'termwise',**changes},
                                      variant=name,scenario=model))
            for lanes in (4,16):
                rows.append(solve(L,1,'star',epsilon,{**configuration(.001),'local_encoding':True,'parallel_lanes':lanes},
                                  variant=f'compact {lanes} lanes',scenario='angle'))
            rows.append(solve(L,1,'star',epsilon,configuration(.001),variant='JW reference',scenario='angle'))
            rows.append(solve(L,1,'ft',epsilon,configuration(.001),variant='individual synthesis',scenario='ft'))
            for batch in config['hwp_batches']:
                if batch > L*L:continue
                c = {**configuration(.001),'rotation_synthesis':'hamming_weight','hwp_batch_size':batch}
                rows.append(solve(L,1,'ft',epsilon,c,variant=f'HWP batch {batch}',scenario='ft'))
            print('algorithm alternatives',epsilon,L,flush=True)
    write('algorithms.json',rows)
    return rows


def capacity_table(rows):
    config = settings(); output = []
    for p in config['physical_errors']:
        for epsilon in config['precisions']:
            for Q in config['qubit_budgets']:
                row = dict(p=p,absolute_error=epsilon,Q=Q)
                for arm,scenario,key in [('nisq','effective','nisq_effective'),('nisq','pec_two_qubit','nisq_pec'),
                        ('nisq','compiled_pec','nisq_compiled_pec'),('star','held','star'),('star','original_count','star_count'),
                        ('hybrid','held','hybrid'),('ft','held','ft')]:
                    group = [v for v in rows if (v['p'],v['absolute_error'],v['arm'],v['scenario'])==(p,epsilon,arm,scenario)]
                    row[key] = capacity(group,Q)
                best = max(row[k] for k in ('nisq_effective','star','hybrid','ft'))
                row['capacity_leaders'] = ','.join(k for k in ('nisq_effective','star','hybrid','ft') if row[k]==best) if best else 'none'
                row['sweep_ceiling_reached'] = best == max(config['sizes'])
                output.append(row)
    return output


def plot(rows, algorithms):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np
    from plot_capacity import COLORS
    config=settings()
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,3,figsize=(13.4,8.0),sharex=True,sharey=True)
    for i,epsilon in enumerate(config['precisions']):
        for j,p in enumerate(config['physical_errors']):
            ax=axes[i,j]; group=[v for v in rows if v['p']==p and v['absolute_error']==epsilon]
            thresholds=[10**v['log10_Q_min'] for v in group if v['log10_Q_min'] is not None and 1<=v['log10_Q_min']<=9]
            Q=sorted(set(np.geomspace(10,1e9,220))|set(thresholds))
            for arm in ARMS:
                scenarios=['effective'] if arm=='nisq' else (['held','scaled'] if arm=='star' else list(config['source_scenarios']))
                values={s:np.array([capacity([v for v in group if v['arm']==arm and v['scenario']==s],q) for q in Q],float) for s in scenarios}
                if arm!='nisq':
                    ax.fill_between(Q,np.minimum.reduce(list(values.values())),np.maximum.reduce(list(values.values())),step='post',alpha=.13,color=COLORS[arm])
                y=values[scenarios[0]].copy();y[y==0]=np.nan
                ax.step(Q,y,where='post',color=COLORS[arm],lw=2,ls='--' if arm=='hybrid' else '-')
            for arm,scenario,color,style in [('nisq','pec_two_qubit',COLORS['nisq'],':'),('nisq','compiled_pec','#777777',':'),
                    ('star','original_count',COLORS['star'],'-.')]:
                y=np.array([capacity([v for v in group if v['arm']==arm and v['scenario']==scenario],q) for q in Q],float)
                y[y==0]=np.nan;ax.step(Q,y,where='post',color=color,ls=style,lw=1.2)
            ax.set_xscale('log');ax.set_xlim(10,1e9);ax.set_ylim(0,50);ax.set_yticks(range(0,49,8));ax.grid(alpha=.2)
            ax.set_title(f'p = {p:g}  ·  ε = {epsilon:g}',loc='left')
            if i==1:ax.set_xlabel('Available physical qubits, Q')
            if j==0:ax.set_ylabel('Largest modeled side L')
    legend=[(COLORS['nisq'],'-','NISQ: assumed effective mitigation'),(COLORS['nisq'],':','NISQ: two-qubit PEC'),
            ('#777777',':','NISQ: compiled PEC proxy'),(COLORS['star'],'-','STAR: angle + floor'),
            (COLORS['star'],'-.','STAR: original count law'),(COLORS['hybrid'],'--','Hybrid'),(COLORS['ft'],'-','FTQC')]
    fig.legend([Line2D([],[],color=c,ls=s,lw=2) for c,s,_ in legend],[t for _,_,t in legend],
               loc='lower center',bbox_to_anchor=(.5,.065),ncol=3,frameon=False,fontsize=9)
    fig.suptitle('Physical error and precision change the architecture comparison',x=.065,ha='left')
    fig.text(.5,.018,'τ = 1 · one week · fixed clocks and finite-patch Trotter envelope · discrete sizes through L = 48\n'
             'Bands vary source/floor assumptions, not confidence. Effective NISQ bias and lower-p source performance are assumed.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.18,1,.95))
    for ext in ('pdf','png'):fig.savefig(ROOT/f'figures/physical_error_capacity.{ext}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.7))
    for ax,epsilon in zip(axes,config['precisions']):
        for variant in ['individual synthesis','HWP batch 8','HWP batch 16','best declared menu']:
            group=[v for v in algorithms if v['arm']=='ft' and v['absolute_error']==epsilon and v['variant']==variant]
            if variant=='best declared menu':
                group=[min([v for v in algorithms if v['arm']=='ft' and v['absolute_error']==epsilon and v['L']==L],
                           key=lambda v:v['log10_Q_min']) for L in config['algorithm_sizes']]
            Q=sorted(set(np.geomspace(1e4,1e9,200))|{10**v['log10_Q_min'] for v in group if 4<=v['log10_Q_min']<=9})
            y=np.array([capacity(group,q) for q in Q],float);y[y==0]=np.nan
            ax.step(Q,y,where='post',label=variant)
        ax.set_xscale('log');ax.set_xlim(1e4,1e9);ax.set_ylim(0,17);ax.set_yticks([0,4,8,16])
        ax.set_xlabel('Available physical qubits, Q');ax.set_title(f'ε = {epsilon:g}',loc='left');ax.grid(alpha=.2)
    axes[0].set_ylabel('Largest tested side L');axes[1].legend(frameon=False,fontsize=9)
    fig.suptitle('Repeated-angle synthesis: T savings must pay for workspace',x=.08,ha='left')
    fig.text(.5,.02,'τ = 1 · p = 10⁻³ · only L = 4, 8, 16 tested · same plaquette accuracy envelope\nBatch workspace and arithmetic charged; steps and error allocations reoptimized.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.07,1,.94))
    for ext in ('pdf','png'):fig.savefig(ROOT/f'figures/hamming_weight_comparison.{ext}',dpi=180)
    plt.close(fig)


def validate(rows, algorithms):
    config=settings()
    assert len(rows)==len(config['physical_errors'])*len(config['precisions'])*len(config['sizes'])*12
    assert len({(v['p'],v['absolute_error'],v['L'],v['arm'],v['scenario']) for v in rows})==len(rows)
    frozen=read(OUT/'frozen_baseline.json');baseline_checks=0
    assert frozen['source_sha256']==hashlib.sha256((ROOT/frozen['source']).read_bytes()).hexdigest()
    for old in frozen['rows']:
        scenario=old['scenario'] if old['arm']=='nisq' else 'held'
        match=[v for v in rows if (v['p'],v['absolute_error'],v['L'],v['arm'],v['scenario'])==(.001,old['absolute_error'],old['L'],old['arm'],scenario)]
        assert len(match)==1
        new=match[0]
        assert new['r']==old['r'] and isclose(new['log10_Q_min'],old['log10_Q_min'],abs_tol=1e-12)
        baseline_checks+=1
    for v in rows+algorithms:
        if v['log10_Q_min'] is None:continue
        b=v['budget_fractions'];assert isclose(sum(b.values()),1,abs_tol=1e-12) and min(b.values())>=0
        assert isclose(v['assumed_trotter_error'],v['absolute_error']*b['trotter'],rel_tol=1e-12)
        c={**configuration(v['p']),**v['model_options']}
        if v['arm']!='nisq':
            tiles=v.get('protected_tiles',v['data_tiles'])
            pL=c['surface_A']*(c['p']/c['surface_threshold'])**((v['distance']+1)/2)
            faults=tiles*ceil(v['seconds_per_attempt']/c['cycle_seconds'])*pL
            assert log(faults)+v['chi']/2<=log(v['absolute_error']*b['logical'])+1e-10
            assert isclose(v['chi'],v['chi_angle']+v['chi_floor']+v.get('chi_count',0),abs_tol=1e-9)
        if v['arm'] in ('ft','hybrid'):
            assert log(v['T_count']*v['magic_error'])+v['chi']/2<=log(v['absolute_error']*b['magic'])+1e-10
            assert isclose(v['input_source_error'],source_error(c),rel_tol=1e-12)
            # Cleanup circuitry is included in the protected tile-time fault ledger.
            assert v['protected_tiles']>=v['data_tiles']+v['factories']*(1+12*v['cleanup_levels'])
    ordering=read(OUT/'ordering.json')
    assert len(ordering['rows'])==8
    for name,digest in ordering['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    for v in ordering['rows']:
        assert v['relative_tail_spread']<.05
        for point in v['measurements']:
            assert point['norm_error']<1e-10
            assert isclose(point['coefficient'],abs(point['signed_error'])*point['r']**2,rel_tol=1e-12)
    record=dict(status='PASS: frozen baseline, error ledgers, source domains, and row census',
                sweep_rows=len(rows),algorithm_rows=len(algorithms),baseline_replays=baseline_checks)
    write('validation.json',record)
    return record


def summaries(rows, algorithms):
    table=capacity_table(rows);write('capacity.json',table);csv_write('capacity.csv',table)
    fields=['p','absolute_error','L','arm','scenario','variant','r','distance','chi','chi_angle','chi_floor','chi_count',
            'two_qubit_gates','rotations','theta','T_count','synthesized_rotations','hwp_toffolis','hwp_workspace',
            'factories','cleanup_levels','qubits_per_copy','seconds_per_attempt','log10_Q_min']
    csv_write('algorithms.csv',[{k:v.get(k) for k in fields} for v in algorithms])
    csv_write('operating_points.csv',[{k:v.get(k) for k in fields} for v in rows if v['L']==4])
    text=['# Physical error, precision, and algorithm alternatives','',
          'The p = 0.001, held-source rows reproduce the frozen added-patch baseline. The main figure leads with 0.01 absolute error; 0.001 remains equally visible.',
          '[Capacity figure](../../figures/physical_error_capacity.pdf) · [Hamming-weight comparison](../../figures/hamming_weight_comparison.pdf) · [Full manuscript](../../tex/main.pdf)','',
          'The source footprint (1024 qubits), nominal preparation cycles (100), and minimum handoff distance (15) stay fixed. Held scenarios retain source error 2e-9 and STAR floor 1e-5. The scaled scenario assumes source error proportional to p³ and floor proportional to p. Slow supply multiplies the source preparation cycles by four. These are scenarios, not confidence bounds or new gadget calibrations.',
          'The original STAR count law uses the exact dephasing inverse with q = 4p/15 per rotation. Its delivery model is held fixed to isolate the injection law. It is not a full reconstruction of the original protocol.',
          'Cleanup has ideal input-error suppression 35 q³; faults in its protected circuitry and stored states are charged separately to the logical tile-time budget. No zero-fault-per-campaign requirement is imposed.','',
          '## Capacity at Q = one million, tau = 1','',
          '| p | absolute error | NISQ effective | NISQ 2q PEC | STAR angle | STAR count | Hybrid | FTQC |',
          '|---:|---:|---:|---:|---:|---:|---:|---:|']
    for v in table:
        if v['Q']==1000000:
            text.append('| '+' | '.join(str(v[k]) for k in ('p','absolute_error','nisq_effective','nisq_pec','star','star_count','hybrid','ft'))+' |')
    text+=['','Entries are lattice sides L. A value of 48 reaches the sampled size ceiling, not a hardware limit. Full discrete staircases and ties are retained in capacity.csv.',
           '','## FT repeated-angle synthesis','',
           '| L | epsilon | Individual T / shot | Selected HWP T / shot | HWP batch | Live ancillas | Q ratio (individual / HWP) |',
           '|---:|---:|---:|---:|---:|---:|---:|']
    for eps in settings()['precisions']:
        for L in settings()['algorithm_sizes']:
            group=[v for v in algorithms if v['arm']=='ft' and v['absolute_error']==eps and v['L']==L]
            a=next(v for v in group if v['variant']=='individual synthesis')
            b=min(group,key=lambda v:v['log10_Q_min'])
            text.append(f'| {L} | {eps} | {a["T_count"]:,} | {b["T_count"]:,} | {b.get("hwp_batch_size",1)} | {b.get("hwp_workspace",0)} | {10**(a["log10_Q_min"]-b["log10_Q_min"]):.2f} |')
    text+=['','HWP batches only equal-angle rotations within the same stage. A b-bit batch uses b-popcount(b) clean ancillas/Toffolis and floor(log2 b)+1 synthesized rotations. Ancillas are reused between batches but charged for the entire shot. T counts include four T per measurement-assisted Toffoli; runtime additionally charges two logical Clifford clocks per Toffoli. Routing and multiple delivery ports remain the reference layout assumptions.',
           'At strict 4x4 neither selected circuit uses cleanup. HWP reduces the magic allocation from 38.5% to 16.0% and increases the statistical allocation from 29.8% to 58.0%. The reduced shot count makes the area gain larger than the T-count reduction.',
           '','## Local encoding and ordering','',
           'Compact encoding uses a 3N-qubit register estimate. The parity-gadget scenario charges four CNOTs per weight-three hopping Pauli and two per onsite Pauli, plus 8N preparation gates. A two-CNOT hopping scenario assumes fusion; a routing scenario doubles gate and placement costs. Native-grid layout, preparation, and 50% gate utilization are uncompiled assumptions. Both the effective-mitigation and two-qubit PEC costs are shown in algorithms.csv.',
           'The compact scenario uses termwise evolution and the existing termwise accuracy envelope. Thus the comparison with the plaquette-calibrated native proxy changes algorithm as well as mapping; counts are separately exposed. Compact STAR keeps that same envelope and compares four with sixteen assumed analog lanes. No parallel speedup is inferred from encoding locality alone.',
           'The separate ordering.json experiment changes only bond ordering, then changes onsite grouping as a second variant, on 2x2 and 3x3 at tau=1. These diagnostic coefficients do not replace the frozen planning envelope.',
           '','The workshop-model comparison is in [../workshop_reconciliation_v1/summary.md](../workshop_reconciliation_v1/summary.md).',
           'Reproduce with `make physics`; audit with `python python/physics_sweep.py --check`.','']
    (OUT/'summary.md').write_text('\n'.join(text))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.plot or args.check:
        rows=read(OUT/'costs.json');algorithms=read(OUT/'algorithms.json')
    else:
        rows=run_sweep();algorithms=run_algorithms()
    print(validate(rows,algorithms),flush=True)
    if not args.check:
        summaries(rows,algorithms);plot(rows,algorithms)


if __name__=='__main__':main()
