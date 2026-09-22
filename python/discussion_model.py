"""Fast analytical discussion curves with explicit, editable scenario knobs.

Experimental circuit counts anchor NISQ; its effective mitigation overhead is
assumed. Step counts transfer a small-patch fit, not an accuracy certificate.
STAR, FT and hybrid use paper-inspired error/delivery laws. PEC references
show the cost of specified channel inversion beside effective mitigation.
"""
import csv
import json
from math import ceil, exp, log, log10, log2, sqrt
from functools import lru_cache
from statistics import median

from model import log_gamma2_depolarizing, minimum_distance, workload
from project_paths import ROOT


def assumptions():
    hardware=json.loads((ROOT/'configs/assumptions.json').read_text())
    discussion=json.loads((ROOT/'configs/discussion_model.json').read_text())
    return {**hardware,**discussion}


def ideal_shots(c,epsilon,rule='gaussian'):
    """Planning count at the declared signal scale, or a Hoeffding bound."""
    if rule=='hoeffding':return ceil(log(2/c['delta'])/(2*(c['stat_fraction']*epsilon)**2))
    if rule!='gaussian':raise ValueError(rule)
    variance=c['signal_reference']*(1-c['signal_reference'])
    return ceil(c['gaussian_z']**2*variance/(c['stat_fraction']*epsilon)**2)


def sampling_cost(qubits,seconds,chi,c,epsilon,shot_rule='gaussian'):
    log_shots=log(ideal_shots(c,epsilon,shot_rule))+chi
    log_copies=max(0.,log_shots+log(seconds/c['week_seconds']))
    if log_copies<30:log_copies=log(ceil(exp(log_copies)))
    return dict(qubits_per_copy=qubits,seconds_per_attempt=seconds,chi=chi,shot_rule=shot_rule,
                log10_total_shots=log_shots/log(10),
                log10_Q_min=(log(qubits)+log_copies)/log(10)
                if seconds<=c['week_seconds'] else None)


def step_count(tau,epsilon,c):
    if c['step_policy']=='fixed_step_size':return max(1,ceil(c['steps_per_tau']*tau))
    if c['step_policy']!='small_patch_fit':raise ValueError(c['step_policy'])
    fraction=c['trotter_budget_fraction']
    if not 0<fraction<=c['trotter_fraction']:raise ValueError('Trotter allocation exceeds the existing reserve')
    return max(1,ceil(sqrt(c['trotter_coefficient']*tau**c['trotter_time_power']/(fraction*epsilon))))


def fit_evidence(c):
    raw=json.loads((ROOT/'results/refinement_v3/small_dynamics.json').read_text())
    points=[dict(patch=v['patch'],steps=v['steps'],error=v['absolute_trotter_error'],
                 coefficient=v['absolute_trotter_error']*v['steps']**2)
            for v in raw['trotter'] if v['patch'] in ('square2','square3') and v['tau']==1 and v['steps']>=16]
    return dict(points=points,median_coefficient=median(v['coefficient'] for v in points),
                chosen_coefficient=c['trotter_coefficient'],time_power_assumed=c['trotter_time_power'],
                scope='Plaquette formula at tau=1 on two small patches; larger lattices, other orderings, and time law are extrapolations')


@lru_cache(maxsize=None)
def compiled_reference(L):
    return json.loads((ROOT/f'results/refinement_v3/parallel_L{L}_observable/result.json').read_text())


def estimate(L,tau,arm,c,epsilon,steps=None):
    if L<2 or L%2 or tau<=0 or epsilon<=0:raise ValueError('Require even L >= 2, tau > 0, epsilon > 0')
    r=step_count(tau,epsilon,c) if steps is None else steps;N=L*L
    if not isinstance(r,int) or r<1:raise ValueError('Require positive integer steps')
    base=dict(L=L,N=N,tau=tau,r=r,arm=arm,absolute_error=epsilon,
              step_policy=c['step_policy'] if steps is None else 'explicit_steps',
              assumed_trotter_error=c['trotter_coefficient']*tau**c['trotter_time_power']/r**2)
    if arm=='nisq':
        # The cubic gate interpolation and linear depth interpolation pass
        # through the published 6x6 / three-step circuit point.
        G=c['nisq_anchor_gates']*(r/c['nisq_anchor_steps'])*(L/c['nisq_anchor_L'])**3
        depth=c['nisq_anchor_depth']*(r/c['nisq_anchor_steps'])*(L/c['nisq_anchor_L'])
        q=ceil(sqrt(2*N+2))**2
        seconds=c['reset_seconds']+depth*c['gate_seconds']+c['measurement_seconds']
        model=c.get('nisq_model','effective');rule='gaussian'
        if model=='effective':
            chi=log(c['nisq_anchor_sampling_overhead'])*(c['p']/c['nisq_anchor_p'])*G/c['nisq_anchor_gates']
        elif model=='pec_two_qubit':
            G=ceil(G);chi=G*log_gamma2_depolarizing(c['p'],2);rule='hoeffding'
        elif model=='compiled_pec':
            anchor=min(L,8);record=compiled_reference(anchor);s=record['schedule'];x=record['observable_exposure']
            volume=(r/4)*(L/anchor)**3
            chi=(x['exposure']['spam']*(L/anchor)**2+(x['chi']-x['exposure']['spam'])*volume)*(c['p']/c['nisq_anchor_p'])
            G=s['two_qubit_gates']*volume;depth=s['depth_layers']*(r/4)*(L/anchor)
            seconds=c['reset_seconds']+depth*c['gate_seconds']+c['measurement_seconds']
            q=s['physical_copy_qubits'] if L<=8 else ceil(sqrt(2*N))**2
            rule='hoeffding'
        else:raise ValueError(model)
        return {**base,**sampling_cost(q,seconds,chi,c,epsilon,rule),
                'algorithm':'native swap-network proxy' if model!='compiled_pec' else 'compiled plaquette proxy',
                'nisq_model':model,'two_qubit_gates':G,'depth':depth,'distance':0}

    lanes=c['parallel_lanes']
    data_tiles=ceil(1.5*(2*N+6))+c['extra_tiles_per_lane']*lanes
    if arm=='star':
        w=workload(L,tau,r,'termwise',c['u_over_j'])
        chi=4*c['star_alpha']*c['p']*w.theta+4*c['star_floor']*w.rotations
        clocks=2+c['star_clocks_per_rotation']*w.rotations/lanes
        def duration(d):
            return clocks*d*c['cycle_seconds']+w.rotations*c['feedforward_seconds']/lanes+c['reset_seconds']+c['measurement_seconds']
        # Retain the existing sufficient Gamma amplification of logical bias.
        solved=minimum_distance(data_tiles,duration,log(epsilon*c['logical_fraction'])-chi/2,c)
        if solved is None:return {**base,'log10_Q_min':None}
        d,seconds,_=solved
        return {**base,**sampling_cost(2*data_tiles*d*d,seconds,chi,c,epsilon),
                'algorithm':'termwise Strang','distance':d,'rotations':w.rotations,
                'theta':w.theta,'lanes':lanes}

    if arm not in ('ft','hybrid'):raise ValueError(arm)
    if c['p']!=c['nisq_anchor_p']:raise ValueError('Cultivation reference is anchored only at p=0.001')
    w=workload(L,tau,r,'plaquette',c['u_over_j'])
    if arm=='ft':
        epsilon_rotation=epsilon*c['synthesis_fraction']/(2*w.rotations)
        T_per_rotation=ceil(3.02*log2(1/epsilon_rotation)+.77)
        nt=w.fixed_T+w.rotations*T_per_rotation;analog=0;chi=0.
    else:
        T_per_rotation=0;nt=w.fixed_T;analog=w.rotations
        chi=4*c['star_alpha']*c['p']*w.theta+4*c['star_floor']*analog
    log_amp=chi/2
    qsource=c['source_error'];cleanup=int(log(nt*qsource)+log_amp>log(epsilon*c['magic_fraction']))
    qout=35*qsource**3 if cleanup else qsource
    if log(nt*qout)+log_amp>log(epsilon*c['magic_fraction']):return {**base,'log10_Q_min':None}
    candidates=[]
    for factories in c['factory_choices']:
        tiles=data_tiles+factories*(1+12*cleanup)
        def duration(d):
            source=c['source_cycles_per_T']*c.get('source_cycle_multiplier',1)+2*d
            service=(15*source+11*d)/(1-15*qsource) if cleanup else source
            supply=ceil(nt/factories)*service
            consume=(c['ft_clocks_per_T']*nt+c['star_clocks_per_rotation']*analog+(2*w.rotations if arm=='ft' else 0))*d/lanes
            return (2*d+max(supply,consume))*c['cycle_seconds']+(nt+analog)*c['feedforward_seconds']/lanes+c['reset_seconds']+c['measurement_seconds']
        solved=minimum_distance(tiles,duration,log(epsilon*c['logical_fraction'])-log_amp,c)
        if solved is None:continue
        d,_,_=solved;d=max(d,c['source_end_distance']);seconds=duration(d)
        q=2*tiles*d*d+factories*c['source_qubits']
        candidates.append({**base,**sampling_cost(q,seconds,chi,c,epsilon),
                           'algorithm':'plaquette Strang','distance':d,'T_count':nt,
                           'theta':w.theta,'analog_rotations':analog,'fixed_T':w.fixed_T,
                           'T_per_rotation':T_per_rotation,'factories':factories,
                           'cleanup_levels':cleanup,'magic_error':qout,'lanes':lanes})
    feasible=[v for v in candidates if v['log10_Q_min'] is not None]
    return min(feasible,key=lambda v:v['log10_Q_min']) if feasible else {**base,'log10_Q_min':None}


def capacity(rows,Q):
    feasible=[v['L'] for v in rows if v['log10_Q_min'] is not None and v['log10_Q_min']<=log10(Q)+1e-12]
    return max(feasible,default=0)


def families(c):
    return [(f'tau_{tau:g}',tau) for tau in c['times']]+[('tau_proportional',None)]


def sweep(c):
    rows=[]
    for epsilon in (.001,.01):
        for family,fixed_tau in families(c):
            for arm in ('nisq','star','ft','hybrid'):
                scenarios=dict(c['scenarios'])
                if arm=='nisq':scenarios.update(pec_two_qubit={'nisq_model':'pec_two_qubit'},compiled_pec={'nisq_model':'compiled_pec'})
                for scenario,changes in scenarios.items():
                    for L in range(2,c['max_L']+1,2):
                        tau=fixed_tau if fixed_tau is not None else c['time_family_scale']*L
                        rows.append({**estimate(L,tau,arm,{**c,**changes},epsilon),'scenario':scenario,'family':family})
    return rows


def plot(rows,c,epsilon=.001):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from plot_capacity import COLORS
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'svg.fonttype':'none'})
    panels=families(c);nrows=ceil(len(panels)/2)
    fig,axes=plt.subplots(nrows,2,figsize=(13.6,3.4*nrows+2.0),sharey=True,squeeze=False)
    labels={'nisq':'NISQ · assumed effective mitigation','star':'STAR · analog rotations',
            'ft':'FTQC · cultivated T states','hybrid':'Hybrid · analog + fixed T states'}
    ymax=4*ceil(max(capacity([v for v in rows if v['family']==family and v['arm']==arm and v['scenario']==scenario and v['absolute_error']==epsilon],1e9)
                   for family,_ in panels for arm in labels for scenario in c['scenarios'])/4)
    ymax=max(ymax,4)
    for ax,(family,tau) in zip(axes.flat,panels):
        for arm in labels:
            subset=[v for v in rows if v['arm']==arm and v['family']==family and v['absolute_error']==epsilon]
            scenarios=list(c['scenarios'])+(['pec_two_qubit','compiled_pec'] if arm=='nisq' else [])
            thresholds=[10**v['log10_Q_min'] for v in subset if v['log10_Q_min'] is not None and 1<=v['log10_Q_min']<=9]
            Q=np.array(sorted(set(np.geomspace(10,1e9,220).tolist()+thresholds)))
            values={scenario:np.array([capacity([v for v in subset if v['scenario']==scenario],q) for q in Q],dtype=float)
                    for scenario in scenarios}
            lower=np.minimum.reduce(list(values.values()));upper=np.maximum.reduce(list(values.values()))
            ax.fill_between(Q,lower,upper,step='post',color=COLORS[arm],alpha=.10)
            center=values['central'].copy();center[center==0]=np.nan
            ax.step(Q,center,where='post',color=COLORS[arm],lw=2.2,
                    ls='--' if arm=='hybrid' else '-',label=labels[arm])
            if arm=='nisq':
                for scenario,style,label in [('pec_two_qubit',':','NISQ PEC · two-qubit noise'),
                                             ('compiled_pec','-.','NISQ PEC · compiled-noise proxy')]:
                    curve=values[scenario].copy();curve[curve==0]=np.nan
                    ax.step(Q,curve,where='post',color=COLORS[arm],lw=1.35,ls=style,label=label)
        ax.set_xscale('log');ax.set_xlim(10,1e9);ax.set_ylim(0,ymax+1)
        ax.set_xticks([1e2,1e4,1e6,1e8]);ax.set_yticks(range(0,ymax+1,4 if ymax<=36 else 8))
        ax.set_xlabel('Available physical qubits, Q')
        title=f'τ = {tau:g}  ·  r = {step_count(tau,epsilon,c)}' if tau is not None else f'τ = {c["time_family_scale"]:g} L  ·  size-dependent r'
        ax.set_title(title,loc='left',fontsize=12);ax.grid(alpha=.18)
        if ymax>=c['max_L']:
            ax.axhspan(c['max_L']-.4,c['max_L']+1,color='#eeeeee',zorder=0)
            ax.text(.98,.975,f'L = {c["max_L"]}: sweep limit',ha='right',va='top',transform=ax.transAxes,fontsize=8,color='#666666')
    for ax in list(axes.flat)[len(panels):]:ax.set_visible(False)
    for ax in axes[:,0]:ax.set_ylabel('Largest modeled side L (L × L sites)')
    handles,legend=axes[0,0].get_legend_handles_labels()
    order=[0,3,4,5,1,2]
    fig.legend([handles[i] for i in order],[legend[i] for i in order],loc='lower center',
               bbox_to_anchor=(.5,.064),ncol=3,frameon=False,fontsize=9)
    fig.suptitle('Discussion model: step refinement, mitigation, and hybrid protection',fontsize=16,x=.07,ha='left',y=.987)
    time_power={2:'²',3:'³'}.get(c['trotter_time_power'],'^'+str(c['trotter_time_power']))
    rule=(f'assumed Trotter error = {c["trotter_coefficient"]:g} τ{time_power}/r², budget = {c["trotter_budget_fraction"]:g} ε'
          if c['step_policy']=='small_patch_fit' else f'r = ceil({c["steps_per_tau"]:g} τ)')
    fig.text(.07,.944,f'One week · p = 10⁻³ · ε = {epsilon:g} · {rule}',fontsize=10,color='#555555')
    fig.text(.5,.012,'Small-patch-informed scenario, not a large-lattice accuracy certificate. Bands are scenario ranges, including PEC references.\n'
             'Effective NISQ bias is assumed affordable. PEC uses Hoeffding; other curves use Gaussian planning counts. Full-lattice, even-L estimates.',
             ha='center',fontsize=8.2,color='#555555')
    fig.subplots_adjust(left=.07,right=.985,bottom=.205,top=.875,wspace=.10,hspace=.34)
    name='discussion_capacity' if epsilon==.001 else 'discussion_capacity_loose'
    for ext in ('pdf','png','svg'):
        path=ROOT/f'figures/{name}.{ext}';fig.savefig(path,dpi=200)
        if ext=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


def run():
    c=assumptions();rows=sweep(c);out=ROOT/'results/discussion_v2';out.mkdir(exist_ok=True)
    (out/'assumptions.json').write_text(json.dumps(c,indent=2)+'\n')
    (out/'costs.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    (out/'trotter_fit.json').write_text(json.dumps(fit_evidence(c),indent=2)+'\n')
    summaries=[]
    for epsilon in (.001,.01):
        for family,tau in families(c):
            for Q in c['qubit_budgets']:
                values={arm:capacity([v for v in rows if v['arm']==arm and v['scenario']=='central' and v['family']==family and v['absolute_error']==epsilon],Q)
                        for arm in ('nisq','star','ft','hybrid')}
                for scenario in ('pec_two_qubit','compiled_pec'):
                    values[scenario]=capacity([v for v in rows if v['arm']=='nisq' and v['scenario']==scenario and v['family']==family and v['absolute_error']==epsilon],Q)
                summaries.append(dict(absolute_error=epsilon,family=family,Q=Q,**values))
    with (out/'capacity.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
    sensitivity=[]
    for epsilon in (.001,.01):
        for L in (4,6):
            for arm in ('nisq','star','ft','hybrid'):
                for steps in (4,step_count(1,epsilon,c)):
                    v=estimate(L,1,arm,c,epsilon,steps=steps)
                    sensitivity.append({key:v.get(key) for key in ('L','arm','absolute_error','r','chi','T_count','log10_total_shots','qubits_per_copy','seconds_per_attempt','log10_Q_min')})
    with (out/'step_sensitivity.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(sensitivity[0]));writer.writeheader();writer.writerows(sensitivity)
    for epsilon in (.001,.01):plot(rows,c,epsilon)
    summary=['# Revised back-of-the-envelope+ discussion example','',
             'One week, p = 0.001, 10 ns gates, 240 ns syndrome cycles. Four architectures; effective mitigation and PEC are shown separately.',
             f'The chosen step policy is {c["step_policy"]}; the transferred fit is {c["trotter_coefficient"]:g} tau^{c["trotter_time_power"]:g}/r^2, with Trotter allowance {c["trotter_budget_fraction"]:g} epsilon.',
             f'At tau = 1 this gives r = {step_count(1,.001,c)} for epsilon = 0.001 and r = {step_count(1,.01,c)} for epsilon = 0.01.',
             'Larger sizes, other product-formula orderings, and the time law are assumptions, not calibrated accuracy.',
             f'With coefficient 0.22, a 25% allocation gives 30/10 steps (strict/loose); a 12.5% allocation gives 42/14.','',
             '## Largest lattice at tau = 1, absolute error 0.001','',
             '| Physical qubits | NISQ effective | NISQ 2q PEC | NISQ compiled PEC | STAR | Hybrid | FTQC |',
             '|---:|---:|---:|---:|---:|---:|---:|']
    def size(L):return f'{L} × {L}' if L else '—'
    for v in summaries:
        if v['absolute_error']==.001 and v['family']=='tau_1' and v['Q'] in (100,10000,100000,1000000,10000000):
            summary.append('| '+f'{v["Q"]:,}'+' | '+' | '.join(size(v[key]) for key in ('nisq','pec_two_qubit','compiled_pec','star','hybrid','ft'))+' |')
    summary+=['','## Central area required at tau = 1','',
              '| Lattice | Absolute precision | Steps | NISQ effective | STAR | Hybrid | FTQC |',
              '|---:|---:|---:|---:|---:|---:|---:|']
    def area(v):
        value=v['log10_Q_min']
        return 'outside model' if value is None else (f'{round(10**value):,}' if value<8 else f'~10^{value:.1f}')
    for L in (4,6):
        for epsilon in (.001,.01):
            selected={arm:next(v for v in rows if v['L']==L and v['family']=='tau_1' and v['arm']==arm and v['absolute_error']==epsilon and v['scenario']=='central') for arm in ('nisq','star','hybrid','ft')}
            summary.append(f'| {size(L)} | {epsilon} | {step_count(1,epsilon,c)} | '+' | '.join(area(selected[arm]) for arm in ('nisq','star','hybrid','ft'))+' |')
    summary+=['','## What step refinement changes at 4 × 4, tau = 1, strict precision','',
              '| Arm | Q at r = 4 | Q at the chosen refined r |','|---|---:|---:|']
    for arm in ('nisq','star','hybrid','ft'):
        pair=[v for v in sensitivity if v['L']==4 and v['arm']==arm and v['absolute_error']==.001]
        summary.append(f'| {arm} | {area(pair[0])} | {area(pair[1])} |')
    summary+=['','The FT jump includes crossing the chosen magic-state cleanup threshold; it is not purely a linear step-count penalty.',
              'The hybrid does not automatically beat termwise STAR, which already avoids fixed T gates.',
              'The effective NISQ line still assumes affordable residual bias and is not a measured performance claim.',
              'The two-qubit PEC reference uses exact inverse-channel overhead on the gate proxy. The compiled PEC edge transfers saved full-noise exposure and timing, linearly in step count.',
              'Both PEC references use Hoeffding shot counts; central effective, STAR, hybrid, and FT curves use Gaussian planning counts.',
              'The hybrid injects variable angles and supplies high-fidelity fixed T states. It shares the FT factory and applies the same Gamma-amplified logical/magic budget as STAR.','',
              f'The fourth panel follows tau = {c["time_family_scale"]:g} L. With the chosen small-patch-fit time power, r grows as L^({c["trotter_time_power"]/2:g}), so the NISQ exponent grows as L^({3+c["trotter_time_power"]/2:g}).',
              'The L^4 law applies instead when the step count grows linearly in time; choose trotter_time_power = 2 for that sensitivity.','',
              'Details: [discussion guide](../../docs/DISCUSSION_MODEL.md).',
              'Figures: [strict precision](../../figures/discussion_capacity.pdf), [0.01 sensitivity](../../figures/discussion_capacity_loose.pdf).',
              'The [original fixed-step example](../discussion_v1/summary.md) is archived.','']
    (out/'summary.md').write_text('\n'.join(summary))
    for v in summaries:
        if v['absolute_error']==.001 and v['family']=='tau_1':print(v)
    print('Wrote',len(rows),'analytical scenario points; no compiler or dynamics jobs run.')


if __name__=='__main__':run()
