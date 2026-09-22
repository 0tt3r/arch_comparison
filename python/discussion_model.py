"""Fast analytical discussion curves with explicit, editable scenario knobs.

Experimental circuit counts anchor NISQ; its effective mitigation overhead is
assumed. STAR and FT use paper-inspired error/delivery laws. All three assume
their allocated residual bias can be met; no Trotter or training validation.
"""
import csv
import json
from math import ceil, exp, log, log10, log2, sqrt

from model import minimum_distance, workload
from project_paths import ROOT


def assumptions():
    hardware=json.loads((ROOT/'configs/assumptions.json').read_text())
    discussion=json.loads((ROOT/'configs/discussion_model.json').read_text())
    return {**hardware,**discussion}


def ideal_shots(c,epsilon):
    """Gaussian planning count for a Bernoulli doublon at the declared scale."""
    variance=c['signal_reference']*(1-c['signal_reference'])
    return ceil(c['gaussian_z']**2*variance/(c['stat_fraction']*epsilon)**2)


def sampling_cost(qubits,seconds,chi,c,epsilon):
    log_shots=log(ideal_shots(c,epsilon))+chi
    log_copies=max(0.,log_shots+log(seconds/c['week_seconds']))
    if log_copies<30:log_copies=log(ceil(exp(log_copies)))
    return dict(qubits_per_copy=qubits,seconds_per_attempt=seconds,chi=chi,
                log10_total_shots=log_shots/log(10),
                log10_Q_min=(log(qubits)+log_copies)/log(10)
                if seconds<=c['week_seconds'] else None)


def estimate(L,tau,arm,c,epsilon):
    if L<2 or L%2 or tau<=0 or epsilon<=0:raise ValueError('Require even L >= 2, tau > 0, epsilon > 0')
    r=ceil(c['steps_per_tau']*tau);N=L*L
    base=dict(L=L,N=N,tau=tau,r=r,arm=arm,absolute_error=epsilon)
    if arm=='nisq':
        # The cubic gate interpolation and linear depth interpolation pass
        # through the published 6x6 / three-step circuit point.
        G=c['nisq_anchor_gates']*(r/c['nisq_anchor_steps'])*(L/c['nisq_anchor_L'])**3
        depth=c['nisq_anchor_depth']*(r/c['nisq_anchor_steps'])*(L/c['nisq_anchor_L'])
        chi=log(c['nisq_anchor_sampling_overhead'])*(c['p']/c['nisq_anchor_p'])*G/c['nisq_anchor_gates']
        q=ceil(sqrt(2*N+2))**2
        seconds=c['reset_seconds']+depth*c['gate_seconds']+c['measurement_seconds']
        return {**base,**sampling_cost(q,seconds,chi,c,epsilon),
                'algorithm':'native swap-network proxy','two_qubit_gates':G,'depth':depth,'distance':0}

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

    if arm!='ft':raise ValueError(arm)
    if c['p']!=c['nisq_anchor_p']:raise ValueError('Cultivation reference is anchored only at p=0.001')
    w=workload(L,tau,r,'plaquette',c['u_over_j'])
    epsilon_rotation=epsilon*c['synthesis_fraction']/(2*w.rotations)
    T_per_rotation=ceil(3.02*log2(1/epsilon_rotation)+.77)
    nt=w.fixed_T+w.rotations*T_per_rotation
    qsource=c['source_error'];cleanup=int(nt*qsource>epsilon*c['magic_fraction'])
    qout=35*qsource**3 if cleanup else qsource
    if nt*qout>epsilon*c['magic_fraction']:return {**base,'log10_Q_min':None}
    candidates=[]
    for factories in c['factory_choices']:
        tiles=data_tiles+factories*(1+12*cleanup)
        def duration(d):
            source=c['source_cycles_per_T']*c.get('source_cycle_multiplier',1)+2*d
            service=(15*source+11*d)/(1-15*qsource) if cleanup else source
            supply=ceil(nt/factories)*service
            consume=(c['ft_clocks_per_T']*nt+2*w.rotations)*d/lanes
            return (2*d+max(supply,consume))*c['cycle_seconds']+nt*c['feedforward_seconds']/lanes+c['reset_seconds']+c['measurement_seconds']
        solved=minimum_distance(tiles,duration,log(epsilon*c['logical_fraction']),c)
        if solved is None:continue
        d,_,_=solved;d=max(d,c['source_end_distance']);seconds=duration(d)
        q=2*tiles*d*d+factories*c['source_qubits']
        candidates.append({**base,**sampling_cost(q,seconds,0.,c,epsilon),
                           'algorithm':'plaquette Strang','distance':d,'T_count':nt,
                           'T_per_rotation':T_per_rotation,'factories':factories,
                           'cleanup_levels':cleanup,'magic_error':qout,'lanes':lanes})
    feasible=[v for v in candidates if v['log10_Q_min'] is not None]
    return min(feasible,key=lambda v:v['log10_Q_min']) if feasible else {**base,'log10_Q_min':None}


def capacity(rows,Q):
    feasible=[v['L'] for v in rows if v['log10_Q_min'] is not None and v['log10_Q_min']<=log10(Q)+1e-12]
    return max(feasible,default=0)


def sweep(c):
    rows=[]
    for epsilon in (.001,.01):
        for scenario,changes in c['scenarios'].items():
            for tau in c['times']:
                for arm in ('nisq','star','ft'):
                    for L in range(2,c['max_L']+1,2):
                        rows.append({**estimate(L,tau,arm,{**c,**changes},epsilon),'scenario':scenario})
    return rows


def plot(rows,c,epsilon=.001):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from plot_capacity import COLORS
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,len(c['times']),figsize=(14.8,5.6),sharey=True)
    labels={'nisq':'NISQ · effective mitigation','star':'STAR · analog rotations','ft':'FTQC · cultivated T states'}
    ymax=4*ceil(max(capacity([v for v in rows if v['tau']==tau and v['arm']==arm and v['scenario']==s and v['absolute_error']==epsilon],1e9)
                   for tau in c['times'] for arm in labels for s in c['scenarios'])/4)
    for ax,tau in zip(axes,c['times']):
        for arm in labels:
            subset=[v for v in rows if v['arm']==arm and v['tau']==tau and v['absolute_error']==epsilon]
            thresholds=[10**v['log10_Q_min'] for v in subset if v['log10_Q_min'] is not None and 1<=v['log10_Q_min']<=9]
            Q=np.array(sorted(set(np.geomspace(10,1e9,220).tolist()+thresholds)))
            values={s:np.array([capacity([v for v in subset if v['scenario']==s],q) for q in Q],dtype=float)
                    for s in c['scenarios']}
            lower=np.minimum.reduce(list(values.values()));upper=np.maximum.reduce(list(values.values()))
            ax.fill_between(Q,lower,upper,step='post',color=COLORS[arm],alpha=.12)
            center=values['central'];center[center==0]=np.nan
            ax.step(Q,center,where='post',color=COLORS[arm],lw=2.2,label=labels[arm])
        ax.set_xscale('log');ax.set_xlim(10,1e9);ax.set_ylim(0,ymax+1)
        ax.set_xticks([1e2,1e4,1e6,1e8]);ax.set_yticks(range(0,ymax+1,4 if ymax<=36 else 8))
        ax.set_xlabel('Available physical qubits, Q')
        ax.set_title(f'τ = {tau:g}  ·  r = {ceil(c["steps_per_tau"]*tau)}',loc='left',fontsize=13)
        ax.grid(alpha=.18)
        if ymax>=c['max_L']:
            ax.axhspan(c['max_L']-.4,c['max_L']+1,color='#eeeeee',zorder=0)
            ax.text(.98,.975,f'L = {c["max_L"]}: sweep limit',ha='right',va='top',transform=ax.transAxes,fontsize=8,color='#666666')
    axes[0].set_ylabel('Largest modeled lattice side L\n(L × L Hubbard sites)')
    handles,labels_=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels_,loc='lower center',bbox_to_anchor=(.5,.095),ncol=3,frameon=False)
    fig.suptitle('Back-of-the-envelope+: what can a fixed qubit budget buy?',fontsize=17,x=.065,ha='left',y=.985)
    precision='ε = 0.001 (1% relative at D* = 0.1)' if epsilon==.001 else 'ε = 0.01 absolute sensitivity'
    fig.text(.065,.91,f'One week · p = 10⁻³ · {precision} · 10 ns gates / 240 ns code cycles',fontsize=10,color='#555555')
    fig.text(.5,.025,'Bands: declared scenarios, not confidence intervals. Effective NISQ overhead is assumed; Trotter and mitigation biases are assumed affordable.\n'
             'STAR uses termwise Strang; FTQC uses plaquettes; NISQ uses a Willow-anchored proxy. Even L only; full-lattice costs, without physical-cone truncation.',
             ha='center',fontsize=8.2,color='#555555')
    fig.subplots_adjust(left=.065,right=.985,bottom=.25,top=.80,wspace=.12)
    name='discussion_capacity' if epsilon==.001 else 'discussion_capacity_loose'
    for ext in ('pdf','png','svg'):
        path=ROOT/f'figures/{name}.{ext}';fig.savefig(path,dpi=200)
        if ext=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


def run():
    c=assumptions();rows=sweep(c);out=ROOT/'results/discussion_v1';out.mkdir(exist_ok=True)
    (out/'assumptions.json').write_text(json.dumps(c,indent=2)+'\n')
    (out/'costs.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    summaries=[]
    for eps in (.001,.01):
        for tau in c['times']:
            for Q in c['qubit_budgets']:
                values={arm:capacity([v for v in rows if v['arm']==arm and v['scenario']=='central' and v['tau']==tau and v['absolute_error']==eps],Q)
                        for arm in ('nisq','star','ft')}
                summaries.append(dict(absolute_error=eps,tau=tau,Q=Q,**values))
    with (out/'capacity.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
    for eps in (.001,.01):plot(rows,c,eps)
    summary=['# Back-of-the-envelope+ discussion example','',
             'Central scenarios, one week, p = 0.001, 10 ns gates, 240 ns syndrome cycles.',
             'These are analytical planning estimates; mitigation and Trotter biases are assumed affordable.',
             'NISQ gate/depth counts are anchored to Willow. Its assumed sampling penalty is 100 at the anchor, not a measured fit.','',
             '## Largest lattice at tau = 1, absolute error 0.001','',
             '| Physical qubits | NISQ | STAR | FTQC |','|---:|---:|---:|---:|']
    def size(L):return f'{L} × {L}' if L else '—'
    for v in summaries:
        if v['absolute_error']==.001 and v['tau']==1 and v['Q'] in (100,10000,100000,1000000,10000000):
            summary.append(f'| {v["Q"]:,} | {size(v["nisq"])} | {size(v["star"])} | {size(v["ft"])} |')
    summary+=['','## Area required for 6 × 6 at tau = 1','',
              '| Absolute precision | NISQ | STAR | FTQC |','|---:|---:|---:|---:|']
    for eps in (.001,.01):
        selected={arm:next(v for v in rows if v['L']==6 and v['tau']==1 and v['arm']==arm and v['absolute_error']==eps and v['scenario']=='central') for arm in ('nisq','star','ft')}
        summary.append('| '+str(eps)+' | '+' | '.join(f'{round(10**selected[arm]["log10_Q_min"]):,}' for arm in ('nisq','star','ft'))+' |')
    summary+=['','The strict FT case uses an extra magic-state cleanup stage; the loose case does not.',
              'At fixed circuit size, sampling scales as inverse precision squared. Integer replicas and code distance add thresholds.','',
              'NISQ: exp(kappa p G2) effective sampling cost; kappa is fixed by the chosen anchor overhead.',
              'STAR: exp(4 alpha p Theta + 4 f R), with four effective concurrent rotation lanes.',
              'FTQC: synthesized T count, with preparation and consumption charged together.','',
              'Bands vary NISQ anchor overhead 10–1000, STAR floor 0–1e-4, delivery lanes 8–1, and factory service time 0.5–4 times central.',
              'They are scenarios, not confidence intervals.','',
              'Details: [discussion guide](../../docs/DISCUSSION_MODEL.md).',
              'Figures: [strict precision](../../figures/discussion_capacity.pdf), [0.01 sensitivity](../../figures/discussion_capacity_loose.pdf).','']
    (out/'summary.md').write_text('\n'.join(summary))
    for v in summaries:
        if v['absolute_error']==.001 and v['Q'] in (100,10000,100000,1000000):print(v)
    print('Wrote',len(rows),'analytical scenario points; no compiler or dynamics jobs run.')


if __name__=='__main__':run()
