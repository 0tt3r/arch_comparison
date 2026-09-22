"""Regenerate optimized analytical curves from saved time calibration."""
import csv
import json
from math import ceil,log10

from budget_optimization import algorithm_for,calibration,coefficient,optimize,settings
from discussion_model import assumptions,capacity,estimate,families,plot
from project_paths import ROOT

OUT=ROOT/'results/discussion_v3'


def write_csv(name,rows):
    with (OUT/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def calibration_plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    data=calibration()['rows']
    rectangle=OUT/'rectangle_checks/calibration.json'
    if rectangle.exists():data=data+json.loads(rectangle.read_text())['rows']
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),sharey=True)
    for ax,algorithm in zip(axes,('plaquette','termwise')):
        for patch,style,color in [('2x2',':','#888888'),('3x3','-','#14667e'),('4x3','none','#c46620')]:
            rows=sorted([v for v in data if v['patch']==patch and v['algorithm']==algorithm],key=lambda v:v['tau'])
            if rows:ax.plot([v['tau'] for v in rows],[v['coefficient'] for v in rows],ls=style,marker='o',ms=4,color=color,label=patch.replace('x',' × '))
        rows=sorted([v for v in data if v['patch']=='3x3' and v['algorithm']==algorithm],key=lambda v:v['tau'])
        ax.plot([v['tau'] for v in rows],[coefficient(v['tau'],algorithm) for v in rows],'--',color='#14667e',lw=1,label='planning envelope')
        times=np.linspace(.1,8,120)
        ax.plot(times,.22*times**3,'--',color='#a25d84',lw=1,label='previous 0.22 τ³')
        ax.set_yscale('log');ax.set_xlabel('Evolution time τ');ax.set_title(algorithm.capitalize(),loc='left')
        ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=9)
    axes[0].set_ylabel('Measured coefficient |Dᵣ − Dexact| r²')
    fig.suptitle('Trotter error has ordering and finite-patch dependence',x=.08,ha='left')
    fig.text(.5,.015,'Tail coefficients use two finest step counts after a convergence check. Finite patches do not bound bulk error.',ha='center',fontsize=9,color='#555555')
    fig.tight_layout(rect=(0,.05,1,.94))
    for ext in ('png','pdf','svg'):
        p=ROOT/f'figures/time_calibration.{ext}';fig.savefig(p,dpi=180)
        if ext=='svg':p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
    plt.close(fig)


def run():
    c=assumptions();options=settings();rows=[]
    max_tau=max(v['tau'] for v in calibration()['rows'])
    for epsilon in (.001,.01):
        for family,tau_fixed in families(c):
            for arm in ('nisq','star','ft','hybrid'):
                scenarios=dict(c['scenarios'])
                if arm=='nisq':scenarios.update(pec_two_qubit={'nisq_model':'pec_two_qubit'},compiled_pec={'nisq_model':'compiled_pec'})
                for scenario,changes in scenarios.items():
                    for L in range(2,c['max_L']+1,2):
                        tau=tau_fixed if tau_fixed is not None else c['time_family_scale']*L
                        if tau>max_tau:continue
                        row=optimize(L,tau,arm,{**c,**changes},epsilon,options)
                        rows.append({**row,'family':family,'scenario':scenario})
                print('completed',epsilon,family,arm,flush=True)
    OUT.mkdir(exist_ok=True)
    (OUT/'assumptions.json').write_text(json.dumps({'hardware_and_delivery':c,'optimization':options},indent=2)+'\n')
    (OUT/'costs.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    summaries=[]
    for epsilon in (.001,.01):
        for family,_ in families(c):
            for Q in c['qubit_budgets']:
                v=dict(absolute_error=epsilon,family=family,Q=Q)
                for arm in ('nisq','star','hybrid','ft'):
                    v[arm]=capacity([x for x in rows if x['arm']==arm and x['scenario']=='central' and x['family']==family and x['absolute_error']==epsilon],Q)
                for scenario in ('pec_two_qubit','compiled_pec'):
                    v[scenario]=capacity([x for x in rows if x['arm']=='nisq' and x['scenario']==scenario and x['family']==family and x['absolute_error']==epsilon],Q)
                summaries.append(v)
    write_csv('capacity.csv',summaries)
    budget_rows=[]
    for v in rows:
        if v.get('budget_fractions'):
            budget_rows.append({**{k:v.get(k) for k in ('L','tau','arm','absolute_error','family','scenario','r','distance','cleanup_levels','factories','log10_Q_min')},**v['budget_fractions']})
    write_csv('budgets.csv',budget_rows)
    ablation=[]
    for epsilon in (.001,.01):
        for L in (4,6):
            for tau in (1,4,8):
                for arm in ('nisq','star','hybrid','ft'):
                    old=estimate(L,tau,arm,c,epsilon)
                    opt_old=optimize(L,tau,arm,c,epsilon,options,C_override=.22*tau**3)
                    pointwise=optimize(L,tau,arm,c,epsilon,{**options,'coefficient_policy':'pointwise'})
                    current=next(v for v in rows if v['absolute_error']==epsilon and v['L']==L and v['tau']==tau and v['family']==f'tau_{tau}' and v['arm']==arm and v['scenario']=='central')
                    ablation.append(dict(L=L,tau=tau,arm=arm,absolute_error=epsilon,
                        fixed_v2_log10_Q=old['log10_Q_min'],optimized_old_coefficient_log10_Q=opt_old['log10_Q_min'],
                        optimized_measured_coefficient_log10_Q=current['log10_Q_min'],
                        optimized_pointwise_coefficient_log10_Q=pointwise['log10_Q_min'],
                        old_r=old['r'],optimized_old_r=opt_old.get('r'),new_r=current.get('r')))
    write_csv('ablation.csv',ablation)
    for epsilon in (.001,.01):plot(rows,{**c,'presentation_version':'v3'},epsilon)
    calibration_plot()
    def area(v):
        q=v['log10_Q_min']
        return 'outside model' if q is None else (f'{round(10**q):,}' if q<9 else f'~10^{q:.1f}')
    def size(L):return f'{L} × {L}' if L else '—'
    text=['# Optimized discussion example','',
        'One week, p = 0.001; measured finite-patch time dependence and per-arm error allocation. The total precision is unchanged.',
        'The planning coefficient is 1.1 times the running maximum of converged 3 × 3 coefficients for each ordering through the requested time. This avoids relying on isolated finite-patch cancellations. Pointwise coefficients are a sensitivity; transfer to larger lattices remains an assumption. No global time power is fitted. L = 2 uses its own measured coefficients.',
        'The coefficient also covers measured finite-step corrections at the minimum allowed r and twice r; selected optimized counts are validated directly.',
        'Statistics receive the budget remaining after Trotter, logical, synthesis, magic, and residual-bias costs. Mitigated arms reserve 10% for uncalibrated residual bias. Its affordability is assumed.',
        'Optimization searches the declared allocation grid, integer-step neighbors, and the shared cultivated-state menu with zero or one cleanup stage. It is not a proof of the global optimum.','',
        '## Central areas at tau = 1','',
        '| Lattice | Precision | Arm | Steps | Q | Cleanup | Statistical share | Trotter share | Magic share |',
        '|---|---:|---|---:|---:|---:|---:|---:|---:|']
    for L in (4,6):
        for epsilon in (.001,.01):
            for arm in ('nisq','star','hybrid','ft'):
                v=next(x for x in rows if x['L']==L and x['family']=='tau_1' and x['absolute_error']==epsilon and x['arm']==arm and x['scenario']=='central')
                b=v['budget_fractions']
                text.append(f'| {size(L)} | {epsilon} | {arm} | {v["r"]} | {area(v)} | {v.get("cleanup_levels",0)} | {b["stat"]:.3f} | {b["trotter"]:.3f} | {b["magic"]:.3f} |')
    text+=['','## Largest lattice at tau = 1, strict precision','',
           '| Q | NISQ effective | NISQ 2q PEC | NISQ compiled PEC | STAR | Hybrid | FTQC |',
           '|---:|---:|---:|---:|---:|---:|---:|']
    for v in summaries:
        if v['absolute_error']==.001 and v['family']=='tau_1':
            text.append('| '+f'{v["Q"]:,}'+' | '+' | '.join(size(v[k]) for k in ('nisq','pec_two_qubit','compiled_pec','star','hybrid','ft'))+' |')
    text+=['','The [ablation table](ablation.csv) separates budget optimization under the old coefficient from replacing the coefficient. The [budget table](budgets.csv) records every selected split.',
           'The time-growing panel stops at L = 32 because tau = L/4 then reaches the measured endpoint tau = 8. The fixed-time sweeps still reach L = 48.',
           'STAR uses the termwise ordering; FT/hybrid use plaquettes. The NISQ gate proxy transfers the plaquette coefficient, without claiming a match to the experimental ordering.',
           'FT and hybrid use the same source quality, service, and cleanup menu. Source-quality/throughput tradeoffs beyond this menu remain uncalibrated.',
           'PEC references retain Hoeffding counts. Other curves use Gaussian planning counts. Bands show delivery/mitigation scenarios, not confidence intervals or bulk-calibration uncertainty.','',
           'Figures: [strict](../../figures/discussion_capacity.pdf), [loose](../../figures/discussion_capacity_loose.pdf), [time calibration](../../figures/time_calibration.pdf).',
           'Full story: [main manuscript](../../tex/main.pdf). Previous [fixed-allocation example](../discussion_v2/summary.md).','']
    (OUT/'summary.md').write_text('\n'.join(text))
    print('Wrote',len(rows),'optimized scenario points.',flush=True)


if __name__=='__main__':run()
