"""Area, consumption, and source supply optimized on the same error contract."""
import csv
import argparse
import hashlib
import json
from math import ceil, exp, isclose, log, log10

from budget_optimization import optimize
from discussion_model import assumptions, capacity
from project_paths import ROOT

OUT=ROOT/'results/delivery_v1'


def settings():
    return json.loads((ROOT/'configs/delivery_sensitivity.json').read_text())


def choices(arm,config):
    yield 'reference',{'data_block':'reference'}
    if arm in ('star','hybrid'):
        for b in config['analog_lane_choices']:
            yield f'lanes_{b}',{'data_block':'reference','parallel_lanes':b}
    if arm in ('ft','hybrid'):
        for block,params in config['blocks'].items():
            # FT whole-block comparison is serial. Hybrid analog lanes are
            # explored separately; its total T rate remains one block's rate.
            for b in config['analog_lane_choices'] if arm=='hybrid' else [1]:
                yield f'{block}_b{b}',{'data_block':block,'ft_clocks_per_T':params['clocks_per_T'],
                                     'parallel_lanes':b,'T_consumption_ports':1}
        yield 'fast_four_ports',{'data_block':'fast','ft_clocks_per_T':1,
                                'parallel_lanes':4,'T_consumption_ports':4}


def key(row):
    return row['log10_Q_min'],log(row['qubits_per_copy'])+log(row['seconds_per_attempt'])+log(10)*row['log10_total_shots']


def run():
    config=settings();c={**assumptions(),'factory_choices':config['factory_choices']}
    rows=[]
    for eps in config['precisions']:
        for L in config['sizes']:
            for arm in ('star','hybrid','ft'):
                for name,changes in choices(arm,config):
                    v=optimize(L,config['tau'],arm,{**c,**changes},eps)
                    rows.append({**v,'layout':name,'delivery_assumptions':changes})
            print('delivery',eps,L,flush=True)
    OUT.mkdir(exist_ok=True)
    (OUT/'costs.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    (OUT/'assumptions.json').write_text(json.dumps({'delivery':config,'base':c},indent=2)+'\n')
    chosen=[]
    for eps in config['precisions']:
        for L in config['sizes']:
            for arm in ('star','hybrid','ft'):
                group=[v for v in rows if v['absolute_error']==eps and v['L']==L and v['arm']==arm and v['log10_Q_min'] is not None]
                if not group:continue
                for menu in ('reference','whole_block_and_lanes','including_four_ports'):
                    eligible=[v for v in group if (v['layout']=='reference' if menu=='reference' else
                         v['layout']!='fast_four_ports' if menu=='whole_block_and_lanes' else True)]
                    # Reference and reference-lane entries remain in the menu:
                    # it is a sensitivity envelope, not a routing certification.
                    winner=min(eligible,key=key)
                    chosen.append({**winner,'menu':menu})
    (OUT/'selected.json').write_text(json.dumps(chosen,indent=2,allow_nan=False)+'\n')
    columns=['absolute_error','L','arm','menu','layout','r','distance','lanes','factories','cleanup_levels',
             'data_tiles','qubits_per_copy','seconds_per_attempt','log10_Q_min','supply_seconds','consumption_seconds','delivery_bottleneck']
    with (OUT/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows({k:v.get(k) for k in columns} for v in chosen)
    plot(chosen,rows)
    validate(rows,chosen)


def plot(chosen,rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from plot_capacity import COLORS
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.9),sharey=True)
    menus=[('reference',':','Existing four-lane reference'),('whole_block_and_lanes','-','Block / analog-lane menu'),
           ('including_four_ports','--','Also allow four fast T ports')]
    for ax,eps in zip(axes,(.001,.01)):
        for arm,label in [('star','STAR'),('hybrid','Hybrid'),('ft','FTQC')]:
            for menu,ls,_ in menus:
                group=[v for v in chosen if v['absolute_error']==eps and v['arm']==arm and v['menu']==menu]
                Q=sorted(set(np.geomspace(1e4,1e9,200))|{10**v['log10_Q_min'] for v in group if 4<=v['log10_Q_min']<=9})
                L=np.array([capacity(group,q) for q in Q],float);L[L==0]=np.nan
                ax.step(Q,L,where='post',color=COLORS[arm],ls=ls,lw=1.8,label=label if menu=='whole_block_and_lanes' else None)
        ax.set_xscale('log');ax.set_xlim(1e4,1e9);ax.set_ylim(0,50);ax.set_yticks(range(0,49,8));ax.grid(alpha=.2)
        ax.set_title(f'ε = {eps:g}',loc='left');ax.set_xlabel('Available physical qubits, Q')
    axes[0].set_ylabel('Largest modeled side L (L × L sites)');axes[0].legend(frameon=False,loc='upper left')
    from matplotlib.lines import Line2D
    fig.legend([Line2D([],[],color='#555555',ls=ls) for _,ls,_ in menus],[label for _,_,label in menus],
               loc='lower center',bbox_to_anchor=(.5,.07),ncol=3,frameon=False,fontsize=9)
    fig.suptitle('Encoded delivery: faster gates need area and state supply',x=.08,ha='left')
    fig.text(.5,.02,'τ = 1 · one week · p = 10⁻³ · budgets and factories reoptimized · sparse L grid (up to 48)\n'
             'Analog-lane and multiple-port rates are assumptions; the small-patch accuracy envelope is unchanged.',ha='center',fontsize=9,color='#555555')
    fig.tight_layout(rect=(0,.17,1,.94))
    for ext in ('pdf','png','svg'):
        p=ROOT/f'figures/delivery_capacity.{ext}';fig.savefig(p,dpi=180)
        if ext=='svg':p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
    plt.close(fig)
    # Direct supply/consumption evidence at the main point, independent of replicas.
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.4))
    for ax,arm in zip(axes,('ft','hybrid')):
        names=['reference','compact_b1','intermediate_b1','fast_b1','fast_four_ports']
        group=[next(v for v in rows if v['absolute_error']==.001 and v['L']==4 and v['arm']==arm and v['layout']==name) for name in names]
        x=np.arange(len(group));width=.36
        ax.bar(x-width/2,[v['supply_seconds'] for v in group],width,label='State supply',color='#ae783d')
        ax.bar(x+width/2,[v['consumption_seconds'] for v in group],width,label='Data-block consumption',color='#426b88')
        ax.set_xticks(x,['Reference\n4 lanes','Compact\n1 port','Intermediate\n1 port','Fast\n1 port','Fast\n4 ports'])
        ax.set_title('FTQC' if arm=='ft' else 'Hybrid',loc='left');ax.set_yscale('log');ax.set_ylabel('Seconds per shot (overlapped stages)')
        ax.grid(axis='y',alpha=.2);ax.legend(frameon=False,fontsize=9)
    fig.suptitle('At 4 × 4, supply limits the fastest FT sensitivity',x=.08,ha='left')
    fig.text(.5,.015,'ε = 0.001 · τ = 1 · max(supply, consumption) sets duration, plus initialization and feedforward.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,.94))
    for ext in ('pdf','png'):
        fig.savefig(ROOT/f'figures/delivery_bottlenecks.{ext}',dpi=180)
    plt.close(fig)


def validate(rows,chosen):
    c=assumptions()
    for v in rows:
        if v['log10_Q_min'] is None:continue
        assert abs(sum(v['budget_fractions'].values())-1)<1e-12
        assert v['seconds_per_attempt']>=max(v.get('supply_seconds',0),v.get('consumption_seconds',0))
        d=v['distance'];sources=v.get('factories',0);eps=v['absolute_error'];budget=v['budget_fractions']
        tiles=(v['qubits_per_copy']-sources*c['source_qubits'])/(2*d*d)
        log_fault=log(tiles*ceil(v['seconds_per_attempt']/c['cycle_seconds'])*c['surface_A'])+(d+1)/2*log(c['p']/c['surface_threshold'])
        assert log_fault+v['chi']/2<=log(eps*budget['logical'])+1e-12
        if v['arm'] in ('ft','hybrid'):
            assert log(v['T_count']*v['magic_error'])+v['chi']/2<=log(eps*budget['magic'])+1e-12
        shots=ceil(c['gaussian_z']**2*c['signal_reference']*(1-c['signal_reference'])/(eps*budget['stat'])**2)
        assert isclose(log10(shots)+v['chi']/log(10),v['log10_total_shots'],abs_tol=1e-10)
        lr=max(0,log(shots)+v['chi']+log(v['seconds_per_attempt']/c['week_seconds']))
        if lr<30:lr=log(ceil(exp(lr)))
        assert isclose((log(v['qubits_per_copy'])+lr)/log(10),v['log10_Q_min'],abs_tol=1e-10)
        if v['layout'].split('_')[0] in ('compact','intermediate','fast') and v['layout']!='fast_four_ports':
            assert v['T_consumption_ports']==1
    for v in chosen:
        reference=next(x for x in chosen if x['L']==v['L'] and x['arm']==v['arm'] and x['absolute_error']==v['absolute_error'] and x['menu']=='reference')
        assert v['log10_Q_min']<=reference['log10_Q_min']+1e-12
        assert v.get('factories',0)<max(settings()['factory_choices']),'Expand the source menu before claiming an interior optimum'
    inputs=['python/discussion_model.py','python/budget_optimization.py','python/delivery_sensitivity.py',
            'configs/delivery_sensitivity.json','configs/discussion_model.json','configs/assumptions.json',
            'configs/budget_optimization.json','results/discussion_v3/calibration.json','results/discussion_v3/finite_step_guard.json']
    report={'status':'PASS: budgets, logical and magic faults, replicas, no duplicated T rate, delivery lower bounds, menu dominance, and source-menu endpoint',
            'rows':len(rows),'selected':len(chosen),
            'sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs}}
    (OUT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.check:
        validate(json.loads((OUT/'costs.json').read_text()),json.loads((OUT/'selected.json').read_text()))
    else:run()
