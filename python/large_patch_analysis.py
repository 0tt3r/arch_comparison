"""Directly measured 4x4 step menu beside transferred planning coefficients."""
import csv
import argparse
import hashlib
import json
from math import isclose

from budget_optimization import algorithm_for, coefficient, optimize, settings as budget_settings
from delivery_sensitivity import choices, key, settings as delivery_settings
from discussion_model import assumptions, estimate
from frozen_baseline import verify_sources
from project_paths import ROOT
from small_dynamics import rectangle

OUT=ROOT/'results/dynamics_4x4'


def measurements():
    files=sorted(OUT.glob('*/*x*_tau*_*.json'))
    records=[(p,json.loads(p.read_text())) for p in files]
    records=[(p,v) for p,v in records if v['contract']['algorithm']!='benchmark']
    result=[]
    build=json.loads((OUT/'build.json').read_text())
    for p,v in records:
        c=v['contract'];data=v['result']
        assert c['source_sha256']==hashlib.sha256((ROOT/'cpp/hubbard_sector.cpp').read_bytes()).hexdigest()
        assert c['driver_sha256']==hashlib.sha256((ROOT/'python/large_dynamics.py').read_bytes()).hexdigest()
        assert c['binary_sha256']==build['sha256']['build/hubbard_sector']
        assert v['coefficient_sha256']==hashlib.sha256(p.with_suffix('.coefficients.txt').read_bytes()).hexdigest()
        assert v['status']=='complete'
        assert abs(data['norm_squared']-1)<1e-8
        nx,ny=c['nx'],c['ny'];algorithm=c['algorithm']
        if algorithm=='exact':continue
        refs=[x for _,x in records if x['contract']['nx']==nx and x['contract']['ny']==ny
              and x['contract']['tau']==c['tau'] and x['contract']['algorithm']=='exact']
        assert len(refs)==1
        reference=refs[0]
        assert reference['chebyshev_state_tail_bound']<1e-12
        sites=rectangle(nx,ny)
        # Save both the compiled square's observed coordinate and the common
        # patch coordinate. Ordered Trotter circuits need not share reflections.
        for xy in sorted(set([(1,1),(nx//2,ny//2)])):
            i=sites.index(xy);D=reference['result']['doublons'][i]
            signed=data['doublons'][i]-D
            result.append(dict(patch=f'{nx}x{ny}',tau=c['tau'],algorithm=algorithm,r=c['r'],
                site=list(xy),exact_doublon=D,doublon=data['doublons'][i],signed_error=signed,
                coefficient=abs(signed)*c['r']**2,norm_error=abs(data['norm_squared']-1),
                source=str(p.relative_to(ROOT))))
    return result


def run():
    rows=measurements()
    large=[v for v in rows if v['patch']=='4x4' and v['site']==[2,2] and v['tau']==1]
    for algorithm in ('plaquette','termwise'):
        group=sorted([v for v in large if v['algorithm']==algorithm],key=lambda v:v['r'])
        assert {8,16,32}<=set(v['r'] for v in group),f'Incomplete {algorithm} step menu'
        a,b=[v['coefficient'] for v in group[-2:]]
        assert abs(a-b)/max(a,b)<.05,f'Unconverged {algorithm} coefficient'
    config=delivery_settings();base={**assumptions(),'factory_choices':config['factory_choices']}
    options=budget_settings();costs=[];selected=[];interpolated=[]
    for eps in (.001,.01):
        for arm in ('nisq','star','hybrid','ft'):
            algorithm=algorithm_for(arm)
            measured=[v for v in large if v['algorithm']==algorithm]
            menu=list(choices(arm,config)) if arm!='nisq' else [('effective',{}),('pec_two_qubit',{'nisq_model':'pec_two_qubit'}),('compiled_pec',{'nisq_model':'compiled_pec'})]
            C=1.1*max([v['coefficient'] for v in measured]+[coefficient(1,algorithm)/1.1])
            for name,changes in menu:
                fitted=optimize(4,1,arm,{**base,**changes},eps,C_override=C)
                interpolated.append({**fitted,'layout':name,'coefficient_patch':'4x4',
                    'coefficient_policy':'1.1 times largest measured coefficient, retaining the old planning envelope as a floor',
                    'scope':'Step interpolation from directly measured 4x4 values; selected r is not necessarily directly checked'})
                group=[]
                for v in measured:
                    # Direct observable error, with a declared 1e-10 numerical
                    # allowance. No small-patch transfer or r interpolation.
                    error=abs(v['signed_error'])+1e-10
                    row=optimize(4,1,arm,{**base,**changes},eps,
                        {**options,'explicit_step_candidates':[v['r']]},C_override=error*v['r']**2)
                    row.update(layout=name,r=v['r'],measured_trotter_error=abs(v['signed_error']),
                               numerical_allowance=1e-10,coefficient_patch='4x4',
                               coefficient_policy='direct measured error at this step count',
                               measurement_source=v['source'])
                    costs.append(row)
                    if row['log10_Q_min'] is not None:group.append(row)
                if group:selected.append(min(group,key=key))
    (OUT/'measurements.json').write_text(json.dumps(rows,indent=2)+'\n')
    (OUT/'costs.json').write_text(json.dumps(costs,indent=2,allow_nan=False)+'\n')
    (OUT/'selected.json').write_text(json.dumps(selected,indent=2,allow_nan=False)+'\n')
    (OUT/'interpolated.json').write_text(json.dumps(interpolated,indent=2,allow_nan=False)+'\n')
    columns=['absolute_error','arm','layout','r','distance','factories','cleanup_levels','log10_Q_min','measured_trotter_error']
    with (OUT/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows({k:v.get(k) for k in columns} for v in selected)
    plot(rows)
    text=['# Direct 4 × 4 calibration at tau = 1','',
          'Central coordinate (2,2); OBC Neel quench at U/J=4. All rows retain the reference signal D*=0.1 for sampling.',
          'This is a measured-step comparison; it does not replace the larger-lattice planning envelope.',
          'NISQ still transfers plaquette dynamics to its experimental circuit-count proxy, and residual mitigation bias remains assumed.','',
          '| Ordering | r | Exact doublon | Signed Trotter error | abs(error) r² |',
          '|---|---:|---:|---:|---:|']
    for v in sorted(large,key=lambda v:(v['algorithm'],v['r'])):
        text.append(f"| {v['algorithm']} | {v['r']} | {v['exact_doublon']:.10f} | {v['signed_error']:.8g} | {v['coefficient']:.6f} |")
    text+=['','Resource allocations below optimize only over directly measured r values and the existing budget grid. A 1e-10 numerical allowance is added to the measured error.','',
           '| Error | Arm | Layout | r | Physical qubits |','|---:|---|---|---:|---:|']
    for v in selected:
        if v['layout'] in ('reference','fast_b1','fast_four_ports','lanes_2','effective','pec_two_qubit','compiled_pec'):
            text.append(f"| {v['absolute_error']} | {v['arm']} | {v['layout']} | {v['r']} | {10**v['log10_Q_min']:.5g} |")
    text+=['','The doubling grid can penalize an arm that strongly prefers intermediate r. The following **interpolated sensitivity** retains the original continuous-in-r error model, now with a 4x4 coefficient and 10% margin. Its selected r values are not all directly checked.','',
           '| Error | Arm | Layout | r | 4x4 planning coefficient | Physical qubits |','|---:|---|---|---:|---:|---:|']
    for v in interpolated:
        if v['layout'] in ('reference','fast_b1','fast_four_ports','effective','pec_two_qubit','compiled_pec'):
            text.append(f"| {v['absolute_error']} | {v['arm']} | {v['layout']} | {v['r']} | {v['coefficient']:.6g} | {10**v['log10_Q_min']:.5g} |")
    (OUT/'summary.md').write_text('\n'.join(text)+'\n')
    inputs=[ROOT/p for p in ['cpp/hubbard_sector.cpp','python/large_dynamics.py','python/large_patch_analysis.py',
            'python/budget_optimization.py','python/discussion_model.py','python/delivery_sensitivity.py',
            'python/resource_extensions.py',
            'configs/assumptions.json','configs/discussion_model.json','configs/budget_optimization.json',
            'configs/delivery_sensitivity.json']]
    inputs+=sorted(OUT.glob('*/*x*_tau*_*.json'))
    for v in selected:
        assert abs(sum(v['budget_fractions'].values())-1)<1e-12
        assert abs(v['absolute_error']*v['budget_fractions']['trotter']-v['measured_trotter_error']-1e-10)<1e-12
    report={'status':'PASS: direct step census, tail convergence, source hashes, norm checks, and measured-error budgets',
            'measured_rows':len(rows),'candidate_rows':len(costs),'selected_rows':len(selected),
            'scope':'4x4 tau=1 only; no bulk or longer-time accuracy claim; NISQ ordering and mitigation remain proxies',
            'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
    (OUT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print('\n'.join(text))


def plot(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    import numpy as np
    small=json.loads((ROOT/'results/discussion_v3/calibration.json').read_text())['rows']
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.5))
    for ax,algorithm in zip(axes,('plaquette','termwise')):
        for patch,color in [('3x3','#888888'),('4x3','#bb772d'),('4x4','#14667e')]:
            if patch=='3x3':
                v=next(v for v in small if v['patch']==patch and v['tau']==1 and v['algorithm']==algorithm)
                r=[m['r'] for m in v['measurements']];error=[abs(m['signed_error']) for m in v['measurements']]
            else:
                site=[2,2] if patch=='4x4' else [1,1]
                group=sorted([v for v in rows if v['patch']==patch and v['algorithm']==algorithm and v['site']==site],key=lambda v:v['r'])
                r=[v['r'] for v in group];error=[abs(v['signed_error']) for v in group]
            ax.loglog(r,error,'o-',color=color,label=patch.replace('x',' × '))
        r=np.geomspace(8,64,100)
        ax.loglog(r,coefficient(1,algorithm)/r**2,'--',color='#333333',label='3 × 3 planning envelope')
        ax.set_title(algorithm.capitalize(),loc='left');ax.set_xlabel('Trotter steps r');ax.set_xticks([8,16,32,64],['8','16','32','64'])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=9)
    axes[0].set_ylabel('Absolute central doublon error')
    fig.suptitle('Direct finite-lattice calibration at τ = 1',x=.08,ha='left')
    fig.text(.5,.01,'4 × 4 observes (2,2); smaller patches observe (1,1). Exact references use the same finite OBC lattice.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,.94))
    for ext in ('pdf','png'):fig.savefig(ROOT/f'figures/large_patch_calibration.{ext}',dpi=180)
    plt.close(fig)


def check():
    assert measurements()==json.loads((OUT/'measurements.json').read_text())
    report=json.loads((OUT/'validation.json').read_text())
    verify_sources(report['sha256'])
    config=delivery_settings();base={**assumptions(),'factory_choices':config['factory_choices']}
    for row in json.loads((OUT/'selected.json').read_text())+json.loads((OUT/'interpolated.json').read_text()):
        menu=dict(choices(row['arm'],config)) if row['arm']!='nisq' else {
            'effective':{},'pec_two_qubit':{'nisq_model':'pec_two_qubit'},'compiled_pec':{'nisq_model':'compiled_pec'}}
        b=row['budget_fractions']
        c={**base,**menu[row['layout']],**{k+'_fraction':b[k] for k in ('stat','logical','magic','synthesis')},
           'forced_cleanup':row.get('cleanup_levels',0)}
        replay=estimate(row['L'],row['tau'],row['arm'],c,row['absolute_error'],steps=row['r'])
        assert isclose(replay['log10_Q_min'],row['log10_Q_min'],abs_tol=1e-12)
    for row in json.loads((OUT/'selected.json').read_text()):
        assert abs(sum(row['budget_fractions'].values())-1)<1e-12
        assert abs(row['absolute_error']*row['budget_fractions']['trotter']-row['measured_trotter_error']-1e-10)<1e-12
    jobs=json.loads((OUT/'jobs.json').read_text())
    resumed=json.loads((OUT/'resume_note.json').read_text())
    for job in jobs:
        assert job['job_state']=='F'
        if job['jobid']==resumed['job']:
            assert job['Exit_status']==271 and resumed['completed_cases_preserved']
        else:assert job['Exit_status']==0
    for v in json.loads((OUT/'convergence.json').read_text())['rows']:
        assert v['relative_tail_coefficient_spread']<=.05
    # Memory placement changes execution, not the finite-lattice calculation.
    numa=[json.loads(p.read_text()) for p in sorted((OUT/'numa_smoke').glob('*/*.json'))]
    assert len(numa)==3 and len({v['host'] for v in numa})==1
    reference=numa[0]['result']['doublons']
    assert all(max(abs(a-b) for a,b in zip(reference,v['result']['doublons']))<1e-10 for v in numa)
    print(report['status'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.check:check()
    else:run()
