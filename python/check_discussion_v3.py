"""Independent budget arithmetic, selected-step physics, and grid sensitivity."""
import argparse
import csv
import hashlib
import json
from math import ceil,exp,isclose,log,log10

from budget_optimization import optimize,settings
from discussion_model import assumptions,capacity
from project_paths import ROOT
from time_calibration import SplitSector

OUT=ROOT/'results/discussion_v3'


def step_guard():
    source=json.loads((OUT/'calibration.json').read_text())
    opts=settings();models={};rows=[];exact={}
    for v in source['rows']:
        patch=v['patch'];tau=v['tau'];algorithm=v['algorithm'];key=(patch,tau)
        if patch not in models:models[patch]=SplitSector(*map(int,patch.split('x')))
        model=models[patch]
        if key not in exact:exact[key]=model.doublon(model.exact(tau),(1,1))
        rmin=max(1,ceil(opts['minimum_steps_per_tau']*tau))
        measurements=[]
        for r in (rmin,2*rmin):
            psi=model.trotter(tau,r,algorithm)
            error=abs(model.doublon(psi,(1,1))-exact[key])
            measurements.append(dict(r=r,error=error,coefficient=error*r*r))
        rows.append(dict(patch=patch,tau=tau,algorithm=algorithm,minimum_steps=rmin,measurements=measurements,
            coefficient=max([x['coefficient'] for x in measurements]+[v['coefficient']])))
    (OUT/'finite_step_guard.json').write_text(json.dumps(dict(rows=rows,
        scope='Measured coefficients at the minimum allowed r and twice r, combined with the asymptotic pair; interpolation in r is checked at selected optimized points'),indent=2)+'\n')
    print('Wrote',len(rows),'finite-step guard cases')


def supplementary_checks():
    c=assumptions();opts=settings();rows=json.loads((OUT/'costs.json').read_text())
    chosen=[v for v in rows if v['scenario']=='central' and v['L'] in (4,6) and v['family'] in ('tau_1','tau_4','tau_8')]
    model=SplitSector(3,3);exact={};physics=[];done=set()
    for v in chosen:
        key=(v['tau'],v['calibration_algorithm'],v['r'])
        if key in done:continue
        done.add(key);tau,algorithm,r=key
        if tau not in exact:exact[tau]=model.doublon(model.exact(tau),(1,1))
        error=abs(model.doublon(model.trotter(tau,r,algorithm),(1,1))-exact[tau])
        predicted=v['coefficient']/r**2
        physics.append(dict(tau=tau,algorithm=algorithm,r=r,measured_error=error,model_error=predicted,passes=error<=predicted))
        assert error<=predicted,(key,error,predicted)
    fine={**opts,
        'trotter_fraction_grid':sorted(set(opts['trotter_fraction_grid']+[.5*(a+b) for a,b in zip(opts['trotter_fraction_grid'],opts['trotter_fraction_grid'][1:])])),
        'logical_fraction_grid':[.0025,.005,.01,.02,.035,.05,.075,.1,.15,.2,.3],
        'synthesis_fraction_grid':[.005,.01,.02,.03,.05,.075,.1,.15,.2]}
    refinement=[]
    for v in chosen:
        if v['tau']!=1:continue
        new=optimize(v['L'],v['tau'],v['arm'],c,v['absolute_error'],fine)
        refinement.append(dict(L=v['L'],arm=v['arm'],absolute_error=v['absolute_error'],
            reported_log10_Q=v['log10_Q_min'],finer_grid_log10_Q=new['log10_Q_min'],
            ratio=10**(v['log10_Q_min']-new['log10_Q_min']),fine_r=new['r']))
    data=dict(selected_step_checks=physics,grid_sensitivity=refinement,
              max_grid_improvement=max(v['ratio'] for v in refinement))
    (OUT/'supplementary_checks.json').write_text(json.dumps(data,indent=2)+'\n')
    return data


def audit():
    c=assumptions();rows=json.loads((OUT/'costs.json').read_text());opts=settings()
    feasible=0
    for v in rows:
        if v['log10_Q_min'] is None:continue
        feasible+=1;b=v['budget_fractions'];eps=v['absolute_error']
        assert min(b.values())>=0 and b['stat']>0
        assert isclose(sum(b.values()),1,abs_tol=2e-12)
        assert isclose(v['coefficient']/v['r']**2,eps*b['trotter'],rel_tol=1e-12)
        changes=c['scenarios'].get(v['scenario'],{})
        p={**c,**changes}
        if v['shot_rule']=='gaussian':
            n=ceil(p['gaussian_z']**2*p['signal_reference']*(1-p['signal_reference'])/(eps*b['stat'])**2)
        else:n=ceil(log(2/p['delta'])/(2*(eps*b['stat'])**2))
        assert isclose(log10(n)+v['chi']/log(10),v['log10_total_shots'],abs_tol=2e-10)
        # Work in logarithms, since intentionally infeasible PEC exposures can be enormous.
        lr=max(0,log(n)+v['chi']+log(v['seconds_per_attempt']/p['week_seconds']))
        if lr<30:lr=log(ceil(exp(lr)))
        expected=(log(v['qubits_per_copy'])+lr)/log(10)
        assert isclose(expected,v['log10_Q_min'],abs_tol=1e-9)
        if v['arm']!='nisq':
            d=v['distance'];sources=v.get('factories',0)
            tiles=(v['qubits_per_copy']-sources*p['source_qubits'])/(2*d*d)
            lf=log(tiles*ceil(v['seconds_per_attempt']/p['cycle_seconds'])*p['surface_A'])+(d+1)/2*log(p['p']/p['surface_threshold'])
            assert lf+v['chi']/2<=log(eps*b['logical'])+1e-12
        if v['arm'] in ('ft','hybrid'):
            assert log(v['T_count']*v['magic_error'])+v['chi']/2<=log(eps*b['magic'])+1e-12
        assert b['residual_bias']==opts['residual_bias_fraction'][v['arm']]
    with (OUT/'capacity.csv').open() as f:
        capacities=list(csv.DictReader(f))
    for row in capacities:
        for key in ('nisq','star','hybrid','ft','pec_two_qubit','compiled_pec'):
            arm=key if key in ('nisq','star','hybrid','ft') else 'nisq'
            scenario='central' if key==arm else key
            subset=[v for v in rows if v['arm']==arm and v['scenario']==scenario and v['family']==row['family'] and v['absolute_error']==float(row['absolute_error'])]
            assert int(row[key])==capacity(subset,int(row['Q']))
    calibration_rows=[]
    for path in (OUT/'calibration.json',OUT/'rectangle_checks/calibration.json'):
        data=json.loads(path.read_text())
        assert data['source_sha256']==hashlib.sha256((ROOT/'python/time_calibration.py').read_bytes()).hexdigest()
        for v in data['rows']:
            assert v['converged'] and v['relative_coefficient_spread']<=.05
            assert v['exact_norm_error']<1e-8
            assert all(x['norm_error']<1e-8 for x in v['measurements'])
        calibration_rows+=data['rows']
    checks=json.loads((OUT/'supplementary_checks.json').read_text())
    guard=json.loads((OUT/'finite_step_guard.json').read_text())['rows']
    assert len(guard)==72
    assert all(v['coefficient']>=max(x['coefficient'] for x in v['measurements']) for v in guard)
    assert all(v['passes'] for v in checks['selected_step_checks'])
    report=dict(scenario_points=len(rows),feasible_model_points=feasible,capacity_rows=len(capacities),
        calibration_cases=len(calibration_rows),finite_step_guard_cases=len(guard),selected_step_checks=len(checks['selected_step_checks']),
        maximum_finer_grid_improvement=checks['max_grid_improvement'],
        status='PASS: model budget arithmetic, finite-patch convergence, selected-step checks, and capacity inversion; bulk transfer and residual bias uncalibrated')
    (OUT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--supplementary',action='store_true');p.add_argument('--guard',action='store_true');args=p.parse_args()
    if args.guard:step_guard()
    elif args.supplementary:print(json.dumps(supplementary_checks(),indent=2))
    else:audit()
