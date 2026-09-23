"""Joint error allocation and source-menu optimization for discussion v3.

The Trotter and magic costs are charged at their model values. Statistics get
what remains after logical, synthesis, and residual-bias reserves. All other
costs use discussion_model, so fixed and optimized estimates share accounting.
"""
import json
from functools import lru_cache
from math import ceil, exp, log, sqrt

from discussion_model import estimate
from model import workload
from project_paths import ROOT
from resource_extensions import analog_exposure, source_error, synthesis_plan


def settings():
    return json.loads((ROOT/'configs/budget_optimization.json').read_text())


@lru_cache(maxsize=1)
def calibration():
    return json.loads((ROOT/'results/discussion_v3/calibration.json').read_text())


@lru_cache(maxsize=1)
def finite_step_guard():
    return json.loads((ROOT/'results/discussion_v3/finite_step_guard.json').read_text())


def coefficient(tau,algorithm,options=None):
    options=options or settings()
    values=[v for v in calibration()['rows'] if v['patch']==options['coefficient_patch']
            and v['algorithm']==algorithm and abs(v['tau']-tau)<1e-10]
    if len(values)!=1:raise ValueError(f'No measured coefficient at tau={tau}, algorithm={algorithm}')
    if not values[0]['converged']:raise ValueError('Unconverged coefficient')
    policy=options.get('coefficient_policy','pointwise')
    if policy=='running_max':
        values=[v for v in calibration()['rows'] if v['patch']==options['coefficient_patch']
                and v['algorithm']==algorithm and v['tau']<=tau and v['converged']]
    elif policy!='pointwise':raise ValueError(policy)
    coefficients=[v['coefficient'] for v in values]
    if options.get('finite_step_guard',False):
        coefficients.extend(v['coefficient'] for v in finite_step_guard()['rows'] if v['algorithm']==algorithm and v['patch']==options['coefficient_patch']
            and (v['tau']<=tau if policy=='running_max' else abs(v['tau']-tau)<1e-10))
    return options['coefficient_margin']*max(coefficients)


def algorithm_for(arm):
    return 'termwise' if arm=='star' else 'plaquette'


def step_candidates(C,tau,epsilon,options):
    minimum=max(1,ceil(options['minimum_steps_per_tau']*tau))
    return sorted({max(minimum,ceil(sqrt(C/(f*epsilon)))) for f in options['trotter_fraction_grid']})


def optimize(L,tau,arm,c,epsilon,options=None,C_override=None,refine=True):
    options=options or settings();algorithm=c.get('workload_algorithm',algorithm_for(arm))
    if L==2 and options.get('exact_2x2_coefficients',True):options={**options,'coefficient_patch':'2x2'}
    C=coefficient(tau,algorithm,options) if C_override is None else C_override
    bias=options['residual_bias_fraction'][arm]
    base=dict(L=L,tau=tau,arm=arm,absolute_error=epsilon,calibration_algorithm=algorithm,
              coefficient=C,coefficient_patch=options['coefficient_patch'],
              coefficient_policy=options.get('coefficient_policy','pointwise') if C_override is None else 'explicit_override',
              optimization='declared grid')
    best=None;count=0
    def evaluate_steps(r):
        nonlocal best,count
        trotter=C/(r*r*epsilon)
        if trotter+bias>=1:return
        w=workload(L,tau,r,algorithm,c['u_over_j'])
        if arm=='nisq':
            chi=0. # The estimator cost itself is computed by estimate below.
        elif arm=='ft':chi=0.
        else:chi=sum(analog_exposure(w,c))
        synthesis_grid=options['synthesis_fraction_grid'] if arm=='ft' else [0.]
        for synthesis in synthesis_grid:
            nt=w.fixed_T
            if arm=='ft':nt=synthesis_plan(w,epsilon,{**c,'synthesis_fraction':synthesis})['T_count']
            for cleanup in options['cleanup_choices'] if arm in ('ft','hybrid') else [0]:
                magic=0.
                if arm in ('ft','hybrid'):
                    q0=source_error(c)
                    q=q0 if cleanup==0 else 35*q0**3
                    lm=log(nt*q/epsilon)+chi/2
                    if lm>=0:continue
                    magic=max(1e-300,exp(lm))*(1+1e-12)
                for logical in options['logical_fraction_grid'] if arm!='nisq' else [0.]:
                    stat=1-trotter-bias-synthesis-magic-logical
                    if stat<=0:continue
                    candidate={**c,'stat_fraction':stat,'synthesis_fraction':synthesis,
                               'magic_fraction':magic,'logical_fraction':logical,
                               'forced_cleanup':cleanup}
                    row=estimate(L,tau,arm,candidate,epsilon,steps=r);count+=1
                    if row['log10_Q_min'] is None:continue
                    row.update(base)
                    row['assumed_trotter_error']=C/r**2
                    row['budget_fractions']=dict(stat=stat,trotter=trotter,logical=logical,
                        synthesis=synthesis,magic=magic,residual_bias=bias)
                    exposure=analog_exposure(w,c) if arm in ('star','hybrid') else (0.,0.)
                    row['chi_angle'],row['chi_floor']=exposure
                    if c.get('star_law')=='original_count' and arm in ('star','hybrid'):
                        row['chi_count']=row['chi_angle'];row['chi_angle']=0.
                    # Break integer-replica ties by smaller expected space-time,
                    # then by more statistical allowance, to avoid gratuitous work.
                    key=(row['log10_Q_min'],log(row['qubits_per_copy'])+log(row['seconds_per_attempt'])+row['log10_total_shots']*log(10),-stat)
                    if best is None or key<best[0]:best=(key,row)
    steps=options.get('explicit_step_candidates')
    if steps is None:steps=step_candidates(C,tau,epsilon,options)
    elif not steps or any(not isinstance(r,int) or r<ceil(options['minimum_steps_per_tau']*tau) for r in steps):
        raise ValueError('Explicit step menu must contain admissible positive integers')
    for r in steps:evaluate_steps(r)
    if best is not None and refine and 'explicit_step_candidates' not in options:
        center=best[1]['r'];minimum=max(1,ceil(options['minimum_steps_per_tau']*tau))
        for r in range(max(minimum,center-2),center+3):
            if r not in steps:evaluate_steps(r)
    if best is None:return {**base,'log10_Q_min':None,'candidates_evaluated':count}
    return {**best[1],'candidates_evaluated':count}
