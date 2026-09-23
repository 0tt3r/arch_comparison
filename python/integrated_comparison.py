"""Integrate component evidence and expose NISQ accuracy/cost requirements.

Historical calculations are read-only. STAR uses the existing exact restricted
instrument, including noisy PCEC, RUS retries and declared detector rejection.
Growth is an explicit composition scenario. No experimental estimator is fitted.
"""
import argparse
import csv
import hashlib
import json
from functools import lru_cache
from math import ceil, exp, expm1, floor, isclose, log, log10, log1p, sqrt

from budget_optimization import settings as budget_settings, step_candidates
from discussion_model import assumptions, capacity, data_block_tiles, estimate, sampling_cost
from model import log_gamma2_dephasing, workload
from physics_sweep import planning_coefficient, solve
from project_paths import ROOT
from resource_extensions import source_error
from star_channel_study import growth_error, preparation, raw_rus_details, rus_channel

OUT = ROOT / 'results/integrated_v1'
ARMS = ('star', 'hybrid', 'ft')


def read(path):
    return json.loads(path.read_text())


@lru_cache(maxsize=1)
def settings():
    return read(ROOT / 'configs/integrated_comparison.json')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def write_csv(name, rows):
    columns = list(dict.fromkeys(k for v in rows for k in v))
    with (OUT / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def configuration(p, scenario='component'):
    cfg = settings()
    return {**assumptions(), 'p': p, 'source_p_policy': 'held_at_anchor',
            'factory_choices': cfg['factory_choices'], **cfg['scenarios'][scenario]}


@lru_cache(maxsize=1)
def growth_pairs():
    return read(ROOT / 'results/frontier_v1/growth_summary.json')['pairs']


def growth_probability(seed, p, scenario, target=25):
    if seed == target:
        return 0.
    if scenario == 'zero_growth':
        if seed != 9:
            raise ValueError('Only d=9 has a zero-growth sensitivity; d=5 is measured noisy')
        return 0.
    return next(v['pauli_total_upper'] for v in growth_pairs()
                if (v['di'], v['df'], v['p']) == (seed, target, p))


CHANNELS = {}


def channel_contract():
    return {'source': sha(ROOT / 'python/star_channel_study.py'),
            'growth': sha(ROOT / 'results/frontier_v1/growth_summary.json'),
            'cap': settings()['rus_cap'], 'targets': settings()['growth_targets']}


def load_channels():
    path = OUT / 'channels.json'
    if path.exists():
        data = read(path)
        if data['contract'] != channel_contract():
            raise ValueError('Channel cache contract changed; use a new result version')
        CHANNELS.update(data['rows'])


def save_channels():
    write('channels.json', {'contract': channel_contract(), 'rows': CHANNELS})


def rotation_channel(theta, seed, p, scenario, target=25):
    """Re-evaluate the existing generator-fault instrument at actual angles.

    The old growth model uses a d=9 anchor internally. Invert that parameter
    conversion so every growth call receives exactly the saved component bound.
    This does not fit or assume a power law for unmeasured injection distances.
    """
    cfg = settings()
    g = growth_probability(seed, p, scenario, target)
    detector = cfg['scenarios'][scenario]['detector_factor']
    key = '|'.join(map(str, (theta, seed, p, g, detector, target)))
    if key not in CHANNELS:
        scale = 0. if seed == target else growth_error(seed, p, 1e-10, target) / 1e-10
        anchor = g / scale if scale else 0.
        if seed != target:
            assert isclose(growth_error(seed, p, anchor, target), g, rel_tol=1e-12, abs_tol=1e-20)
        v = rus_channel(theta, seed, p, mode='noisy_pcec', floor=anchor,
                        data_distance=target, detector_factor=detector, cap=cfg['rus_cap'])
        CHANNELS[key] = dict(theta=theta, seed=seed, p=p, growth_probability=g,
                             detector_factor=detector, **v)
    return CHANNELS[key]


def component_exposure(w, p, seed, scenario, target=25):
    chi = log_norm = seconds = tail = 0.
    for angle, count in w.dynamic_angles:
        v = rotation_channel(angle, seed, p, scenario, target)
        q = v['q_dephasing']
        chi += count * log_gamma2_dephasing(q)
        log_norm += count * log1p(abs(v['coherent_component']) / (1 - 2*q))
        seconds += count * v['mean_seconds']
        tail += count * v['uncorrected_rus_tail']
    return dict(chi=chi, coherent_log_norm=log_norm,
                coherent_bias_bound=.5*expm1(min(log_norm, 700.)),
                rus_tail_bound=tail, serial_analog_seconds=seconds,
                growth_probability=growth_probability(seed, p, scenario, target),
                growth_target=target)


def cost_key(row):
    return (row['log10_Q_min'], log(row['qubits_per_copy']) +
            log(row['seconds_per_attempt']) + log(10)*row['log10_total_shots'])


def encoded_candidate(w, arm, c, epsilon, component, fractions, cleanup, layout):
    """Whole-circuit ledger with joint analog service and cultivated T supply.

    Injection/growth references end at d=13 or 25. Larger data distances use a
    declared protected handoff, charged in time and logical volume. The growth
    census is not represented as a measurement at those larger distances.
    """
    cfg = settings()
    target = component['growth_target']
    chi = component['chi']
    nt = w.fixed_T if arm == 'hybrid' else 0
    local = {**c, **cfg['layouts'][layout]}
    data_tiles = data_block_tiles(w.L, local)
    best = None
    for factories in local['factory_choices'] if nt else [0]:
        tiles = data_tiles + factories*(1 + 12*cleanup)
        minimum_d = max(target, local['source_end_distance'] if nt else 3)
        for d in range(minimum_d, local['max_distance'] + 1, 2):
            # Uniform stretching is an explicit timing sensitivity for the
            # larger-distance receiver, not an inferred higher-d instrument.
            analog = component['serial_analog_seconds'] * d/target / local['parallel_lanes']
            handoff = (2*d*local['cycle_seconds']*w.rotations/local['parallel_lanes']
                       if d > target else 0.)
            supply = consume = 0.
            if nt:
                raw = local['source_cycles_per_T']*local['source_cycle_multiplier'] + 2*d
                service = ((15*raw + 11*d)/(1 - 15*source_error(local)) if cleanup else raw)
                supply = ceil(nt/factories)*service*local['cycle_seconds']
                consume = local['ft_clocks_per_T']*nt*d*local['cycle_seconds']/local['T_consumption_ports']
            seconds = (2*d*local['cycle_seconds'] + max(analog + handoff + consume, supply)
                       + (w.rotations/local['parallel_lanes'] + nt/local['T_consumption_ports'])*local['feedforward_seconds']
                       + local['reset_seconds'] + local['measurement_seconds'])
            cycles = ceil(seconds/local['cycle_seconds'])
            log_fault = (log(tiles*cycles*local['surface_A'])
                         + (d+1)/2*log(local['p']/local['surface_threshold']))
            if log_fault <= log(epsilon*fractions['logical']) - chi/2:
                break
        else:
            continue
        if seconds > local['week_seconds']:
            continue
        qubits = 2*tiles*d*d + factories*local['source_qubits']
        row = sampling_cost(qubits, seconds, chi, {**local, 'stat_fraction': fractions['stat']}, epsilon)
        row.update(distance=d, injection_distance=component['seed'],
                   factories=factories, cleanup_levels=cleanup, T_count=nt,
                   data_tiles=data_tiles, protected_tiles=tiles, layout=layout,
                   supply_seconds=supply, consumption_seconds=consume,
                   analog_seconds=analog, extra_handoff_seconds=handoff,
                   log_logical_fault_bound=log_fault,
                   **component)
        if best is None or cost_key(row) < cost_key(best):
            best = row
    return best


def optimize_analog(L, tau, arm, epsilon, p, scenario):
    c = configuration(p, scenario)
    options = budget_settings()
    algorithm = 'termwise' if arm == 'star' else 'plaquette'
    C = planning_coefficient(L, tau, algorithm)
    seeds = [9, 25] if scenario == 'zero_growth' else settings()['injection_seeds']
    candidates = []
    visited = set()

    def at_step(r):
        if r in visited:
            return
        visited.add(r)
        w = workload(L, tau, r, algorithm, c['u_over_j'])
        trotter = C/r**2/epsilon
        for seed,target in [(seed,target) for seed in seeds for target in settings()['growth_targets'] if seed <= target]:
            part = dict(component_exposure(w, p, seed, scenario, target), seed=seed)
            residual = options['residual_bias_fraction'][arm]
            instrument = (part['coherent_bias_bound'] + part['rus_tail_bound'])/epsilon
            for cleanup in options['cleanup_choices'] if arm == 'hybrid' else [0]:
                magic = 0.
                if arm == 'hybrid':
                    q0 = source_error(c)
                    q = 35*q0**3 if cleanup else q0
                    log_magic = log(w.fixed_T*q/epsilon) + part['chi']/2
                    if log_magic >= 0:
                        continue
                    magic = exp(log_magic)*(1 + 1e-12)
                for logical in options['logical_fraction_grid']:
                    stat = 1 - trotter - residual - instrument - magic - logical
                    if stat <= 0:
                        continue
                    fractions = dict(trotter=trotter, residual_bias=residual, instrument=instrument,
                                     magic=magic, logical=logical, stat=stat, synthesis=0.)
                    for layout in settings()['layouts'] if arm == 'hybrid' else ['compact']:
                        row = encoded_candidate(w, arm, c, epsilon, part, fractions, cleanup, layout)
                        if row is not None:
                            row.update(L=L, tau=tau, arm=arm, absolute_error=epsilon, p=p,
                                       scenario=scenario, algorithm=algorithm, r=r, rotations=w.rotations,
                                       theta=w.theta, coefficient=C, assumed_trotter_error=C/r**2,
                                       budget_fractions=fractions)
                            candidates.append(row)
    for r in step_candidates(C, tau, epsilon, options):
        at_step(r)
    if candidates:
        center = min(candidates, key=cost_key)['r']
        for r in range(max(ceil(options['minimum_steps_per_tau']*tau), center-2), center+3):
            at_step(r)
    if not candidates:
        return dict(L=L, tau=tau, arm=arm, absolute_error=epsilon, p=p, scenario=scenario,
                    log10_Q_min=None, status='No feasible candidate in the declared menu')
    return dict(min(candidates, key=cost_key), candidate_count=len(candidates),
                status='Conditional component composition; aggregate layout')


def optimize_ft(L, tau, epsilon, p, scenario):
    cfg = settings()
    candidates = []
    for layout, block in cfg['layouts'].items():
        for batch in cfg['hwp_batches']:
            if batch > L*L:
                continue
            synthesis = 'individual' if batch == 1 else 'hamming_weight'
            c = {**configuration(p, scenario), **block, 'rotation_synthesis': synthesis,
                 'hwp_batch_size': batch}
            row = solve(L, tau, 'ft', epsilon, c, scenario=scenario, layout=layout,
                        synthesis=synthesis, batch=batch)
            if row['log10_Q_min'] is not None:
                candidates.append(row)
    return dict(min(candidates, key=cost_key), candidate_count=len(candidates),
                status='Synthesis counts checked; aggregate delivery assumption')


def capacity_table(rows):
    cfg = settings()
    table = []
    for p in cfg['physical_errors']:
        for eps in cfg['precisions']:
            for tau in cfg['times']:
                for Q in cfg['qubit_budgets']:
                    row = dict(p=p, epsilon=eps, tau=tau, Q=Q)
                    for arm in ARMS:
                        values = []
                        for scenario in cfg['scenarios']:
                            selected = [v for v in rows if (v['p'], v['absolute_error'], v['tau'], v['arm'], v['scenario'])
                                        == (p, eps, tau, arm, scenario)]
                            value = capacity(selected, Q)
                            row[arm + '_' + scenario] = value
                            values.append(value)
                        row[arm + '_min'], row[arm + '_max'] = min(values), max(values)
                    row['robust_capacity_leader'] = ','.join(arm for arm in ARMS
                        if row[arm+'_min'] > max(row[a+'_max'] for a in ARMS if a != arm)) or 'unresolved'
                    row['at_size_ceiling'] = max(row[a+'_max'] for a in ARMS) == max(cfg['sizes'])
                    table.append(row)
    return table


def planned_runtime(row, Q):
    """Runtime of the selected implementation; not a fresh fixed-Q optimizer."""
    copies = int(Q//row['qubits_per_copy'])
    if copies == 0 or row['log10_Q_min'] is None:
        return None
    log_time = log(row['seconds_per_attempt']) + log(10)*row['log10_total_shots'] - log(copies)
    # A capped value cannot affect a one-week requirement. It prevents float
    # overflow in astronomically expensive, but algebraically valid scenarios.
    return max(row['seconds_per_attempt'], exp(min(700., log_time)))


def threshold(available_attempts, remaining_error, variance=.09, z=1.96):
    if available_attempts < 1 or remaining_error <= 0:
        return None
    return available_attempts*remaining_error**2/(z*z*variance)


def offset_limit(available_attempts, error_after_trotter, F, variance=.09, z=1.96):
    if available_attempts < 1 or F <= 0:
        return None
    value = error_after_trotter - z*sqrt(F*variance/available_attempts)
    return value if value >= 0 else None


def nisq_attempt(L, r, c):
    q = ceil(sqrt(2*L*L + 2))**2
    depth = c['nisq_anchor_depth']*(r/c['nisq_anchor_steps'])*(L/c['nisq_anchor_L'])
    seconds = c['reset_seconds'] + depth*c['gate_seconds'] + c['measurement_seconds']
    gates = c['nisq_anchor_gates']*(r/c['nisq_anchor_steps'])*(L/c['nisq_anchor_L'])**3
    return q, seconds, gates


def break_even(rows):
    """Requirements on an unknown estimator, without assigning it a bias.

    F is variance per attempted circuit divided by D*(1-D*); rejection is
    already included in F. Training shots are charged separately and are
    assumed to have the same duration as evaluation shots. The accuracy
    statement is a Gaussian planning model, not an estimator certificate.
    """
    cfg = settings()
    table = []
    for p in cfg['physical_errors']:
      c = configuration(p)
      for eps in cfg['precisions']:
       for tau in cfg['times']:
        for L in cfg['break_even_sizes']:
         C = planning_coefficient(L, tau, 'plaquette')
         max_r = ceil(sqrt(8*C/(eps*(1-max(cfg['bias_fractions'])))))
         menu = [(r, C/r**2, *nisq_attempt(L, r, c))
                 for r in range(max(ceil(8*tau), ceil(sqrt(C/eps))), max_r+1)]
         for Q in cfg['break_even_qubit_budgets']:
          refs = [v for v in rows if v['arm'] in ARMS and v['log10_Q_min'] is not None
                  and (v['L'],v['tau'],v['absolute_error'],v['p']) == (L,tau,eps,p)]
          for scenario in cfg['scenarios']:
           feasible = [(planned_runtime(v,Q),v) for v in refs if v['scenario']==scenario]
           feasible = [(t,v) for t,v in feasible if t is not None]
           reference = min(feasible, key=lambda x:x[0]) if feasible else None
           for speedup in cfg['speedups']:
            if reference is None and speedup != 1:
                continue
            # Beating an encoded implementation and meeting the week cutoff
            # are both required. When no encoded copy fits, only the latter applies.
            allowed = min(c['week_seconds'], reference[0]/speedup) if reference else c['week_seconds']
            for training in cfg['training_attempts']:
             available = [(r,e,q,t,g, max(0,int(Q//q)*floor(allowed/t)-training)) for r,e,q,t,g in menu]
             base = dict(p=p,epsilon=eps,tau=tau,L=L,Q=Q,scenario=scenario,speedup=speedup,
                         allowed_seconds=allowed,reference_arm=reference[1]['arm'] if reference else 'none_fits',
                         encoded_seconds=reference[0] if reference else None,training_attempts=training)
             for f in cfg['bias_fractions']:
                candidates = [(threshold(M,eps*(1-f)-e),r,e,q,t,g,M) for r,e,q,t,g,M in available]
                candidates = [v for v in candidates if v[0] is not None]
                if not candidates:
                    table.append(dict(base,kind='F_max',bias_fraction=f,F_max=None,status='no_budget'))
                    continue
                F,r,e,q,t,g,M = max(candidates,key=lambda x:x[0])
                table.append(dict(base,kind='F_max',bias_fraction=f,bias=eps*f,F_max=F,r=r,
                                  trotter_error=e,qubits_per_copy=q,attempt_seconds=t,two_qubit_gates=g,
                                  evaluation_attempts=M,kappa_max=log(F)/(p*g),status='requirement_only'))
             for F in cfg['variance_multipliers']:
                candidates = [(offset_limit(M,eps-e,F),r,e,q,t,g,M) for r,e,q,t,g,M in available]
                candidates = [v for v in candidates if v[0] is not None]
                if not candidates:
                    table.append(dict(base,kind='bias_max',F=F,bias_max=None,status='no_budget'))
                    continue
                b,r,e,q,t,g,M = max(candidates,key=lambda x:x[0])
                table.append(dict(base,kind='bias_max',F=F,bias_max=b,r=r,trotter_error=e,
                                  qubits_per_copy=q,attempt_seconds=t,two_qubit_gates=g,
                                  evaluation_attempts=M,status='requirement_only'))
    return table


def input_hashes():
    names = ['python/integrated_comparison.py', 'configs/integrated_comparison.json',
             'python/star_channel_study.py','python/physics_sweep.py','python/discussion_model.py',
             'python/budget_optimization.py','python/resource_extensions.py','python/model.py',
             'python/dynamics_followup_analysis.py','configs/assumptions.json',
             'configs/discussion_model.json','configs/budget_optimization.json',
             'results/discussion_v3/calibration.json','results/discussion_v3/finite_step_guard.json',
             'results/dynamics_followup_v1/measurements.json','results/frontier_v1/growth_summary.json',
             'results/implementation_v1/star_channels.json','results/physics_sweep_v1/costs.json']
    return {name:sha(ROOT/name) for name in names}


def validate(rows, thresholds):
    for v in rows:
        if v['log10_Q_min'] is None:
            continue
        b = v['budget_fractions']
        assert isclose(sum(b.values()),1,abs_tol=1e-12) and min(b.values()) >= 0
        assert isclose(v['assumed_trotter_error'],v['coefficient']/v['r']**2,rel_tol=1e-12)
        if v['arm'] in ('star','hybrid'):
            w = workload(v['L'],v['tau'],v['r'],v['algorithm'])
            part = component_exposure(w,v['p'],v['injection_distance'],v['scenario'],v['growth_target'])
            for k in ('chi','coherent_bias_bound','rus_tail_bound','serial_analog_seconds'):
                assert isclose(v[k],part[k],rel_tol=1e-12,abs_tol=1e-18),k
            assert v['distance'] >= v['growth_target']
            assert v['log_logical_fault_bound'] + v['chi']/2 <= log(v['absolute_error']*b['logical']) + 1e-10
        expected = sampling_cost(v['qubits_per_copy'],v['seconds_per_attempt'],v['chi'],
                                 {**configuration(v['p'],v['scenario']), 'stat_fraction':b['stat']},v['absolute_error'])
        assert isclose(expected['log10_Q_min'],v['log10_Q_min'],abs_tol=1e-11)
    for v in thresholds:
        if v['status']=='no_budget':
            continue
        if v['kind']=='F_max':
            error = v['bias']+v['trotter_error']+1.96*sqrt(v['F_max']*.09/v['evaluation_attempts'])
        else:
            error = v['bias_max']+v['trotter_error']+1.96*sqrt(v['F']*.09/v['evaluation_attempts'])
        assert isclose(error,v['epsilon'],rel_tol=1e-11)
        assert v['evaluation_attempts']+v['training_attempts'] <= int(v['Q']//v['qubits_per_copy'])*floor(v['allowed_seconds']/v['attempt_seconds'])
    return dict(status='PASS: component composition, budgets, replicas and break-even identities',
                resource_rows=len(rows),threshold_rows=len(thresholds),channel_points=len(CHANNELS),
                scope=settings()['scope'])


def run(check=False):
    load_channels()
    if check:
        provenance = read(OUT/'inputs.json')
        assert provenance == input_hashes(), 'Integrated inputs changed'
        rows, thresholds = read(OUT/'costs.json'), read(OUT/'break_even.json')
        assert read(OUT/'capacity.json') == capacity_table(rows)
        assert read(OUT/'validation.json') == validate(rows,thresholds)
        print('Integrated accounting verified')
        return
    contract = input_hashes()
    if (OUT/'inputs.json').exists() and read(OUT/'inputs.json') != contract:
        raise ValueError('Inputs changed; remove only integrated costs.json and inputs.json to recompute')
    write('inputs.json',contract)
    rows = read(OUT/'costs.json') if (OUT/'costs.json').exists() else []
    done = {(v['L'],v['tau'],v['absolute_error'],v['p'],v['arm'],v['scenario']) for v in rows}
    cfg = settings()
    for tau in cfg['times']:
     for p in cfg['physical_errors']:
      for eps in cfg['precisions']:
       for L in cfg['sizes']:
        for scenario in cfg['scenarios']:
         for arm in ARMS:
            key = (L,tau,eps,p,arm,scenario)
            if key in done:
                continue
            if arm=='ft':
                # Its two fast-supply scenarios are identical by definition.
                prior = next((v for v in rows if (v['L'],v['tau'],v['absolute_error'],v['p'],v['arm'],v['scenario'])
                              == (L,tau,eps,p,'ft','component')),None)
                row = dict(prior,scenario=scenario) if scenario=='zero_growth' and prior else optimize_ft(L,tau,eps,p,scenario)
            else:
                row = optimize_analog(L,tau,arm,eps,p,scenario)
            rows.append(row)
            write('costs.json',rows)
            save_channels()
            # Persistent scalar channel records suffice between work items;
            # the exact instrument's internal mixtures can be much larger.
            preparation.cache_clear()
            raw_rus_details.cache_clear()
        print('integrated',tau,p,eps,L,len(rows),'rows',flush=True)
    thresholds = break_even(rows)
    write('break_even.json',thresholds)
    write_csv('break_even.csv',thresholds)
    caps = capacity_table(rows)
    write('capacity.json',caps)
    write_csv('capacity.csv',caps)
    write('validation.json',validate(rows,thresholds))
    save_channels()
    print('Completed',len(rows),'resource cases;',len(thresholds),'NISQ requirements',flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    run(parser.parse_args().check)
