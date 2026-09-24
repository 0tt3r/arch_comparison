"""Bounded optimization round: which choices change the capacity leader.

Five switches from OPUS55_REVIEW.md section 0 are tested alone and together:
raw_rus      raw repeat-until-success STAR channels beside noisy PCEC;
lanes        optimized analog lanes, charged per lane and capped by concurrency;
ft_delivery  joint T-consumption ports, factories and Hamming-weight batches,
             without STAR lane tiles in the FT block;
distances    separate source and data distances with a charged transfer;
richardson   two-branch step extrapolation with branch-specific costs.
The v1_menu configuration reproduces results/integrated_v1 at tau=1; --check
enforces that replay. Channels and tile charges are declared models, not
compiled schedules or new gadget calibrations.
"""
import argparse
import csv
import hashlib
import json
import multiprocessing as mp
from dataclasses import replace
from functools import lru_cache
from math import ceil, exp, expm1, isclose, log, log10, log1p, pi, sqrt

import numpy as np

import integrated_comparison as v1
from budget_optimization import coefficient, optimize, settings as budget_settings, step_candidates
from discussion_model import capacity, data_block_tiles, estimate, ideal_shots, sampling_cost
from dynamics_followup_analysis import extended_coefficient
from model import log_gamma2_dephasing, workload
from physics_sweep import configuration as sweep_configuration, planning_coefficient
from project_paths import ROOT
from resource_extensions import plaquette_rotation_stages, source_error, synthesis_plan
from star_channel_study import preparation, raw_rus_details, rus_channel

ROUNDS = {'v1': ('configs/optimization_round.json', 'results/optimization_round_v1'),
          'v2': ('configs/optimization_round_v2.json', 'results/optimization_round_v2'),
          'v3': ('configs/optimization_round_v3.json', 'results/optimization_round_v3')}
CONFIG, OUT = (ROOT/x for x in ROUNDS['v1'])
ENCODED = ('star', 'hybrid', 'ft')
NISQ = ('nisq_effective', 'nisq_pec')
MAX_SEARCH_DISTANCE = 301


def read(path):
    return json.loads(path.read_text())


@lru_cache(maxsize=1)
def settings():
    return read(CONFIG)


def use_round(name):
    """Select a result version: v1 (reference surface law) or v2 (calibrated)."""
    global CONFIG, OUT
    CONFIG, OUT = (ROOT/x for x in ROUNDS[name])
    settings.cache_clear()
    surface_law.cache_clear()


@lru_cache(maxsize=None)
def surface_law(p):
    """Calibrated per-p law, or None for the reference A=0.1, p_th=0.01."""
    law = settings().get('surface_law')
    if not law:
        return None
    fit = read(ROOT/law['source'])['fits'][str(p)]
    return dict(surface_A=fit[law['amplitude']], surface_threshold=fit['effective_p_threshold'])


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT/name
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temp.replace(path)


def write_csv(name, rows):
    columns = list(dict.fromkeys(k for v in rows for k in v))
    with (OUT/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def switches(name):
    base = dict(rus_modes=['noisy_pcec'], optimize_lanes=False, ft_delivery=False,
                separate_distances=False, richardson=False)
    return {**base, **settings()['configurations'][name]}


def family_tau(family, L):
    f = settings()['families'][family]
    return f['tau'] if f['kind'] == 'fixed' else f['scale']*L


# ---------------------------------------------------------------- Trotter data

@lru_cache(maxsize=1)
def _measurements():
    return read(ROOT/'results/dynamics_followup_v1/measurements.json')


def trotter_coefficient(L, tau, algorithm):
    """integrated_v1 coefficient where it exists; otherwise the same envelope rule.

    For times outside the 4x4 checks this is the running maximum of all measured
    patches through tau, exactly the extended_coefficient rule without its
    restriction to tau in {1, 2, 4}.
    """
    if L == 2 or tau in read(ROOT/'configs/patch_sensitivity.json')['times']:
        return planning_coefficient(L, tau, algorithm)
    extra = [v['coefficient'] for v in _measurements() if v['tau'] <= tau and v['algorithm'] == algorithm]
    rectangles = read(ROOT/'results/discussion_v3/rectangle_checks/calibration.json')['rows']
    extra.extend(abs(m['coefficient']) for v in rectangles
                 if v['tau'] <= tau and v['algorithm'] == algorithm for m in v['measurements'])
    return max(coefficient(tau, algorithm), 1.1*max(extra))


def richardson_pairs(algorithm, patches=None):
    """Measured |(4E(2r)-E(r))/3| r^4 from every same-patch (r, 2r) signed pair.

    Sources: 4x4/4x3 dynamics, 2x2/3x3 calibration and 4x3 rectangles. The
    finite-step guard stores unsigned errors; its only valid bound,
    (4|E(2r)|+|E(r)|)/3, cannot cancel the leading r^-2 term and is not used.
    """
    groups = {}
    for v in _measurements():
        if v['algorithm'] == algorithm:
            groups.setdefault((v['patch'], v['tau']), {})[v['r']] = v['signed_error']
    for name in ('results/discussion_v3/calibration.json', 'results/discussion_v3/rectangle_checks/calibration.json'):
        for v in read(ROOT/name)['rows']:
            if v['algorithm'] == algorithm:
                for m in v['measurements']:
                    groups.setdefault((v['patch'], v['tau']), {})[m['r']] = m['signed_error']
    pairs = []
    for (patch, tau), errors in groups.items():
        if patches is not None and patch not in patches:
            continue
        for r, e in errors.items():
            if 2*r in errors:
                pairs.append(dict(patch=patch, tau=tau, r=r, K4=abs(4*errors[2*r]-e)/3*r**4))
    return pairs


@lru_cache(maxsize=None)
def richardson_coefficient(L, tau, algorithm):
    """Envelope of measured extrapolation residuals through tau, with margin."""
    patches = {'2x2'} if L == 2 else None
    values = [v['K4'] for v in richardson_pairs(algorithm, patches) if v['tau'] <= tau + 1e-12]
    if not values:
        raise ValueError(f'No measured Richardson pair through tau={tau}')
    return settings()['richardson']['residual_margin']*max(values)


def branch_plan(r, richardson):
    if not richardson:
        return ((r, 1.),)
    c1, c2 = settings()['richardson']['coefficients']
    return ((r, c1), (settings()['richardson']['step_ratio']*r, c2))


def richardson_steps(K4, tau, epsilon, options):
    minimum = max(1, ceil(options['minimum_steps_per_tau']*tau))
    return sorted({max(minimum, ceil((K4/(f*epsilon))**.25)) for f in options['trotter_fraction_grid']})


def trotter_fraction(L, tau, algorithm, r, richardson, epsilon):
    if richardson:
        return richardson_coefficient(L, tau, algorithm)/r**4/epsilon
    return trotter_coefficient(L, tau, algorithm)/r**2/epsilon


# ---------------------------------------------------------------- STAR channels

CHANNELS = {}
NEW_CHANNELS = {}


def channel_contract():
    return {'source': sha(ROOT/'python/star_channel_study.py'),
            'growth': sha(ROOT/'results/frontier_v1/growth_summary.json'),
            'cap': v1.settings()['rus_cap'], 'targets': v1.settings()['growth_targets']}


def load_channels():
    """Seed from integrated_v1 (noisy PCEC keys) and this round's own cache."""
    old = v1.OUT/'channels.json'
    if old.exists():
        data = read(old)
        if data['contract'] == channel_contract():
            CHANNELS.update({k+'|noisy_pcec': v for k, v in data['rows'].items()})
    for path in dict.fromkeys([ROOT/ROUNDS['v1'][1]/'channels.json', OUT/'channels.json']):
        if not path.exists():
            continue
        data = read(path)
        if data['contract'] != channel_contract():
            raise ValueError('Channel cache contract changed; use a new result version')
        CHANNELS.update(data['rows'])


def save_channels():
    write('channels.json', {'contract': channel_contract(), 'rows': CHANNELS})


def rotation_channel(theta, seed, p, scenario, target, mode):
    """integrated_v1.rotation_channel with the RUS correction mode as a key."""
    cfg = v1.settings()
    g = v1.growth_probability(seed, p, scenario, target)
    detector = settings()['scenarios'][scenario]['detector_factor']
    key = '|'.join(map(str, (theta, seed, p, g, detector, target, mode)))
    if key not in CHANNELS:
        from star_channel_study import growth_error
        scale = 0. if seed == target else growth_error(seed, p, 1e-10, target)/1e-10
        anchor = g/scale if scale else 0.
        value = rus_channel(theta, seed, p, mode=mode, floor=anchor, data_distance=target,
                            detector_factor=detector, cap=cfg['rus_cap'])
        CHANNELS[key] = NEW_CHANNELS[key] = dict(theta=theta, seed=seed, p=p, growth_probability=g,
                                                 detector_factor=detector, mode=mode, **value)
    return CHANNELS[key]


def component_exposure(w, p, seed, scenario, target, mode):
    chi = log_norm = seconds = tail = 0.
    for angle, count in w.dynamic_angles:
        v = rotation_channel(angle, seed, p, scenario, target, mode)
        q = v['q_dephasing']
        chi += count*log_gamma2_dephasing(q)
        log_norm += count*log1p(abs(v['coherent_component'])/(1-2*q))
        seconds += count*v['mean_seconds']
        tail += count*v['uncorrected_rus_tail']
    return dict(chi=chi, coherent_log_norm=log_norm, coherent_bias_bound=.5*expm1(min(log_norm, 700.)),
                rus_tail_bound=tail, serial_analog_seconds=seconds, rus_mode=mode, seed=seed,
                growth_probability=v1.growth_probability(seed, p, scenario, target), growth_target=target)


# ---------------------------------------------------------------- shared ledgers

def apply_delivery(c, sw):
    """Declared delivery sensitivities (round v3); absent keys keep the reference values."""
    changes = {k: sw[v] for k, v in (('extra_tiles_per_lane', 'tiles_per_lane'), ('source_error', 'source_error')) if v in sw}
    return {**c, **changes}


def source_menu(c, sw, cleanup_choices):
    """T sources: cultivation with optional 15-to-1 cleanup, plus Litinski (2019) factories at matching p
    (published, or re-optimized under the calibrated surface law by factory_calibration.py).

    tiles: protected tiles charged per factory (delivery, plus cleanup tiles);
    qubits: physical qubits per factory outside the protected block.
    """
    q0 = source_error(c)
    menu = [dict(name='cultivation', cleanup=k, error=35*q0**3 if k else q0, tiles=1+12*k,
                 qubits=c['source_qubits']) for k in cleanup_choices]
    if sw.get('source_menu') in ('litinski', 'calibrated'):
        table = (settings()['litinski_factories'] if sw['source_menu'] == 'litinski'
                 else read(ROOT/settings()['calibrated_factories'])['litinski_factories'])
        for name, error, qubits, cycles, outputs in table.get(repr(c['p']), []):
            menu.append(dict(name=name, cleanup=0, error=error, tiles=1, qubits=qubits,
                             cycles_per_T=cycles/outputs))
    return menu


def source_service(src, c, d, transfer):
    """Code cycles between T outputs of one factory."""
    if src['name'] != 'cultivation':
        return src['cycles_per_T'] + 2*transfer
    raw = c['source_cycles_per_T']*c['source_cycle_multiplier'] + 2*transfer
    return (15*raw+11*d)/(1-15*source_error(c)) if src['cleanup'] else raw


def analog_fixed_T(w):
    """Plaquette workload whose fixed pi/8 rotations are injected by the STAR instrument."""
    return replace(w, dynamic_angles=w.dynamic_angles+((pi/8, w.fixed_T),), fixed_T=0) if w.fixed_T else w


def port_tile_charge(sw):
    return sw.get('tiles_per_port', settings()['extra_tiles_per_port'])


def base_config(p, scenario):
    return {**v1.configuration(p, 'component'), **settings()['scenarios'][scenario], **(surface_law(p) or {})}


def sampling_log_q(qubits, log_P, stat, c, epsilon, rule='gaussian'):
    """Vectorized discussion_model.sampling_cost; log_P = log(variance x seconds)."""
    n0 = ideal_shots({**c, 'stat_fraction': stat}, epsilon, rule)
    log_copies = np.maximum(0., log(n0)+log_P-log(c['week_seconds']))
    small = log_copies < 30
    log_copies = np.where(small, np.log(np.ceil(np.exp(np.minimum(log_copies, 30.)))), log_copies)
    return (np.log(qubits)+log_copies)/log(10)


def combined_log_P(chis, seconds, weights):
    """Optimal shot allocation between branches: (sum |c_i| sqrt(F_i t_i))^2."""
    terms = [log(abs(w)) + .5*(np.asarray(x)+np.log(s)) for w, x, s in zip(weights, chis, seconds)]
    return 2*np.logaddexp.reduce(np.stack(terms), axis=0) if len(terms) > 1 else 2*terms[0]


def final_sampling(qubits, seconds, chis, weights, stat, c, epsilon, rule='gaussian'):
    """Exact sampling_cost call; Richardson enters as an effective exponent."""
    s_max = max(seconds)
    if len(seconds) == 1:
        chi = chis[0]
    else:
        chi = float(combined_log_P(chis, seconds, weights)) - log(s_max)
    return sampling_cost(qubits, s_max, chi, {**c, 'stat_fraction': stat}, epsilon, rule)


def cost_key(row):
    return (row['log10_Q_min'], log(row['qubits_per_copy'])+log(row['seconds_per_attempt'])
            + log(10)*row['log10_total_shots'])


def odd_distances(minimum):
    start = minimum if minimum % 2 else minimum+1
    return np.arange(max(3, start), MAX_SEARCH_DISTANCE+1, 2, dtype=float)


def first_feasible(feasible, axis=-1):
    """Index of the first True along the distance axis, or -1."""
    index = np.argmax(feasible, axis=axis)
    return np.where(np.take_along_axis(feasible, np.expand_dims(index, axis), axis).squeeze(axis), index, -1)


# ---------------------------------------------------------------- STAR and hybrid

def lane_menu(L, algorithm, sw):
    cfg = settings()
    if not sw['optimize_lanes']:
        return [cfg['reference_lanes']]
    width = sw.get('lane_cap', cfg['concurrency_width'][algorithm])
    cap = L if width == 'L' else L*L
    return [n for n in cfg['lane_choices'] if n <= cap] or [1]


def port_menu(L, sw, layout='fast'):
    """Several T-consumption ports need the fast block's exposed boundaries."""
    if not sw['ft_delivery'] or layout == 'compact':
        return [1]
    return [n for n in settings()['port_choices'] if n <= L*L]


def serial_T_layers(w, synthesis, ports):
    """Serial T-consumption layers for synthesized rotations on parallel ports.

    Each synthesized rotation is a serial T sequence. Rotations in one plaquette
    stage are independent; a port consumes one rotation's sequence at a time.
    Returns sum over stages of repeat*ceil(rotations_in_stage/ports), which equals
    the synthesized-rotation count at one port.
    """
    ports = np.asarray(ports, dtype=float)
    batch = synthesis.get('hwp_batch_size', 1)
    layers = 0.
    for width, repeat in plaquette_rotation_stages(w):
        if batch == 1:
            n = width
        else:
            full, rest = divmod(width, batch)
            n = full*batch.bit_length() + (rest.bit_length() if rest else 0)
        layers = layers + repeat*np.ceil(n/ports)
    return layers


def analog_layouts(arm, lanes, sw):
    if arm == 'star' and not sw['optimize_lanes']:
        return ['compact']
    return [name for name in v1.settings()['layouts']
            if name == 'fast' or lanes <= settings()['compact_lane_limit']]


def encoded_candidate(L, arm, c, epsilon, branches, fractions, src, layout, lanes, sw):
    """Vectorized integrated_v1.encoded_candidate over ports, factories and distance.

    branches: sequence of (weight, workload, component). The hardware, including
    distance, is shared; every branch must meet the logical budget divided by
    the extrapolation one-norm, after its own Gamma amplification.
    """
    cfg = v1.settings()
    local = {**c, **cfg['layouts'][layout], 'parallel_lanes': lanes}
    norm = sum(abs(w) for w, _, _ in branches)
    target = branches[0][2]['growth_target']
    nts = [w.fixed_T if arm == 'hybrid' else 0 for _, w, _ in branches]
    nt_any = max(nts) > 0
    ports = np.array(port_menu(L, sw, layout) if nt_any else [1], dtype=float)[:, None, None]
    port_tiles = port_tile_charge(sw)*(ports-1) if sw['ft_delivery'] else 0.*ports
    factories = np.array(local['factory_choices'] if nt_any else [0], dtype=float)[None, :, None]
    minimum = max(target, local['source_end_distance'] if nt_any and not sw['separate_distances'] else 3)
    d = odd_distances(minimum)[None, None, :]
    shape = np.broadcast(ports, factories, d).shape
    tiles = np.broadcast_to(data_block_tiles(L, local) + port_tiles + factories*src['tiles'], shape)
    cycle = local['cycle_seconds']
    transfer = np.maximum(d, local['source_end_distance']) if sw['separate_distances'] else d
    budget = log(epsilon*fractions['logical']/norm)
    seconds, feasible = [], np.ones(shape, bool)
    for (weight, w, part), nt in zip(branches, nts):
        analog = part['serial_analog_seconds']*d/target/lanes
        handoff = np.where(d > target, 2*d*cycle*w.rotations/lanes, 0.)
        supply = consume = 0.
        if nt:
            service = source_service(src, local, d, transfer)
            supply = np.ceil(nt/factories)*service*cycle
            consume = local['ft_clocks_per_T']*nt*d*cycle/ports
        s = np.broadcast_to(2*d*cycle + np.maximum(analog+handoff+consume, supply)
                            + (w.rotations/lanes + nt/ports)*local['feedforward_seconds']
                            + local['reset_seconds'] + local['measurement_seconds'], shape)
        log_fault = (np.log(tiles*np.ceil(s/cycle)*local['surface_A'])
                     + (d+1)/2*log(local['p']/local['surface_threshold']))
        feasible &= log_fault <= budget - part['chi']/2
        seconds.append(s)
    index = first_feasible(feasible)
    pi, fi = np.nonzero(index >= 0)
    if not pi.size:
        return None
    di = index[pi, fi]
    s_sel = [s[pi, fi, di] for s in seconds]
    d_sel = d[0, 0, di]
    q = 2*tiles[pi, fi, di]*d_sel**2 + factories[0, fi, 0]*src['qubits']
    ok = np.all([s <= local['week_seconds'] for s in s_sel], axis=0)
    if not ok.any():
        return None
    chis = [part['chi'] for _, _, part in branches]
    weights = [w for w, _, _ in branches]
    log_P = combined_log_P(chis, s_sel, weights)
    logq = np.where(ok, sampling_log_q(q, log_P, fractions['stat'], local, epsilon), np.inf)
    k = min(range(len(q)), key=lambda i: (logq[i], log(q[i])+log_P[i]))
    if not np.isfinite(logq[k]):
        return None
    seconds_best = [float(s[k]) for s in s_sel]
    row = final_sampling(float(q[k]), seconds_best, chis, weights, fractions['stat'], local, epsilon)
    first = branches[0][2]
    row.update(distance=int(d_sel[k]), injection_distance=first['seed'], growth_target=target,
               factories=int(factories[0, fi[k], 0]), cleanup_levels=src['cleanup'], T_count=nts[0],
               T_source=src['name'] if nt_any else None,
               data_tiles=int(data_block_tiles(L, local) + port_tiles[pi[k], 0, 0]),
               protected_tiles=float(tiles[pi[k], fi[k], di[k]]), layout=layout,
               lanes=lanes, T_consumption_ports=int(ports[pi[k], 0, 0]), rus_mode=first['rus_mode'],
               branch_seconds=seconds_best, branch_chi=chis,
               branch_T=nts, growth_probability=first['growth_probability'],
               coherent_bias_bound=[p['coherent_bias_bound'] for _, _, p in branches],
               rus_tail_bound=[p['rus_tail_bound'] for _, _, p in branches],
               serial_analog_seconds=[p['serial_analog_seconds'] for _, _, p in branches])
    return row


def optimize_analog(L, tau, arm, epsilon, p, scenario, sw):
    c = apply_delivery(base_config(p, scenario), sw)
    options = budget_settings()
    algorithm = 'termwise' if arm == 'star' else 'plaquette'
    seeds = [9, 25] if scenario == 'zero_growth' else v1.settings()['injection_seeds']
    pairs = [(s, t) for s in seeds for t in v1.settings()['growth_targets'] if s <= t]
    residual = options['residual_bias_fraction'][arm]
    candidates, visited = [], set()
    prune = settings()['chi_prune']

    def parts_for(w, seed, target):
        """Components by mode; PCEC is skipped where raw already exceeds the prune."""
        raw = component_exposure(w, p, seed, scenario, target, 'raw')
        out = {'raw': raw} if 'raw' in sw['rus_modes'] else {}
        if 'noisy_pcec' in sw['rus_modes'] and raw['chi'] <= prune:
            out['noisy_pcec'] = component_exposure(w, p, seed, scenario, target, 'noisy_pcec')
        return out

    def at_step(r, richardson):
        if (r, richardson) in visited:
            return
        visited.add((r, richardson))
        plan = branch_plan(r, richardson)
        works = [(weight, workload(L, tau, steps, algorithm, c['u_over_j'])) for steps, weight in plan]
        if arm == 'hybrid' and sw.get('fixed_T_analog'):
            works = [(weight, analog_fixed_T(w)) for weight, w in works]
        trotter = trotter_fraction(L, tau, algorithm, r, richardson, epsilon)
        norm = sum(abs(w) for w, _ in works)
        for seed, target in pairs:
            parts = [parts_for(w, seed, target) for _, w in works]
            for mode in sw['rus_modes']:
                if any(mode not in x for x in parts):
                    continue
                branches = [(weight, w, x[mode]) for (weight, w), x in zip(works, parts)]
                if max(b[2]['chi'] for b in branches) > prune:
                    continue
                instrument = sum(abs(weight)*(b['coherent_bias_bound']+b['rus_tail_bound'])
                                 for weight, _, b in branches)/epsilon
                digital = arm == 'hybrid' and any(w.fixed_T for _, w, _ in branches)
                for src in source_menu(c, sw, options['cleanup_choices'] if digital else [0]):
                    magic = 0.
                    if digital:
                        logs = [log(abs(weight)*w.fixed_T*src['error']/epsilon) + b['chi']/2 for weight, w, b in branches]
                        log_magic = float(np.logaddexp.reduce(logs))
                        if log_magic >= 0:
                            continue
                        magic = exp(log_magic)*(1+1e-12)
                    for logical in options['logical_fraction_grid']:
                        stat = 1-trotter-residual-instrument-magic-logical
                        if stat <= 0:
                            continue
                        fractions = dict(trotter=trotter, residual_bias=residual, instrument=instrument,
                                         magic=magic, logical=logical, stat=stat, synthesis=0.)
                        for lanes in lane_menu(L, algorithm, sw):
                            for layout in (analog_layouts(arm, lanes, sw) if arm == 'hybrid' or sw['optimize_lanes'] else ['compact']):
                                row = encoded_candidate(L, arm, c, epsilon, branches, fractions,
                                                        src, layout, lanes, sw)
                                if row is not None:
                                    row.update(L=L, tau=tau, arm=arm, absolute_error=epsilon, p=p,
                                               scenario=scenario, algorithm=algorithm, r=r,
                                               branches=[[s, w] for s, w in plan], extrapolation='richardson' if richardson else 'single',
                                               extrapolation_norm=norm, rotations=works[0][1].rotations,
                                               theta=works[0][1].theta, budget_fractions=fractions,
                                               trotter_coefficient=(richardson_coefficient(L, tau, algorithm) if richardson
                                                                    else trotter_coefficient(L, tau, algorithm)))
                                    candidates.append(row)
    for richardson in ([False, True] if sw['richardson'] else [False]):
        C = richardson_coefficient(L, tau, algorithm) if richardson else trotter_coefficient(L, tau, algorithm)
        menu = richardson_steps(C, tau, epsilon, options) if richardson else step_candidates(C, tau, epsilon, options)
        for r in menu:
            at_step(r, richardson)
        own = [v for v in candidates if (v['extrapolation'] == 'richardson') == richardson]
        if own:
            center = min(own, key=cost_key)['r']
            for r in range(max(ceil(options['minimum_steps_per_tau']*tau), center-2), center+3):
                at_step(r, richardson)
    if not candidates:
        return dict(L=L, tau=tau, arm=arm, absolute_error=epsilon, p=p, scenario=scenario,
                    log10_Q_min=None, status='No feasible candidate in the declared menu')
    best = min(candidates, key=cost_key)
    return dict(best, candidate_count=len(candidates), status='Conditional component composition; declared delivery model')


# ---------------------------------------------------------------- full FT

def optimize_ft(L, tau, epsilon, p, scenario, sw):
    c = apply_delivery(base_config(p, scenario), sw)
    options = budget_settings()
    cfg = settings()
    batches = [b for b in (cfg['extended_hwp_batches'] if sw['ft_delivery'] else cfg['reference_hwp_batches'])
               if b <= L*L]
    factories = np.array(c['factory_choices'], dtype=float)[None, :, None]
    d = odd_distances(3)[None, None, :]
    cycle = c['cycle_seconds']
    lanes_tiles = 0 if sw['ft_delivery'] else cfg['reference_lanes']
    best, count, visited, groups = None, 0, set(), {}

    def at_step(r, richardson):
        nonlocal best, count
        if (r, richardson) in visited:
            return
        visited.add((r, richardson))
        plan = branch_plan(r, richardson)
        norm = sum(abs(w) for _, w in plan)
        trotter = trotter_fraction(L, tau, 'plaquette', r, richardson, epsilon)
        if trotter >= 1:
            return
        works = [(weight, workload(L, tau, steps, 'plaquette', c['u_over_j'])) for steps, weight in plan]
        for layout in v1.settings()['layouts']:
            block = v1.settings()['layouts'][layout]
            ports = np.array(port_menu(L, sw, layout), dtype=float)[:, None, None]
            port_tiles = port_tile_charge(sw)*(ports-1) if sw['ft_delivery'] else 0.
            for batch in batches:
                style = 'individual' if batch == 1 else 'hamming_weight'
                for synthesis in options['synthesis_fraction_grid']:
                    local = {**c, **block, 'rotation_synthesis': style, 'hwp_batch_size': batch,
                             'synthesis_fraction': synthesis, 'parallel_lanes': lanes_tiles}
                    plans = [synthesis_plan(w, epsilon/norm, local) for _, w in works]
                    workspace = max(x.get('hwp_workspace', 0) for x in plans)
                    data_tiles = data_block_tiles(L, {**local, 'extra_work_qubits': workspace}) + port_tiles
                    for src in source_menu(c, sw, options['cleanup_choices']):
                        qout, cleanup = src['error'], src['cleanup']
                        lm = float(np.logaddexp.reduce([log(abs(weight)*x['T_count']*qout/epsilon)
                                                        for (weight, _), x in zip(works, plans)]))
                        if lm >= 0:
                            continue
                        magic = exp(lm)*(1+1e-12)
                        tiles = data_tiles + factories*src['tiles']
                        transfer = np.maximum(d, c['source_end_distance']) if sw['separate_distances'] else d
                        service = source_service(src, c, d, transfer)
                        clifford_ports = ports if sw['ft_delivery'] else cfg['reference_lanes']
                        durations, serial = [], []
                        for (_, w), x in zip(works, plans):
                            nt = x['T_count']
                            rotations = x.get('synthesized_rotations', w.rotations)
                            # Rotation sequences are serial per port; fixed and Toffoli T parallelize.
                            t_serial = (serial_T_layers(w, x, ports)*x['T_per_rotation']
                                        + (nt - rotations*x['T_per_rotation'])/ports)
                            supply = np.ceil(nt/factories)*service
                            consume = (block['ft_clocks_per_T']*t_serial + x['clifford_clock_units']/clifford_ports)*d
                            durations.append((2*d+np.maximum(supply, consume))*cycle + t_serial*c['feedforward_seconds']
                                             + c['reset_seconds'] + c['measurement_seconds'])
                            serial.append(t_serial)
                        shape = np.broadcast(ports, factories, d).shape
                        durations = [np.broadcast_to(s, shape) for s in durations]
                        tiles_b = np.broadcast_to(tiles, shape)
                        for logical in options['logical_fraction_grid']:
                            stat = 1-trotter-options['residual_bias_fraction']['ft']-synthesis-magic-logical
                            if stat <= 0:
                                continue
                            budget = log(epsilon*logical/norm)
                            feasible = np.ones(shape, bool)
                            for s in durations:
                                feasible &= (np.log(tiles_b*np.ceil(s/cycle)*c['surface_A'])
                                             + (d+1)/2*log(c['p']/c['surface_threshold'])) <= budget
                            index = first_feasible(feasible)
                            if not (index >= 0).any():
                                continue
                            pi, fi = np.nonzero(index >= 0)
                            first = d[0, 0, index[pi, fi]]
                            chosen = first if sw['separate_distances'] else np.maximum(first, c['source_end_distance'])
                            di = ((chosen-3)/2).astype(int)
                            s_sel = [s[pi, fi, di] for s in durations]
                            q = 2*tiles_b[pi, fi, di]*chosen**2 + factories[0, fi, 0]*src['qubits']
                            ok = np.all([s <= c['week_seconds'] for s in s_sel], axis=0)
                            if not ok.any():
                                continue
                            weights = [w for _, w in plan]
                            log_P = combined_log_P([0.]*len(plan), s_sel, weights)
                            logq = np.where(ok, sampling_log_q(q, log_P, stat, c, epsilon), np.inf)
                            k = min(range(len(q)), key=lambda i: (logq[i], log(q[i])+log_P[i]))
                            count += len(q)
                            if not np.isfinite(logq[k]):
                                continue
                            seconds = [float(s[k]) for s in s_sel]
                            row = final_sampling(float(q[k]), seconds, [0.]*len(plan), weights, stat, c, epsilon)
                            if row['log10_Q_min'] is None:
                                continue
                            nt0 = plans[0]['T_count']
                            supply0 = float(np.ceil(nt0/factories[0, fi[k], 0])*service[0, 0, di[k]]*cycle)
                            consume0 = float((block['ft_clocks_per_T']*np.broadcast_to(serial[0], ports.shape)[pi[k], 0, 0]
                                              + plans[0]['clifford_clock_units']/(ports[pi[k], 0, 0] if sw['ft_delivery'] else cfg['reference_lanes']))
                                             * chosen[k]*cycle)
                            row.update(L=L, tau=tau, arm='ft', absolute_error=epsilon, p=p, scenario=scenario,
                                       algorithm='plaquette', r=r, branches=[[s, w] for s, w in plan],
                                       extrapolation='richardson' if richardson else 'single', extrapolation_norm=norm,
                                       layout=layout, synthesis=style, hwp_batch_size=batch,
                                       hwp_workspace=workspace, T_count=nt0, branch_T=[x['T_count'] for x in plans],
                                       T_per_rotation=plans[0]['T_per_rotation'], cleanup_levels=cleanup, T_source=src['name'],
                                       magic_error=qout, factories=int(factories[0, fi[k], 0]),
                                       T_consumption_ports=int(ports[pi[k], 0, 0]), distance=int(chosen[k]),
                                       first_feasible_distance=int(first[k]), data_tiles=data_tiles if np.isscalar(data_tiles) else int(data_tiles[pi[k], 0, 0]),
                                       protected_tiles=float(tiles_b[pi[k], fi[k], di[k]]),
                                       branch_seconds=seconds, branch_chi=[0.]*len(plan),
                                       supply_seconds=supply0, consumption_seconds=consume0,
                                       delivery_bottleneck='supply' if supply0 > consume0 else 'consumption',
                                       budget_fractions=dict(stat=stat, trotter=trotter, logical=logical, synthesis=synthesis,
                                                             magic=magic, residual_bias=options['residual_bias_fraction']['ft'],
                                                             instrument=0.),
                                       trotter_coefficient=(richardson_coefficient(L, tau, 'plaquette') if richardson
                                                            else trotter_coefficient(L, tau, 'plaquette')))
                            if best is None or cost_key(row) < cost_key(best):
                                best = row
                            group = (layout, batch, richardson)
                            if group not in groups or cost_key(row) < cost_key(groups[group]):
                                groups[group] = row
    for richardson in ([False, True] if sw['richardson'] else [False]):
        C = richardson_coefficient(L, tau, 'plaquette') if richardson else trotter_coefficient(L, tau, 'plaquette')
        menu = richardson_steps(C, tau, epsilon, options) if richardson else step_candidates(C, tau, epsilon, options)
        for r in menu:
            at_step(r, richardson)
        # integrated_v1 refined the step count separately for every layout and batch.
        centers = sorted({v['r'] for g, v in groups.items() if g[2] == richardson})
        for center in centers:
            for r in range(max(ceil(options['minimum_steps_per_tau']*tau), center-2), center+3):
                at_step(r, richardson)
    if best is None:
        return dict(L=L, tau=tau, arm='ft', absolute_error=epsilon, p=p, scenario=scenario,
                    log10_Q_min=None, status='No feasible candidate in the declared menu')
    return dict(best, candidate_count=count, status='Synthesis counts checked; declared delivery model')


# ---------------------------------------------------------------- NISQ

def optimize_nisq(L, tau, epsilon, p, model, sw):
    c = {**sweep_configuration(p), 'nisq_model': {'nisq_effective': 'effective', 'nisq_pec': 'pec_two_qubit'}[model]}
    options = budget_settings()
    C = trotter_coefficient(L, tau, 'plaquette')
    single = optimize(L, tau, 'nisq', c, epsilon, options, C_override=C)
    single.update(arm=model, extrapolation='single', branches=[[single.get('r'), 1.]], p=p)
    rows = [single] if single.get('log10_Q_min') is not None else []
    if sw['richardson']:
        K4 = richardson_coefficient(L, tau, 'plaquette')
        bias = options['residual_bias_fraction']['nisq']
        visited = set()

        def at_step(r):
            if r in visited:
                return
            visited.add(r)
            trotter = K4/r**4/epsilon
            stat = 1-trotter-bias
            if stat <= 0:
                return
            plan = branch_plan(r, True)
            parts = [estimate(L, tau, 'nisq', {**c, 'stat_fraction': stat}, epsilon, steps=s) for s, _ in plan]
            rule = parts[0]['shot_rule']
            row = final_sampling(parts[0]['qubits_per_copy'], [v['seconds_per_attempt'] for v in parts],
                                 [v['chi'] for v in parts], [w for _, w in plan], stat, c, epsilon, rule)
            if row['log10_Q_min'] is None:
                return
            row.update(L=L, tau=tau, arm=model, absolute_error=epsilon, p=p, r=r, extrapolation='richardson',
                       branches=[[s, w] for s, w in plan], branch_chi=[v['chi'] for v in parts],
                       branch_seconds=[v['seconds_per_attempt'] for v in parts],
                       two_qubit_gates=[v['two_qubit_gates'] for v in parts], nisq_model=c['nisq_model'],
                       budget_fractions=dict(stat=stat, trotter=trotter, residual_bias=bias),
                       trotter_coefficient=K4, extrapolation_norm=5/3)
            rows.append(row)
        for r in richardson_steps(K4, tau, epsilon, options):
            at_step(r)
        own = [v for v in rows if v['extrapolation'] == 'richardson']
        if own:
            center = min(own, key=lambda v: v['log10_Q_min'])['r']
            for r in range(max(ceil(options['minimum_steps_per_tau']*tau), center-2), center+3):
                at_step(r)
    if not rows:
        return dict(L=L, tau=tau, arm=model, absolute_error=epsilon, p=p, log10_Q_min=None,
                    status='No feasible candidate')
    best = min(rows, key=lambda v: (v['log10_Q_min'], v['seconds_per_attempt']))
    return dict(best, scenario='component', status='NISQ effective overhead assumed; PEC is a specified reference')


# ---------------------------------------------------------------- driver

def work_items():
    cfg = settings()
    items = []
    for config in cfg['configurations']:
        scenarios = ['component'] + [s for s, names in cfg['scenario_configurations'].items() if config in names]
        for family, f in cfg['families'].items():
            for p in cfg['physical_errors']:
                for eps in cfg['precisions']:
                    for L in f['sizes']:
                        for scenario in scenarios:
                            arms = ENCODED + NISQ if scenario == 'component' else ('star', 'hybrid')
                            for arm in arms:
                                items.append((config, family, L, p, eps, scenario, arm))
    return items


def item_key(v):
    return (v['configuration'], v['family'], v['L'], v['p'], v['absolute_error'], v['scenario'], v['arm'])


def evaluate(item):
    config, family, L, p, eps, scenario, arm = item
    sw = switches(config)
    tau = family_tau(family, L)
    NEW_CHANNELS.clear()
    if arm == 'ft':
        row = optimize_ft(L, tau, eps, p, scenario, sw)
    elif arm in NISQ:
        row = optimize_nisq(L, tau, eps, p, arm, sw)
    else:
        row = optimize_analog(L, tau, arm, eps, p, scenario, sw)
    row.update(configuration=config, family=family, switches=sw, scenario=scenario)
    preparation.cache_clear()
    raw_rus_details.cache_clear()
    return row, dict(NEW_CHANNELS)


def _init(channels):
    CHANNELS.update(channels)


def phases(items):
    """Configurations sharing expensive PCEC angles run after their cache exists."""
    order = [['v1_menu'], ['raw_rus', 'lanes', 'ft_delivery', 'distances'], ['richardson'], ['all']]
    listed = {n for names in order for n in names}
    order.append(sorted({i[0] for i in items} - listed))   # sensitivity configurations of later rounds
    return [[i for i in items if i[0] in names] for names in order]


def run(workers):
    load_channels()
    contract = input_hashes()
    if (OUT/'inputs.json').exists() and read(OUT/'inputs.json') != contract:
        raise ValueError('Inputs changed; remove this version\'s costs.json and inputs.json to recompute')
    write('inputs.json', contract)
    rows = read(OUT/'costs.json') if (OUT/'costs.json').exists() else []
    done = {item_key(v) for v in rows}
    for phase in phases(work_items()):
        todo = [i for i in phase if i not in done]
        # Expensive items first keep the pool busy at the end.
        todo.sort(key=lambda i: (-i[2], i[6] != 'hybrid'))
        if not todo:
            continue
        with mp.get_context('fork').Pool(workers, initializer=_init, initargs=(dict(CHANNELS),)) as pool:
            for n, (row, channels) in enumerate(pool.imap_unordered(evaluate, todo, chunksize=1)):
                rows.append(row)
                CHANNELS.update(channels)
                if n % 25 == 0 or n == len(todo)-1:
                    write('costs.json', rows)
                    save_channels()
                    print('optimization round', len(rows), 'rows', flush=True)
    finish(rows)


# ---------------------------------------------------------------- analysis

def capacities(rows):
    cfg = settings()
    table = []
    for config in cfg['configurations']:
        scenarios = ['component'] + [s for s, names in cfg['scenario_configurations'].items() if config in names]
        for scenario in scenarios:
            for family in cfg['families']:
                for p in cfg['physical_errors']:
                    for eps in cfg['precisions']:
                        group = [v for v in rows if (v['configuration'], v['family'], v['p'], v['absolute_error']) == (config, family, p, eps)]
                        for Q in cfg['qubit_budgets']:
                            row = dict(configuration=config, scenario=scenario, family=family, p=p, epsilon=eps, Q=Q)
                            for arm in ENCODED + NISQ:
                                arm_scenario = scenario if arm in ('star', 'hybrid') else 'component'
                                row[arm] = capacity([v for v in group if v['arm'] == arm and v['scenario'] == arm_scenario], Q)
                            best = max(row[a] for a in ENCODED)
                            leaders = [a for a in ENCODED if row[a] == best]
                            row['encoded_leader'] = 'none' if best == 0 else ('tie:'+'+'.join(leaders) if len(leaders) > 1 else leaders[0])
                            row['nisq_effective_at_or_above_encoded'] = row['nisq_effective'] >= best and row['nisq_effective'] > 0
                            row['at_size_ceiling'] = max(row[a] for a in ENCODED + NISQ) == max(cfg['families'][family]['sizes'])
                            table.append(row)
    return table


def changes(table):
    """Cells where a switch changes the encoded leader relative to v1_menu."""
    reference = settings().get('baseline_configuration', 'v1_menu')
    base = {(v['family'], v['p'], v['epsilon'], v['Q']): v for v in table
            if v['configuration'] == reference and v['scenario'] == 'component'}
    out = []
    for v in table:
        if v['configuration'] == reference or v['scenario'] != 'component':
            continue
        b = base[(v['family'], v['p'], v['epsilon'], v['Q'])]
        moved = {a: [b[a], v[a]] for a in ENCODED + NISQ if b[a] != v[a]}
        if moved or b['encoded_leader'] != v['encoded_leader']:
            out.append(dict(configuration=v['configuration'], family=v['family'], p=v['p'], epsilon=v['epsilon'],
                            Q=v['Q'], leader_before=b['encoded_leader'], leader_after=v['encoded_leader'],
                            leader_changed=b['encoded_leader'] != v['encoded_leader'], capacity_moves=moved))
    return out


def input_hashes():
    names = ['python/optimization_round.py', str(CONFIG.relative_to(ROOT)), 'python/integrated_comparison.py',
             'configs/integrated_comparison.json', 'python/star_channel_study.py', 'python/physics_sweep.py',
             'python/discussion_model.py', 'python/budget_optimization.py', 'python/resource_extensions.py',
             'python/model.py', 'python/dynamics_followup_analysis.py', 'configs/assumptions.json',
             'configs/discussion_model.json', 'configs/budget_optimization.json', 'configs/patch_sensitivity.json',
             'results/discussion_v3/calibration.json', 'results/discussion_v3/finite_step_guard.json',
             'results/discussion_v3/rectangle_checks/calibration.json',
             'results/dynamics_followup_v1/measurements.json', 'results/frontier_v1/growth_summary.json']
    if settings().get('surface_law'):
        names.append(settings()['surface_law']['source'])
    if settings().get('calibrated_factories'):
        names.append(settings()['calibrated_factories'])
    return {name: sha(ROOT/name) for name in names}


def validate(rows):
    """Replay integrated_v1 and physics_sweep_v1 (v1) or round-v1 NISQ rows (v2), and audit every ledger."""
    old = read(v1.OUT/'costs.json')
    replays = 0
    calibrated = bool(settings().get('surface_law'))
    if calibrated:
        # The surface-code law cannot affect physical-qubit rows: they must equal round v1.
        base = {item_key(v): v for v in read(ROOT/ROUNDS['v1'][1]/'costs.json')}
        for v in rows:
            if v['arm'] in NISQ and item_key(v) in base:
                a, b = base[item_key(v)]['log10_Q_min'], v['log10_Q_min']
                assert (a is None) == (b is None) and (a is None or isclose(a, b, abs_tol=1e-12)), item_key(v)
                replays += 1
    replay = settings().get('replay_round')
    if replay:
        # The reference configuration of a sensitivity round must equal the round it extends.
        prior = {item_key(v): v for v in read(ROOT/ROUNDS[replay][1]/'costs.json')}
        for v in rows:
            if v['configuration'] == 'all':
                a, b = prior[item_key(v)]['log10_Q_min'], v['log10_Q_min']
                assert (a is None) == (b is None) and (a is None or isclose(a, b, abs_tol=1e-12)), item_key(v)
                replays += 1
    for v in rows:
        if calibrated or v['configuration'] != 'v1_menu' or v['family'] != 'tau_1':
            continue
        if v['arm'] in ENCODED:
            match = [o for o in old if (o['L'], o['tau'], o['absolute_error'], o['p'], o['arm'], o['scenario'])
                     == (v['L'], 1, v['absolute_error'], v['p'], v['arm'], v['scenario'])]
        else:
            sweep = read(ROOT/'results/physics_sweep_v1/costs.json')
            model = {'nisq_effective': 'effective', 'nisq_pec': 'pec_two_qubit'}[v['arm']]
            match = [o for o in sweep if (o['L'], o['absolute_error'], o['p'], o['arm'], o['scenario'])
                     == (v['L'], v['absolute_error'], v['p'], 'nisq', model)]
        if not match:
            continue  # sizes added after integrated_v1 have nothing to replay
        assert len(match) == 1, (v['arm'], v['L'])
        a, b = match[0]['log10_Q_min'], v['log10_Q_min']
        if a is not None and a > 10 and b is None:
            # The exposure prune may drop a design v1 costed above 10^10 qubits.
            continue
        assert (a is None) == (b is None) and (a is None or isclose(a, b, abs_tol=1e-9)), (v['arm'], v['L'], a, b)
        replays += 1
    for tau in (1, 2, 4):
        for algorithm in ('termwise', 'plaquette'):
            assert trotter_coefficient(4, tau, algorithm) == extended_coefficient(tau, algorithm, _measurements())
    for v in rows:
        if v['log10_Q_min'] is None:
            continue
        b = v['budget_fractions']
        assert isclose(sum(b.values()), 1, abs_tol=1e-12) and min(b.values()) >= 0, item_key(v)
        if v['arm'] in ('star', 'hybrid', 'ft'):
            c = base_config(v['p'], v['scenario'])
            pL = (v['distance']+1)/2*log(c['p']/c['surface_threshold']) + log(c['surface_A'])
            for s, chi in zip(v['branch_seconds'], v['branch_chi']):
                fault = log(v['protected_tiles']*ceil(s/c['cycle_seconds'])) + pL
                assert fault + chi/2 <= log(v['absolute_error']*b['logical']/v['extrapolation_norm']) + 1e-9, item_key(v)
            sw = v['switches']
            if v['arm'] in ('star', 'hybrid'):
                assert v['distance'] >= v['growth_target']
                if v['T_count'] and not sw['separate_distances']:
                    assert v['distance'] >= c['source_end_distance']
            if v['arm'] == 'ft' and not sw['separate_distances']:
                assert v['distance'] >= c['source_end_distance']
            if not sw['optimize_lanes'] and v['arm'] != 'ft':
                assert v['lanes'] == settings()['reference_lanes']
            if not sw['ft_delivery']:
                assert v['T_consumption_ports'] == 1
    pairs = {}
    for key, value in CHANNELS.items():
        stem, mode = key.rsplit('|', 1)
        pairs.setdefault(stem, {})[mode] = value['q_dephasing']
    compared = [x for x in pairs.values() if len(x) == 2]
    # Justifies skipping PCEC once raw RUS already exceeds the exposure prune.
    assert all(x['noisy_pcec'] >= x['raw']*(1-1e-9) for x in compared)
    return dict(status=('PASS: round-v1 NISQ replay, calibrated surface law, budgets, fault ledgers, switch domains, PCEC >= raw'
                        if calibrated else
                        'PASS: integrated_v1/physics_sweep_v1 replay, budgets, fault ledgers, switch domains, PCEC >= raw'),
                resource_rows=len(rows), v1_replays=replays, channel_pairs_compared=len(compared),
                channel_points=len(CHANNELS), scope=settings()['scope'])


def finish(rows):
    table = capacities(rows)
    write('capacity.json', table)
    write_csv('capacity.csv', table)
    delta = changes(table)
    write('changes.json', delta)
    write('validation.json', validate(rows))
    save_channels()
    print('Completed', len(rows), 'rows;', sum(v['leader_changed'] for v in delta), 'leader changes', flush=True)


def check():
    load_channels()
    assert read(OUT/'inputs.json') == input_hashes(), 'Optimization-round inputs changed'
    rows = read(OUT/'costs.json')
    assert {item_key(v) for v in rows} == set(work_items()), 'Incomplete or extra work items'
    assert read(OUT/'capacity.json') == capacities(rows)
    assert read(OUT/'changes.json') == changes(capacities(rows))
    assert read(OUT/'validation.json') == validate(rows)
    print('Optimization round verified')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--round', choices=sorted(ROUNDS), default='v1')
    parser.add_argument('--finish', action='store_true', help='Recompute tables from saved rows')
    args = parser.parse_args()
    use_round(args.round)
    if args.check:
        check()
    elif args.finish:
        load_channels()
        finish(read(OUT/'costs.json'))
    else:
        run(args.workers)
