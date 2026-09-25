"""QUOPS capability for physical qubits, STAR, and full fault tolerance on one array.

QUOPS circuits (Proctor et al., arXiv:2609.12146) of width w and size s alternate
a CNOT layer on floor(w/2) uniformly random pairs with a layer of R_P(theta)
rotations, P in {X,Y,Z} and theta uniform on [0, 2 pi). Size counts one per
rotation and two per CNOT. A shape passes when the mean process polarization is
at least 1/sqrt(e) (standard) or, in the mitigated mode, when the 1/gamma^2
sampling overhead of the mirror test fits in one week.

For every shape and arm we compute the fewest physical qubits Q_min(w, s); the
capability at a budget Q is the largest s with Q_min <= Q inside w^2 <= s <= w^3.
The arms reuse the Hubbard study's instruments: the STAR generator-fault channel
at the circuit's random angles, the calibrated surface-code law, and FT
cultivation, synthesis, ports and factories. Routing and lattice-surgery timing
are declared models (configs/quops_study.json).
"""
import argparse
import csv
import hashlib
import json
import multiprocessing as mp
from functools import lru_cache
from math import ceil, exp, floor, log, log2, pi, sqrt

import numpy as np

from discussion_model import assumptions
from project_paths import ROOT
from star_channel_study import growth_error, preparation, raw_rus_details, rus_channel, wrap

import os
CONFIG = ROOT/os.environ.get('QUOPS_CONFIG', 'configs/quops_study.json')
OUT = ROOT/json.loads(CONFIG.read_text()).get('output', 'results/quops_v1')
ARMS = ('nisq', 'star', 'ft', 'ft_better_source')
MODES = ('standard', 'mitigated')


def read(path):
    return json.loads(path.read_text())


@lru_cache(maxsize=1)
def settings():
    return read(CONFIG)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT/name
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temp.replace(path)


@lru_cache(maxsize=None)
def hardware(p):
    c = assumptions()
    law = settings()['surface_law']
    fit = read(ROOT/law['source'])['fits'][str(p)]
    return {**c, 'p': p, 'surface_A': fit[law['amplitude']], 'surface_threshold': fit['effective_p_threshold']}


# ---------------------------------------------------------------- circuit shape

def shape(w, s):
    """Gate counts of a size-s QUOPS circuit of width w (fractional final layer allowed)."""
    pairs = w//2
    per_layer = w + 2*pairs
    layers = s/per_layer
    return dict(w=w, s=s, layers=layers, rotations=layers*w, cnots=layers*pairs, pairs=pairs)


def embed(w):
    """a_w = 4^w/(4^w-1): a gate's process infidelity embedded in the w-qubit polarization."""
    return 1/(1-4.**-w) if w < 60 else 1.


def in_cone(w, s):
    return w*w <= s <= w**3


def budget(mode, copies, shots, seconds, B):
    """Allowed -log(gamma): 1/2 for the standard test, or the week-limited mitigated overhead."""
    if mode == 'standard':
        return 0.5 if shots*seconds <= copies*B else -1.
    # Mitigated: shots/gamma^2 attempts must fit, and gamma must not exceed the standard threshold's meaning.
    return 0.5*log(copies*B/(shots*seconds)) if copies*B > shots*seconds else -1.


def test_shots():
    t = settings()['test']
    return t['circuits']*t['shots_per_circuit']*t['mirror_factor']


# ---------------------------------------------------------------- physical qubits

@lru_cache(maxsize=1)
def routing_data():
    """Measured per-layer means (area 1) and fits from results/quops_routing_v1."""
    import csv
    src = ROOT/settings()['routing']['source']
    rows = list(csv.DictReader((src/'raw.csv').open()))
    means = {}
    for router in ('cqc', 'constructive'):
        for w in sorted({int(r['width']) for r in rows}):
            sel = [r for r in rows if r['router'] == router and int(r['width']) == w and r['area_multiplier'] == '1']
            means[router, w] = (np.mean([float(r['fused_cnots_per_layer']) for r in sel]),
                                np.mean([float(r['fused_cnot_depth_per_layer']) for r in sel]))
    fits = read(src/'fit.json')['fits']
    return means, fits


def routed_layer(w, router):
    """(two-qubit gates, two-qubit depth) per benchmark layer on the grid.

    Measured means are log-interpolated for 4 <= w <= 144; the fitted laws
    apply beyond; two and three qubits need no routing.
    """
    if w <= 3 or settings().get('connectivity') == 'all_to_all':
        return float(w//2), 1.
    means, fits = routing_data()
    widths = sorted(v for r, v in means if r == router)
    if w <= widths[-1]:
        lo = max(v for v in widths if v <= w); hi = min(v for v in widths if v >= w)
        (g0, d0), (g1, d1) = means[router, lo], means[router, hi]
        t = 0. if hi == lo else (log(w)-log(lo))/(log(hi)-log(lo))
        return exp((1-t)*log(g0)+t*log(g1)), exp((1-t)*log(d0)+t*log(d1))
    f = fits[router+'_area1']
    g = f['fused_cnots_per_layer']; d = f['fused_cnot_depth_per_layer']
    return g['a']*w**g['b'], d['a']*w**d['b']


def nisq_cost(w, s, p, mode, Q):
    """Q_min for the physical-qubit arm, or None. Replicas only help the mitigated mode."""
    c = hardware(p)
    x = shape(w, s)
    side = ceil(sqrt(w))
    best = None
    for router in ('cqc', 'constructive'):
        g2, d2 = routed_layer(w, router)
        two_q = x['layers']*g2
        depth = x['layers']*(d2 + 1)                  # plus the rotation layer
        slots = w*depth - (2*two_q + x['rotations'])
        # The study's Stim convention (PLAN.md C6; results/surface_law_v1): DEPOLARIZE2(p) per
        # two-qubit gate, DEPOLARIZE1(p/10) per one-qubit gate, DEPOLARIZE1(p) per idle 10 ns slot.
        # A Pauli channel with total error probability q has process infidelity q.
        idle = max(0., slots)
        loss = embed(w)*(two_q*p + x['rotations']*p*settings().get('one_qubit_error_ratio', .1) + idle*p)
        if best is None or loss < best[0]:
            best = (loss, two_q, depth, router)
    loss, two_q, depth, router = best
    seconds = depth*c['gate_seconds'] + c['measurement_seconds'] + c['reset_seconds']
    q = side*side
    shots = test_shots()
    if mode == 'standard':
        copies = max(1, ceil(shots*seconds/c['week_seconds']))
        return dict(qubits=q*copies, loss=loss, seconds=seconds, copies=copies,
                    two_qubit_gates=two_q, depth=depth, router=router) if loss <= 0.5 else None
    if 2*loss > 69:
        return None      # beyond 10^30 attempts
    need = shots*exp(2*loss)*seconds/c['week_seconds']
    copies = max(1, ceil(need))
    if copies is None:
        return None
    return dict(qubits=q*copies, loss=loss, seconds=seconds, copies=copies, two_qubit_gates=two_q, depth=depth, router=router)


# ---------------------------------------------------------------- encoded arms

def surgery_clocks(w):
    """Code clocks per layer of long-range CNOTs; optional congestion growth with width."""
    g = settings()['surgery']
    base = g['clocks_per_cnot_layer']
    if 'congestion_width' in g:
        base *= max(1., (w/g['congestion_width'])**g['congestion_exponent'])
    return base


def tiles_block(n):
    """Fast Litinski block: every logical tile has an exposed boundary; one tile each with all-to-all links."""
    if settings().get('connectivity') == 'all_to_all':
        return n
    return ceil(2*n + sqrt(8*n) + 1)


def p_logical(d, c):
    return c['surface_A']*(c['p']/c['surface_threshold'])**((d+1)/2)


@lru_cache(maxsize=None)
def star_channel(p, seed, target, mode, g, angle_max=pi/8):
    """Mean over the injected-angle distribution of the accepted STAR rotation channel.

    For QUOPS, theta ~ U[0, 2 pi) and Clifford parts are absorbed, so the injected
    angle is uniform on [-pi/8, pi/8]. The angle sweep narrows this to
    [-angle_max, angle_max]. By symmetry we average |angle|.
    """
    k = settings()['star']['angle_samples']
    angles = [(i+.5)*angle_max/k for i in range(k)]
    scale = 0. if seed == target else growth_error(seed, p, 1e-10, target)/1e-10
    anchor = g/scale if scale else 0.
    rows = [rus_channel(a, seed, p, mode=mode, floor=anchor, data_distance=target, cap=48) for a in angles]
    preparation.cache_clear(); raw_rus_details.cache_clear()
    q = float(np.mean([v['q_dephasing'] for v in rows]))
    return dict(q=q, loss=q, coherent=float(np.mean([abs(v['coherent_component']) for v in rows])),
                seconds=float(np.mean([v['mean_seconds'] for v in rows])))


def growth_upper(seed, target, p):
    if seed == target:
        return 0.
    pairs = read(ROOT/settings().get('growth_summary', 'results/frontier_v1/growth_summary.json'))['pairs']
    return next(v['pauli_total_upper'] for v in pairs if (v['di'], v['df'], v['p']) == (seed, target, p))


def encoded_search(w, s, p, mode, tiles_fn, layer_fn, loss_fn, d_min):
    """Smallest area over distance for an encoded arm; tiles/time/loss given as functions of d."""
    c = hardware(p)
    x = shape(w, s)
    shots = test_shots()
    best = None
    for d in range(d_min if d_min % 2 else d_min+1, 101, 2):
        tiles, extra = tiles_fn(d)
        seconds = x['layers']*layer_fn(d) + 2*d*c['cycle_seconds'] + c['measurement_seconds'] + c['reset_seconds']
        cycles = seconds/c['cycle_seconds']
        loss = loss_fn(d) + tiles*cycles*p_logical(d, c)
        q = 2*tiles*d*d + extra
        if mode == 'standard':
            if loss > 0.5:
                continue
            copies = max(1, ceil(shots*seconds/c['week_seconds']))
        else:
            need = shots*exp(min(2*loss, 600.))*seconds/c['week_seconds']
            if need > 1e30:
                continue
            copies = max(1, ceil(need))
        row = dict(qubits=q*copies, loss=loss, seconds=seconds, copies=copies, distance=d, tiles=tiles)
        if best is None or row['qubits'] < best['qubits']:
            best = row
        elif mode == 'standard' and best is not None:
            break        # area grows with d once the loss budget is met
    return best


def star_cost(w, s, p, mode, Q):
    c, cfg = hardware(p), settings()['star']
    x = shape(w, s)
    best = None
    for seed, target in cfg['seed_targets']:
        g = growth_upper(seed, target, p)
        for rus in cfg['rus_modes']:
            ch = star_channel(p, seed, target, rus, g, cfg.get('angle_max', pi/8))
            rot_loss = embed(w)*x['rotations']*ch['loss']
            if rot_loss > (0.5 if mode == 'standard' else 60):
                continue
            for lanes in [n for n in cfg['lanes'] if n <= w]:
                def tiles_fn(d, lanes=lanes):
                    return tiles_block(w) + cfg['tiles_per_lane']*lanes, 0
                def layer_fn(d, lanes=lanes, ch=ch):
                    surgery = surgery_clocks(w)*d*c['cycle_seconds']
                    rotate = ch['seconds']*d/target*ceil(w/lanes)
                    return surgery + rotate + c['feedforward_seconds']*ceil(w/lanes)
                row = encoded_search(w, s, p, mode, tiles_fn, layer_fn, lambda d: rot_loss, target)
                if row and (best is None or row['qubits'] < best['qubits']):
                    best = dict(row, seed=seed, target=target, rus_mode=rus, lanes=lanes, rotation_loss=rot_loss)
    return best


def t_sources(c, q0):
    """Cultivation with optional 15-to-1 cleanup; optionally the Litinski (2019) factories at matching p."""
    menu = [dict(name='cultivation', cleanup=k, error=35*q0**3 if k else q0, tiles=1+12*k,
                 qubits=c['source_qubits']) for k in (0, 1)]
    cfg = settings()['ft']
    if cfg.get('source_menu') == 'litinski':
        table = read(ROOT/cfg['litinski_factories'])['litinski_factories']
        for name, error, qubits, cycles, outputs in table.get(repr(c['p']), []):
            menu.append(dict(name=name, cleanup=0, error=error, tiles=1, qubits=qubits, cycles_per_T=cycles/outputs))
    return menu


def ft_cost(w, s, p, mode, Q, source=None):
    c, cfg = hardware(p), settings()['ft']
    x = shape(w, s)
    if source is None:
        source = cfg.get('cultivation_error', c['source_error'])
    q0 = source
    best = None
    caps = [0.5] if mode == 'standard' else cfg['mitigated_loss_caps']
    for share, allowed in [(f, a) for f in cfg['synthesis_shares'] for a in caps]:
        # The mitigated search includes every standard design, so it cannot do worse.
        zeta = share*allowed/x['rotations']
        # Paper convention: diamond-distance synthesis error eps costs about eps^2 in polarization.
        n_T = ceil(3.02*log2(1/sqrt(zeta)) + .77)
        N_T = x['rotations']*n_T
        for src in t_sources(c, q0):
            cleanup = src['cleanup']
            magic = N_T*src['error']
            if magic + share*allowed > allowed:
                continue
            for ports in [n for n in cfg['ports'] if n <= w]:
                for factories in cfg['factories']:
                    def tiles_fn(d, ports=ports, factories=factories, src=src):
                        return (tiles_block(w) + cfg['tiles_per_port']*(ports-1)
                                + factories*src['tiles'], factories*src['qubits'])
                    def layer_fn(d, ports=ports, factories=factories, src=src):
                        surgery = surgery_clocks(w)*d
                        consume = ceil(w/ports)*n_T*d
                        source = src.get('cycles_per_T', c['source_cycles_per_T']) + 2*max(d, c['source_end_distance'])
                        service = (15*source+11*d)/(1-15*q0) if src['cleanup'] else source
                        supply = ceil(w*n_T/factories)*service
                        return (surgery + max(consume, supply))*c['cycle_seconds'] + ceil(w/ports)*n_T*c['feedforward_seconds']
                    row = encoded_search(w, s, p, mode, tiles_fn, layer_fn,
                                         lambda d, z=zeta*x['rotations'], m=magic: z + m, 3)
                    if row and (best is None or row['qubits'] < best['qubits']):
                        best = dict(row, synthesis_share=share, T_per_rotation=n_T, T_count=N_T,
                                    cleanup_levels=cleanup, T_source=src['name'], ports=ports, factories=factories)
    return best


def ft_better_source(w, s, p, mode, Q):
    """Sensitivity: a declared better cultivation source, same footprint and timing."""
    return ft_cost(w, s, p, mode, Q, source=settings()['ft']['better_source_error'])


def pinnacle_cost(w, s, p, mode, Q):
    """Pinnacle qLDPC architecture (python/pinnacle.py): one processing unit of ceil(w/k) blocks plus a
    magic engine per sampler; Cliffords are absorbed by Pauli-based computation, so the circuit takes
    one logical cycle (or one magic-engine attempt, if slower) per T state. Bounded-range connectivity only."""
    if settings().get('connectivity') != 'all_to_all':
        return None
    import pinnacle as pn
    c, cfg = hardware(p), settings()['ft']
    x = shape(w, s)
    shots = test_shots()
    best = None
    caps = [0.5] if mode == 'standard' else cfg['mitigated_loss_caps']
    for share, allowed in [(f, a) for f in cfg['synthesis_shares'] for a in caps]:
        zeta = share*allowed/x['rotations']
        n_T = ceil(3.02*log2(1/sqrt(zeta)) + .77)
        N_T = x['rotations']*n_T
        for engine in [e for e in pn.ENGINES if abs(e.p-p) < 1e-15]:
            magic = N_T*engine.output_error
            if magic + share*allowed > allowed:
                continue
            for code in pn.CODES:
                rounds = max(code.cycle_rounds, engine.cycles_per_attempt/engine.acceptance)
                seconds = N_T*rounds*c['cycle_seconds'] + c['measurement_seconds'] + c['reset_seconds']
                loss = share*allowed + magic + w*(N_T*rounds/code.cycle_rounds)*pn.logical_error(code, p)
                q = pn.unit_qubits(code, w, engine)
                if mode == 'standard':
                    if loss > 0.5:
                        continue
                    copies = max(1, ceil(shots*seconds/c['week_seconds']))
                else:
                    need = shots*exp(min(2*loss, 600.))*seconds/c['week_seconds']
                    if need > 1e30:
                        continue
                    copies = max(1, ceil(need))
                row = dict(qubits=q*copies, loss=loss, seconds=seconds, copies=copies, code=code.name,
                           engine=engine.name, T_per_rotation=n_T, T_count=N_T, synthesis_share=share)
                if best is None or row['qubits'] < best['qubits']:
                    best = row
    return best


COST = dict(nisq=nisq_cost, star=star_cost, ft=ft_cost, ft_better_source=ft_better_source, pinnacle=pinnacle_cost)


def arms():
    """Arms of the configured study; Pinnacle needs bounded-range (here all-to-all) connectivity."""
    return ARMS + (('pinnacle',) if settings().get('connectivity') == 'all_to_all' else ())


# ---------------------------------------------------------------- sweep

def widths():
    g = settings()['grid']
    small = list(range(2, g['integer_widths_to']+1))
    large = np.unique(np.round(np.geomspace(g['integer_widths_to']+1, g['max_width'], g['large_width_points'])).astype(int))
    return small + [int(v) for v in large if v > g['integer_widths_to']]


def sizes(w):
    g = settings()['grid']
    lo, hi = w*w, w**3
    return sorted({int(round(v)) for v in np.geomspace(lo, hi, g['sizes_per_width'])})


def sweep_arm(arm, p, mode):
    rows = []
    for w in widths():
        failed = 0
        for s in sizes(w):
            row = COST[arm](w, s, p, mode, None)
            if row is None:
                failed += 1
                if failed >= 2:
                    break      # cost is monotone in s at fixed w
                continue
            rows.append(dict(arm=arm, p=p, mode=mode, w=w, s=s, log10_Q_min=log(row['qubits'])/log(10),
                             **{k: v for k, v in row.items() if k != 'qubits'}, qubits=row['qubits']))
        if not any(v['w'] == w for v in rows) and w > 64:
            break              # wider circuits cannot pass once no size at this width does
    return rows


def capability(rows, Q):
    fit = [v['s'] for v in rows if v['qubits'] <= Q]
    return max(fit, default=0)


def run():
    cfg = settings()
    tasks = [(arm, p, mode) for p in cfg['physical_errors'] for mode in MODES for arm in arms()]
    with mp.get_context('fork').Pool(len(tasks)) as pool:
        parts = pool.starmap(sweep_arm, tasks)
    rows = [row for part in parts for row in part]
    for task, part in zip(tasks, parts):
        print('quops', *task, len(part), 'shapes', flush=True)
    write('costs.json', rows)
    finish(rows)


def finish(rows):
    cfg = settings()
    table = []
    for p in cfg['physical_errors']:
        for mode in MODES:
            for Q in np.geomspace(settings().get('q_min', 1e2), 1e9, int(round(10*np.log10(1e9/settings().get('q_min', 1e2))))+1):
                row = dict(p=p, mode=mode, Q=float(Q))
                for arm in arms():
                    # Mitigation is optional: the mitigated capability includes every standard design.
                    modes = ('standard',) if mode == 'standard' else MODES
                    row[arm] = capability([v for v in rows if v['arm'] == arm and v['p'] == p and v['mode'] in modes], Q)
                table.append(row)
    write('capability.json', table)
    write('inputs.json', {n: sha(ROOT/n) for n in ('python/quops_study.py', 'configs/quops_study.json',
                                                    'python/star_channel_study.py', 'configs/assumptions.json',
                                                    'configs/discussion_model.json', settings()['surface_law']['source'])})
    print('done', len(rows))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--finish', action='store_true', help='Recompute tables from saved costs')
    if parser.parse_args().finish:
        finish(read(OUT/'costs.json'))
    else:
        run()
