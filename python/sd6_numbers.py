"""Numbers for the SD6 manuscript pass (rounds v4/v5), printed from saved results only."""
import json
import numpy as np
from project_paths import ROOT

Qs = np.geomspace(1e1, 1e10, 1801)


def load(r):
    return json.loads((ROOT/f'results/optimization_round_{r}/costs.json').read_text())


def cost(rows, cfg, arm, p, eps, family='tau_L_over_4', scen='component'):
    return {v['L']: v['log10_Q_min'] for v in rows if v['configuration'] == cfg and v['arm'] == arm and v['p'] == p
            and v['absolute_error'] == eps and v['scenario'] == scen and v['family'] == family and v['log10_Q_min'] is not None}


def cap(c):
    return np.array([max([L for L, q in c.items() if q <= np.log10(Q)], default=0) for Q in Qs])


def crossings(rows, cfg, p, eps, family='tau_L_over_4'):
    ft, st, hy = (cap(cost(rows, cfg, a, p, eps, family)) for a in ('ft', 'star', 'hybrid'))
    enc = np.maximum(ft, np.maximum(st, hy)); an = np.maximum(st, hy)
    out = {}
    for n in ('nisq_pec', 'nisq_effective'):
        ni = cap(cost(rows, cfg, n, p, eps, family)); b = np.nonzero(enc <= ni)[0]; out[n] = Qs[b[-1]+1]
    b = np.nonzero(ft < an)[0]; out['never_behind'] = Qs[b[-1]+1] if len(b) else Qs[0]
    s = np.nonzero(ft <= an)[0]; out['strict'] = Qs[s[-1]+1]
    m = Qs <= 1e9
    out['star_sat'] = int(st[m].max()); out['hyb_sat'] = int(hy[m].max())
    f = cost(rows, cfg, 'ft', p, eps, family); out['ft_L32'] = 10**f[32] if 32 in f else None
    out['physical_eff'] = cost(rows, cfg, 'nisq_effective', p, eps, family)
    return out


if __name__ == '__main__':
    v5 = load('v5')
    for cfg in ('calibrated_factories', 'all', 'litinski_sources', 'cultivation_exact', 'pcec_only',
                'lanes_N', 'lane_tiles_3', 'lane_tiles_12', 'port_tiles_8', 'star_optimistic', 'star_plaquette', 'better_source'):
        print('==', cfg)
        for p in (1e-3, 1e-4):
            for eps in (1e-2, 1e-3):
                x = crossings(v5, cfg, p, eps)
                print(f"  p={p:g} eps={eps:g}: enc>phys PEC {x['nisq_pec']:.2g} eff {x['nisq_effective']:.2g} | FT never behind {x['never_behind']:.2g} "
                      f"strict {x['strict']:.2g} | STAR L{x['star_sat']} hyb L{x['hyb_sat']} | FT L32 {x['ft_L32'] and '%.2g'%x['ft_L32']}")
    print('physical effective costs p=1e-3 eps=.01', {L: round(q, 1) for L, q in cost(v5, 'all', 'nisq_effective', 1e-3, .01).items()})
    print('physical PEC costs p=1e-3 eps=.01', {L: round(q, 1) for L, q in cost(v5, 'all', 'nisq_pec', 1e-3, .01).items()})
    print('star_plaquette L12', round(cost(v5, 'star_plaquette', 'hybrid', 1e-3, .01)[12], 1), 'STAR', round(cost(v5, 'all', 'star', 1e-3, .01)[12], 1),
          'hybrid', round(cost(v5, 'all', 'hybrid', 1e-3, .01)[12], 1))
    ch = json.loads((ROOT/'results/optimization_round_v5/changes.json').read_text())
    from collections import Counter
    print('v5 leader changes', Counter(c['configuration'] for c in ch if c['leader_changed']))
    t = json.loads((ROOT/'results/optimization_round_v5/capacity.json').read_text())
    for cfg in ('calibrated_factories', 'all'):
        print(cfg, 'non-FT leaders:', [(v['p'], v['epsilon'], '%.0e' % v['Q'], v['encoded_leader']) for v in t
               if v['configuration'] == cfg and v['scenario'] == 'component' and v['encoded_leader'] not in ('ft', 'none')])
    for cfg in ('calibrated_factories',):
        for p in (1e-3, 1e-4):
            for eps in (1e-2, 1e-3):
                ft, st, hy = (cap(cost(v5, cfg, a, p, eps)) for a in ('ft', 'star', 'hybrid'))
                solo = (st > np.maximum(ft, hy)) & (Qs <= 1e9)
                if solo.any():
                    i = np.nonzero(solo)[0]; print('STAR alone', p, eps, '%.2g-%.2g' % (Qs[i[0]], Qs[i[-1]]))
