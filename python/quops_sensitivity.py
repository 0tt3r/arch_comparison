"""QUOPS angle sweep and delivery sensitivity, standard polarization threshold.

Angle sweep: STAR's injected angles are drawn uniformly on [-a, a] for a from
0.003 to pi/8 (the QUOPS value); physical-qubit and FT costs do not depend on
the angles, so this isolates STAR's per-radian exposure under one criterion.
Delivery: the lattice-surgery CNOT-layer time, a congestion penalty growing as
sqrt(w), and the lane and port tile charges, one at a time. Source: the
published Litinski (2019) factories offered alongside cultivation.
"""
import copy
import json
import multiprocessing as mp
from math import log, pi

import numpy as np

import quops_study as q
from project_paths import ROOT

OUT = ROOT/'results/quops_sensitivity_v1'
BASE = copy.deepcopy(q.settings())
ANGLES = [0.003, 0.01, 0.03, 0.1, pi/8]
DELIVERY = {
    'baseline': {},
    'surgery_x2': {'surgery': {'clocks_per_cnot_layer': 4}},
    'surgery_congestion': {'surgery': {'clocks_per_cnot_layer': 2, 'congestion_width': 16, 'congestion_exponent': .5}},
    'lane_tiles_3': {'star': {'tiles_per_lane': 3}},
    'lane_tiles_12': {'star': {'tiles_per_lane': 12}},
    'port_tiles_8': {'ft': {'tiles_per_port': 8}},
}
# T source: the published distillation factories of Litinski (2019) alongside cultivation.
SOURCES = {'baseline': {}, 'litinski_sources': {'ft': {'source_menu': 'litinski'}}}


def merged(changes):
    cfg = copy.deepcopy(BASE)
    for k, v in changes.items():
        cfg[k] = {**cfg[k], **v}
    return cfg


def task(args):
    name, changes, arm, p = args
    cfg = merged(changes)
    q.settings = lambda: cfg          # this worker only
    q.star_channel.cache_clear()
    rows = q.sweep_arm(arm, p, 'standard')
    return name, arm, p, rows


def tasks():
    out = [(f'angle_{a:.4g}', {'star': {'angle_max': a}}, 'star', p) for a in ANGLES for p in (1e-3, 1e-4)]
    out += [(name, ch, arm, 1e-3) for name, ch in DELIVERY.items() for arm in ('star', 'ft')]
    out += [(name, ch, 'ft', 1e-4) for name, ch in SOURCES.items()] + [('litinski_sources', SOURCES['litinski_sources'], 'ft', 1e-3)]
    return out


def main():
    todo = tasks()
    with mp.get_context('fork').Pool(len(todo)) as pool:
        results = pool.map(task, todo)
    Q = np.geomspace(1e2, 1e9, 71)
    table = []
    for name, arm, p, rows in results:
        for Qb in Q:
            table.append(dict(variant=name, arm=arm, p=p, Q=float(Qb), score=q.capability(rows, Qb)))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'capability.json').write_text(json.dumps(table, indent=1)+'\n')
    (OUT/'config.json').write_text(json.dumps(dict(angles=ANGLES, delivery=DELIVERY, sources=SOURCES, base=BASE), indent=1)+'\n')
    print('done', len(table))


if __name__ == '__main__':
    main()
