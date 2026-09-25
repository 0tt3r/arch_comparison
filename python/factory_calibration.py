"""Distillation factories under the Stim-calibrated surface-code law.

Litinski's factory model (Quantum 3, 205, 2019; vendored in
third_party/litinski_magicstates) prices every patch, merge, and storage step
with one per-patch, per-cycle logical error, 0.1(100p)^((d+1)/2), which is our
reference law (A = 0.1, Lambda = 10 at p = 1e-3). This script replaces it with
the calibrated law of results/surface_law_v1 (A = 0.077, Lambda = 7.4 at 1e-3;
A = 0.125, Lambda = 89 at 1e-4), keeps his injection model (faulty pi/8
rotations with error p/3 of each type), re-optimizes all distances over a grid,
and writes a Pareto menu of (output error, qubits, cycles per output) for the
optimizer and the QUOPS study.

    python factory_calibration.py validate          # reproduce the published rows (reference law)
    python factory_calibration.py run --workers 120 # calibrated grid (resumable)
    python factory_calibration.py menu              # Pareto menu from saved evaluations
"""
import argparse
import itertools
import json
import multiprocessing as mp
import sys
from pathlib import Path

from project_paths import ROOT

VENDOR = ROOT/'python/third_party/litinski_magicstates'
OUT = ROOT/'results/factory_calibration_v1'
FIT = ROOT/'results/surface_law_v1/fit.json'
PHYSICAL = (1e-3, 1e-4)
DPS = 40

# Published rows (Tables 1-2) with the arguments that generate them.
PUBLISHED = [
    ('one', 1e-3, (17, 7, 7)), ('one', 1e-4, (7, 3, 3)), ('one', 1e-4, (9, 3, 3)), ('one', 1e-4, (11, 5, 5)),
    ('two15', 1e-3, (11, 5, 5, 25, 11, 11, 6)),
    ('two15', 1e-3, (13, 5, 5, 29, 11, 13, 6)), ('two15', 1e-3, (17, 7, 7, 41, 17, 17, 6)),
    ('two15', 1e-4, (9, 3, 3, 25, 9, 9, 4)),
    ('two20', 1e-3, (13, 5, 5, 23, 11, 13, 6)), ('two20', 1e-3, (13, 5, 5, 27, 13, 15, 4)),
    ('two20', 1e-4, (9, 3, 3, 15, 7, 9, 4)),
]


def law(name, p):
    """Per-patch, per-cycle logical error p_L(d) as a function of d."""
    if name == 'reference':
        return lambda pphys, d: 0.1*(100*pphys)**((d+1)/2)
    fit = json.loads(FIT.read_text())['fits'][str(p)]
    A, pth = fit['A_patch_union_bound'], fit['effective_p_threshold']
    return lambda pphys, d: A*(pphys/pth)**((d+1)/2)


def modules(name, p):
    """Import the vendored model with its logical-error law replaced."""
    if str(VENDOR) not in sys.path:
        sys.path.insert(0, str(VENDOR))
    import mpmath
    mpmath.mp.dps = DPS
    import definitions
    import onelevel15to1
    import twolevel15to1
    import twolevel20to4
    f = law(name, p)
    for m in (definitions, onelevel15to1, twolevel15to1, twolevel20to4):
        m.plog = f
    twolevel15to1.one_level_15to1_state_memoized.cache_clear()
    return onelevel15to1, twolevel15to1, twolevel20to4


def evaluate(case):
    lawname, kind, p, args = case
    one, two15, two20 = modules(lawname, p)
    f = {'one': one.cost_of_one_level_15to1, 'two15': two15.cost_of_two_level_15to1,
         'two20': two20.cost_of_two_level_20to4}[kind]
    try:
        r = f(p, *args)
    except Exception as error:         # e.g. a level-1 output too poor for level 2
        return dict(law=lawname, kind=kind, p=p, args=list(args), status=f'failed: {error}')
    return dict(law=lawname, kind=kind, p=p, args=list(args), status='complete',
                name=r.name.split(' with ')[0] + ' ' + str(tuple(args)),
                error=float(r.distilled_magic_state_error_rate), qubits=int(r.qubits),
                cycles=float(r.distillation_time_in_cycles), outputs=int(r.n_t_gates_produced_per_distillation))


def odd(a, b):
    return range(a, b+1, 2)


def grid():
    cases = []
    for p in PHYSICAL:
        for dx, dz, dm in itertools.product(odd(5, 31), odd(3, 19), odd(3, 19)):
            if dz <= dx and dm <= dx:
                cases.append(('calibrated', 'one', p, (dx, dz, dm)))
        level1 = [(9, 3, 3), (9, 5, 5), (11, 5, 5), (13, 5, 5), (13, 7, 7), (15, 7, 7), (17, 7, 7), (19, 9, 9)] \
            if p == 1e-3 else [(7, 3, 3), (9, 3, 3), (9, 5, 5), (11, 5, 5), (13, 5, 5)]
        dx2s = odd(17, 55) if p == 1e-3 else odd(11, 35)
        dz2s = odd(7, 25) if p == 1e-3 else odd(5, 15)
        for l1, dx2, dz2, extra, nl1 in itertools.product(level1, dx2s, dz2s, (0, 2), (4, 6)):
            if dz2 <= dx2:
                cases.append(('calibrated', 'two15', p, (*l1, dx2, dz2, dz2+extra, nl1)))
        # 20-to-4 is two orders slower: search around the published designs only.
        for kind, q, args in PUBLISHED:
            if kind == 'two20' and q == p:
                *l1, dx2, dz2, dm2, nl1 = args
                for a, b, c in itertools.product((0, 2, 4, 6, 8), (0, 2, 4), (0, 2, 4)):
                    cases.append(('calibrated', 'two20', p, (*l1, dx2+a, dz2+b, dm2+c, nl1)))
    return cases


def key(v):
    return (v['law'], v['kind'], v['p'], tuple(v['args']))


def run(cases, path, workers):
    OUT.mkdir(parents=True, exist_ok=True)
    done = {}
    if path.exists():
        for line in path.read_text().splitlines():
            v = json.loads(line)
            if v['status'] == 'complete':      # failed designs are retried
                done[key(v)] = v
    todo = [c for c in cases if (c[0], c[1], c[2], tuple(c[3])) not in done]
    # Slowest first: 20-to-4, then two-level 15-to-1.
    todo.sort(key=lambda c: {'two20': 0, 'two15': 1, 'one': 2}[c[1]])
    print(len(done), 'saved;', len(todo), 'to evaluate', flush=True)
    with mp.get_context('fork').Pool(workers) as pool, path.open('a') as out:
        for i, v in enumerate(pool.imap_unordered(evaluate, todo, chunksize=1)):
            out.write(json.dumps(v)+'\n')
            out.flush()
            if i % 200 == 0:
                print(i, 'of', len(todo), flush=True)


def pareto(rows):
    """Nondominated in (error, qubits, cycles per output)."""
    pts = sorted(rows, key=lambda v: (v['error'], v['qubits'], v['cycles']/v['outputs']))
    front = []
    for v in pts:
        c = v['cycles']/v['outputs']
        if not any(u['qubits'] <= v['qubits'] and u['cycles']/u['outputs'] <= c for u in front):
            front.append(v)
    return front


def menu():
    rows = {key(v): v for v in map(json.loads, (OUT/'evaluations.jsonl').read_text().splitlines())
            if v['status'] == 'complete'}.values()
    rows = [v for v in rows if v['status'] == 'complete' and 0 < v['error'] < 1e-6]
    result = {}
    for p in PHYSICAL:
        front = pareto([v for v in rows if v['p'] == p])
        # One decade per bin: the smallest footprint and the smallest space-time per output.
        chosen = {}
        for v in front:
            b = int(-__import__('math').log10(v['error']))
            for tag, score in (('qubits', v['qubits']), ('spacetime', v['qubits']*v['cycles']/v['outputs'])):
                if (b, tag) not in chosen or score < chosen[b, tag][0]:
                    chosen[b, tag] = (score, v)
        picks = {id(v): v for _, v in chosen.values()}.values()
        result[repr(p)] = [[v['name'], v['error'], v['qubits'], round(v['cycles'], 2), v['outputs']]
                           for v in sorted(picks, key=lambda v: -v['error'])]
        print(p, len(front), 'Pareto designs;', len(result[repr(p)]), 'in the menu')
    (OUT/'menu.json').write_text(json.dumps({'litinski_factories': result, 'surface_law': str(FIT.relative_to(ROOT))},
                                            indent=1)+'\n')


def main():
    global FIT, OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('validate', 'run', 'menu'))
    parser.add_argument('--workers', type=int, default=120)
    parser.add_argument('--fit', default=str(FIT.relative_to(ROOT)), help='surface-code law (fit.json)')
    parser.add_argument('--out', default=str(OUT.relative_to(ROOT)), help='output folder')
    args = parser.parse_args()
    FIT, OUT = ROOT/args.fit, ROOT/args.out
    if args.action == 'validate':
        cases = [(lawname, kind, p, a) for kind, p, a in PUBLISHED for lawname in ('reference', 'calibrated')]
        run(cases, OUT/'published.jsonl', min(args.workers, len(cases)))
    elif args.action == 'run':
        run(grid(), OUT/'evaluations.jsonl', args.workers)
    else:
        menu()


if __name__ == '__main__':
    main()
