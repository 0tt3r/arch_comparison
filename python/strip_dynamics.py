"""Exact Trotter coefficients on strips beyond 4x4 (U/J = 4, Neel state).

A fermionic MPS cannot reach 5x5 or 6x6 at the needed precision: the quench
builds near-volume-law entanglement by tau = 1 (4x3 at chi = 256 already
discards 43% of the weight). The sector solver is exact up to 18 sites, so we
extend the measured patches along strips of width 2 and 3 (up to 8x2 and 6x3)
to test whether the termwise coefficient converges in length at fixed width.
Uses large_dynamics.solve with cpp/hubbard_strip.cpp, a copy of the
hash-pinned solver with the site limit raised from 16 to 18.

    python strip_dynamics.py --threads 128        # on a PBS node; resumable
    python strip_dynamics.py --analyze
"""
import argparse
import json
import os

import large_dynamics as ld
from project_paths import ROOT

ld.SOURCE = ROOT/'cpp/hubbard_strip.cpp'
ld.BINARY = ROOT/'build/hubbard_strip'
OUT = ROOT/'results/strip_dynamics_v1'
PATCHES = ['2x2', '3x2', '4x2', '5x2', '6x2', '7x2', '8x2', '3x3', '4x3', '5x3', '6x3', '4x4']
STEPS = {1.: [8, 16, 32], 2.: [16, 32]}


def cases():
    out = []
    for tau, steps in STEPS.items():
        for patch in PATCHES:
            for algorithm in ('exact', 'termwise', 'plaquette'):
                for r in ([0] if algorithm == 'exact' else steps):
                    out.append((patch, tau, algorithm, r))
    # Smallest first, so the campaign produces usable rows early.
    return sorted(out, key=lambda c: (eval(c[0].replace('x', '*')), c[1]))


def sites(patch):
    nx, ny = map(int, patch.split('x'))
    return [(x if y % 2 == 0 else nx-1-x, y) for y in range(ny) for x in range(nx)]


def analyze():
    rows = []
    for patch in PATCHES:
        nx, ny = map(int, patch.split('x'))
        order = sites(patch)
        centers = sorted({(nx//2, ny//2), ((nx-1)//2, (ny-1)//2)})
        for tau, steps in STEPS.items():
            path = OUT/f'{patch}_tau{tau:g}_exact_r0_t{THREADS}.json'
            if not path.exists():
                continue
            exact = json.loads(path.read_text())['result']['doublons']
            for algorithm in ('termwise', 'plaquette'):
                for r in steps:
                    p = OUT/f'{patch}_tau{tau:g}_{algorithm}_r{r}_t{THREADS}.json'
                    if not p.exists():
                        continue
                    D = json.loads(p.read_text())['result']['doublons']
                    err = [abs(a-b)*r*r for a, b in zip(D, exact)]
                    rows.append(dict(patch=patch, nx=nx, ny=ny, tau=tau, algorithm=algorithm, r=r,
                                     central=max(err[order.index(c)] for c in centers),
                                     centers=[list(c) for c in centers], maximum=max(err)))
    (OUT/'coefficients.json').write_text(json.dumps(rows, indent=1)+'\n')
    for tau in STEPS:
        for algorithm in ('termwise', 'plaquette'):
            r = STEPS[tau][1]
            print(f'tau={tau:g} {algorithm:9s} r={r}:',
                  '  '.join(f"{v['patch']} {v['central']:.3f}" for v in rows
                            if v['tau'] == tau and v['algorithm'] == algorithm and v['r'] == r))
    return rows


THREADS = int(os.getenv('STRIP_THREADS', '128'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analyze', action='store_true')
    parser.add_argument('--compile', action='store_true')
    args = parser.parse_args()
    if args.analyze:
        analyze()
        return
    if args.compile or not ld.BINARY.exists():
        print(json.dumps(ld.compile_solver(), indent=2), flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    for patch, tau, algorithm, r in cases():
        nx, ny = map(int, patch.split('x'))
        ld.solve(nx, ny, tau, algorithm, r, OUT, THREADS)
    analyze()


if __name__ == '__main__':
    main()
