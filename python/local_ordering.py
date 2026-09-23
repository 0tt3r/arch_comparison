"""Bounded ordering experiment: same Hubbard Hamiltonian, state, and mapping.

Only <=9-site patches run here. No change to the frozen calibration drivers.
The colored permutation keeps the original onsite split, isolating bond order;
the onsite-grouped variant changes that split as a separately named algorithm.
"""
import argparse
import hashlib
import json

import numpy as np

from project_paths import ROOT
from time_calibration import SplitSector


class OrderedSector(SplitSector):
    def colored_edges(self):
        groups = [[] for _ in range(4)]
        for edge in self.edges:
            (x, y), (xx, yy) = [self.sites[i] for i in edge]
            color = x % 2 if y == yy else 2+y % 2
            groups[color].append(edge)
        assert all(len({i for e in group for i in e}) == 2*len(group) for group in groups)
        return groups

    def trotter(self, tau, steps, algorithm='colored'):
        if algorithm in ('plaquette', 'termwise'):
            return super().trotter(tau, steps, algorithm)
        if algorithm not in ('colored', 'onsite_grouped') or self.n > 9:
            raise ValueError('Bounded local-ordering study supports <=9 sites')
        edges = [edge for group in self.colored_edges() for edge in group]
        dt = tau/steps
        sweeps = []
        for spin, masks in enumerate((self.up, self.down)):
            forward = np.eye(len(masks), dtype=complex)
            for edge in edges:
                self._left_bond(forward, edge, dt/2, spin)
            sweeps.append(forward)
        potential = self.potential.reshape(self.shape)
        first_potential = self.local_potential[0] if algorithm == 'colored' else potential
        first = np.exp(-.5j*dt*first_potential)
        middle = np.exp(-1j*dt*(potential-first_potential))
        psi = self.initial.reshape(self.shape).copy()
        for _ in range(steps):
            psi *= first
            psi = sweeps[0]@psi@sweeps[1].T
            psi *= middle
            psi = sweeps[0].T@psi@sweeps[1]
            psi *= first
        return psi.ravel()


def run():
    out = ROOT/'results/physics_sweep_v1'; out.mkdir(parents=True, exist_ok=True)
    rows = []
    for side in (2, 3):
        model = OrderedSector(side, side)
        exact = model.doublon(model.exact(1), (1, 1))
        for algorithm in ('termwise', 'colored', 'onsite_grouped', 'plaquette'):
            values = []
            for r in (8, 16, 32, 64):
                psi = model.trotter(1, r, algorithm)
                error = model.doublon(psi, (1, 1))-exact
                values.append(dict(r=r, signed_error=error, coefficient=abs(error)*r*r,
                                   norm_error=float(abs(np.linalg.norm(psi)-1))))
            a, b = [v['coefficient'] for v in values[-2:]]
            rows.append(dict(patch=f'{side}x{side}', algorithm=algorithm, tau=1,
                             exact_doublon=exact, measurements=values,
                             tail_coefficient=b, relative_tail_spread=abs(a-b)/max(a,b,1e-12)))
            print('ordering', side, algorithm, b, flush=True)
    record = dict(rows=rows, scope='Finite patches at tau=1; no resource curve is recalibrated from this diagnostic',
                  source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in [ROOT/'python/local_ordering.py', ROOT/'python/time_calibration.py', ROOT/'python/small_dynamics.py']})
    (out/'ordering.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    return record


if __name__ == '__main__':
    run()
