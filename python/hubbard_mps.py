"""Termwise Trotter circuits of the U/J = 4 Hubbard model on larger lattices, by fermionic MPS.

The exact sector solver (cpp/hubbard_sector.cpp) stops at 4x4. This evolves
the same termwise product formula as an MPS: one site per lattice site in the
solver's row-snake order (local basis 0, down, up, updown; index 2 n_up + n_down),
fermionic swaps to reach vertical bonds, and each bond's two spin hops applied
as one exact two-site gate (they commute). The step is exactly the solver's:
onsite(dt/2, site 0), bonds forward (dt/2), onsite(dt, other sites), bonds
backward (dt/2), onsite(dt/2, site 0). H = -sum(c^dag c + h.c.) + (U/4) sum Z_up Z_dn,
U = 4, which differs from U n_up n_dn by conserved terms.

The Trotter coefficient C = |D_r - D| r^2 at the central site is taken from two
step counts with the solver-validated r^-2 law: C = (D_r - D_2r) * 4 r^2 / 3,
checked against direct exact values on 3x3 and 4x4.
"""
import argparse
import json
import os
import time
from math import cos, sin

import numpy as np
from scipy.linalg import expm, svd

from project_paths import ROOT

OUT = ROOT/'results/hubbard_mps_v1'
U_OVER_J = 4.


# ---------------------------------------------------------------- local operators

def two_site_operators():
    """4-mode JW (a_up, a_dn, b_up, b_dn); two-site index 4*i_a + i_b with i = 2 n_up + n_dn."""
    a, Z, I = np.array([[0, 1], [0, 0]]), np.diag([1, -1]), np.eye(2)

    def c(k):
        ops = [Z]*k + [a] + [I]*(3-k)
        out = np.eye(1)
        for o in ops:
            out = np.kron(out, o)
        return out
    hop = sum(c(s).T @ c(2+s) + c(2+s).T @ c(s) for s in (0, 1))
    parity = np.array([0, 1, 1, 0])
    fswap = np.zeros((16, 16))
    for ia in range(4):
        for ib in range(4):
            fswap[4*ib+ia, 4*ia+ib] = -1. if parity[ia] and parity[ib] else 1.
    return hop, fswap


HOP, FSWAP = two_site_operators()
ZZ = np.array([1., -1., -1., 1.])             # Z_up Z_dn for index 2 n_up + n_dn
DOUBLON = np.array([0., 0., 0., 1.])


def snake(nx, ny):
    return [(x if y % 2 == 0 else nx-1-x, y) for y in range(ny) for x in range(nx)]


def edges(nx, ny):
    sites = snake(nx, ny)
    index = {s: i for i, s in enumerate(sites)}
    out = []
    for i, (x, y) in enumerate(sites):
        for nb in ((x+1, y), (x, y+1)):
            if nb in index:
                out.append((i, index[nb]))
    return out


# ---------------------------------------------------------------- MPS

class MPS:
    def __init__(self, nx, ny, chi, cutoff):
        self.n = nx*ny
        self.chi, self.cutoff = chi, cutoff
        sites = snake(nx, ny)
        self.A = []
        for x, y in sites:
            t = np.zeros((1, 4, 1), complex)
            t[0, 2 if (x+y) % 2 == 0 else 1, 0] = 1.       # up on even sites, as the solver
            self.A.append(t)
        self.order = list(range(self.n))                 # lattice site at each MPS position
        self.pos = list(range(self.n))
        self.center = 0
        self.discarded = 0.
        self.max_bond = 1

    def move_center(self, k):
        while self.center < k:
            i = self.center
            l, d, r = self.A[i].shape
            q, rr = np.linalg.qr(self.A[i].reshape(l*d, r))
            self.A[i] = q.reshape(l, d, -1)
            self.A[i+1] = np.tensordot(rr, self.A[i+1], axes=(1, 0))
            self.center += 1
        while self.center > k:
            i = self.center
            l, d, r = self.A[i].shape
            q, rr = np.linalg.qr(self.A[i].reshape(l, d*r).T)
            self.A[i] = q.T.reshape(-1, d, r)
            self.A[i-1] = np.tensordot(self.A[i-1], rr.T, axes=(2, 0))
            self.center -= 1

    def gate(self, k, G, direction):
        """Apply a 16x16 gate to MPS positions (k, k+1); leave the center on k+1 (direction > 0) or k."""
        self.move_center(k)
        l, r = self.A[k].shape[0], self.A[k+1].shape[2]
        theta = np.tensordot(self.A[k], self.A[k+1], axes=(2, 0)).reshape(l, 16, r)
        theta = np.einsum('ab,lbr->lar', G, theta).reshape(l*4, 4*r)
        try:
            u, s, vh = svd(theta, full_matrices=False, lapack_driver='gesdd')
        except np.linalg.LinAlgError:
            u, s, vh = svd(theta, full_matrices=False, lapack_driver='gesvd')
        w = s**2/float(s @ s)
        tail = np.append(np.cumsum(w[::-1])[::-1], 0.)  # tail[j]: weight discarded when keeping j values
        keep = min(self.chi, max(1, int(np.argmax(tail <= self.cutoff))))
        self.discarded += float(tail[keep])
        u, s, vh = u[:, :keep], s[:keep]/np.sqrt(float(s[:keep] @ s[:keep])), vh[:keep]
        if direction > 0:
            self.A[k] = u.reshape(l, 4, keep)
            self.A[k+1] = (s[:, None]*vh).reshape(keep, 4, r)
            self.center = k+1
        else:
            self.A[k] = (u*s).reshape(l, 4, keep)
            self.A[k+1] = vh.reshape(keep, 4, r)
            self.center = k
        self.max_bond = max(self.max_bond, keep)

    def swap(self, k, direction):
        self.gate(k, FSWAP, direction)
        a, b = self.order[k], self.order[k+1]
        self.order[k], self.order[k+1] = b, a
        self.pos[a], self.pos[b] = k+1, k

    def bond(self, i, j, G):
        """Apply a symmetric bond gate between lattice sites i and j, restoring the site order."""
        a, b = sorted((self.pos[i], self.pos[j]))
        moves = []
        for k in range(b-1, a, -1):                   # bring position b next to a
            self.swap(k, -1)
            moves.append(k)
        self.gate(a, G, +1)
        for k in reversed(moves):
            self.swap(k, +1)

    def onsite(self, sites, phase):
        for i in sites:
            self.A[self.pos[i]] = self.A[self.pos[i]]*phase[None, :, None]

    def local(self, i, diag):
        k = self.pos[i]
        self.move_center(k)
        t = self.A[k]
        return float(np.real(np.einsum('lar,a,lar->', t.conj(), diag, t)))


def termwise(nx, ny, tau, r, chi, cutoff, log=None):
    psi = MPS(nx, ny, chi, cutoff)
    dt = tau/r
    bonds = edges(nx, ny)
    G = expm(1j*(dt/2)*HOP)
    half = np.exp(-1j*(dt/2)*(U_OVER_J/4)*ZZ)
    full = np.exp(-1j*dt*(U_OVER_J/4)*ZZ)
    rest = range(1, nx*ny)
    start = time.time()
    for step in range(r):
        psi.onsite([0], half)
        for i, j in bonds:
            psi.bond(i, j, G)
        psi.onsite(rest, full)
        for i, j in reversed(bonds):
            psi.bond(i, j, G)
        psi.onsite([0], half)
        if log:
            print(f'{nx}x{ny} r={r} chi={chi} step {step+1}/{r} bond {psi.max_bond} '
                  f'discarded {psi.discarded:.2e} {time.time()-start:.0f}s', file=log, flush=True)
    sites = snake(nx, ny)
    doublons = {f'{x},{y}': psi.local(i, DOUBLON) for i, (x, y) in enumerate(sites)}
    x, y = nx//2, ny//2
    return dict(nx=nx, ny=ny, tau=tau, r=r, chi=chi, cutoff=cutoff, site=[x, y],
                doublon=doublons[f'{x},{y}'], doublons=doublons, discarded_weight=psi.discarded,
                max_bond=psi.max_bond, seconds=time.time()-start)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--patch', default='4x4')
    p.add_argument('--tau', type=float, default=1.)
    p.add_argument('--steps', type=int, nargs='+', default=[8, 16])
    p.add_argument('--chi', type=int, default=512)
    p.add_argument('--cutoff', type=float, default=1e-12)
    args = p.parse_args()
    nx, ny = map(int, args.patch.split('x'))
    OUT.mkdir(parents=True, exist_ok=True)
    for r in args.steps:
        name = OUT/f'{nx}x{ny}_tau{args.tau:g}_r{r}_chi{args.chi}_cut{args.cutoff:g}.json'
        if name.exists():
            continue
        with open(name.with_suffix('.log'), 'w') as log:
            row = termwise(nx, ny, args.tau, r, args.chi, args.cutoff, log)
        row.update(threads=os.getenv('OMP_NUM_THREADS'), pbs_jobid=os.getenv('PBS_JOBID'))
        tmp = name.with_suffix('.tmp')
        tmp.write_text(json.dumps(row, indent=1)+'\n')
        tmp.replace(name)
        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
