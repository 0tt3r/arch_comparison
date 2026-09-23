"""Free-fermion (U = 0) check of the Trotter coefficient's dependence on lattice size.

At U = 0 every gate of the termwise and plaquette circuits of hubbard_circuits is a
Majorana bilinear rotation, so the circuits are Gaussian and run exactly in
O(N) per gate on the 2n x 2n Majorana covariance matrix. This scans the doublon
coefficient C = |D_r - D_exact| r^2 at the central site versus L, the quantity
the manuscript calibrates at U/J = 4 on 3x3, 4x3, and 4x4 patches.

Scope: hopping-only Trotter error and Jordan-Wigner ordering; the interaction
commutators that dominate at U/J = 4 are absent. Validated against a state
vector on 2x2 and 3x3.
"""
import argparse
import json
import multiprocessing as mp
from functools import lru_cache
from math import pi

import numpy as np
from scipy.linalg import expm

import hubbard_circuits as hc
from project_paths import ROOT

OUT = ROOT/'results/free_fermion_trotter_v1'


# ---------------------------------------------------------------- circuits (any L)

def circuit(L, tau, steps, algorithm):
    """hubbard_circuits.build_workload at U = 0, without its even-L restriction.

    Returns the rotation/Z operations (onsite rotations have zero angle and are
    dropped) and the occupied modes of the Neel state.
    """
    modes, dt = hc.mode_map(L), tau/steps
    occupied = [modes[x, y, (x+y) % 2] for x, y in hc.snake_sites(L)]
    if algorithm == 'termwise':
        terms = []
        for x, y in hc.snake_sites(L):
            for nx, ny in ((x+1, y), (x, y+1)):
                if nx < L and ny < L:
                    for spin in (0, 1):
                        for endpoint in ('XX', 'YY'):
                            terms.append(hc.parity_rotation(modes[x, y, spin], modes[nx, ny, spin], endpoint, -.5))
        # Zero-angle onsite terms keep the builder's merge pattern at the Strang centers.
        onsite = [hc.Operation('rotate', (2*i, 2*i+1), 0., 'ZZ', 'dynamic') for i in range(L*L)]
        terms = onsite[:1] + terms + onsite[1:]
        half = [hc.Operation(op.name, op.qubits, op.angle*dt/2, op.paulis, op.role) for op in terms]
        ops = hc._merge((half + half[::-1])*steps)
    else:
        ops = []
        for _ in range(steps):
            for name, duration in (('pink', dt/2), ('gold', dt), ('pink', dt/2)):
                for tile in hc.tiles(L, int(name == 'gold')):
                    for spin in (0, 1):
                        idx = [modes[x, y, spin] for x, y in tile]
                        if len(idx) == 4:
                            ops.extend(hc.plaquette_evolution(idx, duration))
                        else:
                            ops.extend(hc.parity_rotation(*idx, e, -duration/2) for e in ('XX', 'YY'))
    return [op for op in ops if op.name == 'z' or (op.name == 'rotate' and op.paulis != 'ZZ')], occupied


def check_against_builder(L, tau, steps, algorithm):
    """For even L the operations match build_workload(u_over_j=0) exactly."""
    ops, occupied = circuit(L, tau, steps, algorithm)
    built, _ = hc.build_workload(L, tau, steps, algorithm, u_over_j=0)
    kept = [op for op in built if op.name == 'z' or (op.name == 'rotate' and op.paulis != 'ZZ')]
    if algorithm == 'plaquette':
        # The builder merges adjacent pink half steps across step boundaries.
        return sum(op.angle for op in kept if op.role == 'dynamic') == sum(op.angle for op in ops if op.role == 'dynamic')
    return kept == ops and sorted(op.qubits[0] for op in built if op.name == 'x') == sorted(occupied)


# ---------------------------------------------------------------- Majorana algebra

@lru_cache(None)
def bilinear(letters):
    """P = i s g_a g_b for a two-endpoint JW string; returns (offset_a, offset_b, s).

    Found numerically on three modes: g_2k = Z_<k X_k, g_2k+1 = Z_<k Y_k.
    """
    I, X, Y, Z = np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])
    single = dict(I=I, X=X, Y=Y, Z=Z)

    def kron(ps):
        out = np.eye(1)
        for p in ps:
            out = np.kron(out, single[p])
        return out
    n = 3
    gam = []
    for k in range(n):
        for last in 'XY':
            gam.append(kron('Z'*k + last + 'I'*(n-k-1)))
    target = kron(letters[0] + 'Z' + letters[1]) if len(letters) == 2 else kron('I' + letters + 'I')
    lo, hi = (0, 2) if len(letters) == 2 else (1, 1)
    for a in (2*lo, 2*lo+1):
        for b in (2*hi, 2*hi+1):
            if a >= b:
                continue
            for s in (1, -1):
                if np.allclose(target, 1j*s*gam[a]@gam[b]):
                    return a-2*lo, b-2*hi, s
    raise ValueError(letters)


def apply(M, op):
    """M -> R M R^T for exp(-i theta P), P = i s g_a g_b: a rotation by 2 s theta."""
    if op.name == 'z':
        j = op.qubits[0]
        a, b, s, theta = 2*j, 2*j+1, *bilinear('Z')[2:], pi/2
    else:
        lo, hi = op.qubits[0], op.qubits[-1]
        oa, ob, s = bilinear(op.paulis[0] + op.paulis[-1])
        a, b, theta = 2*lo+oa, 2*hi+ob, op.angle
    c, sn = np.cos(2*s*theta), np.sin(2*s*theta)
    ra, rb = M[a].copy(), M[b].copy()
    M[a], M[b] = c*ra + sn*rb, c*rb - sn*ra
    ca, cb = M[:, a].copy(), M[:, b].copy()
    M[:, a], M[:, b] = c*ca + sn*cb, c*cb - sn*ca


def fock_covariance(n, occupied):
    """M_pq = i<g_p g_q> (p != q); Z_k = -i g_2k g_2k+1 so M_2k,2k+1 = -<Z_k> = 2 n_k - 1."""
    M = np.zeros((2*n, 2*n))
    for k in range(n):
        v = 1. if k in occupied else -1.
        M[2*k, 2*k+1], M[2*k+1, 2*k] = v, -v
    return M


def doublon(M, up, down):
    """<n_up n_down> by Wick's theorem, n_k = (1 + i g_2k g_2k+1)/2."""
    G = lambda p, q: (1. if p == q else 0.) - 1j*M[p, q]
    a, b, c, d = 2*up, 2*up+1, 2*down, 2*down+1
    four = G(a, b)*G(c, d) - G(a, c)*G(b, d) + G(a, d)*G(b, c)
    return float(np.real(.25*(1 + 1j*G(a, b) + 1j*G(c, d) - four)))


def trotter_doublons(L, tau, steps, algorithm):
    ops, occupied = circuit(L, tau, steps, algorithm)
    n = 2*L*L
    M = fock_covariance(n, set(occupied))
    for op in ops:
        apply(M, op)
    modes = hc.mode_map(L)
    return np.array([[doublon(M, modes[x, y, 0], modes[x, y, 1]) for x in range(L)] for y in range(L)])


def exact_doublons(L, tau):
    """Single-particle evolution; spin sectors are independent at U = 0."""
    N = L*L
    index = {(x, y): x + L*y for y in range(L) for x in range(L)}
    h = np.zeros((N, N))
    for (x, y), i in index.items():
        for nb in ((x+1, y), (x, y+1)):
            if nb in index:
                h[i, index[nb]] = h[index[nb], i] = -1.
    U = expm(-1j*tau*h)
    dens = {}
    for spin in (0, 1):
        occ = [i for (x, y), i in index.items() if (x+y) % 2 == spin]
        dens[spin] = (np.abs(U[:, occ])**2).sum(axis=1)
    D = dens[0]*dens[1]
    return np.array([[D[index[x, y]] for x in range(L)] for y in range(L)])


# ---------------------------------------------------------------- validation

def statevector_doublon(L, tau, steps, algorithm):
    ops, occupied = circuit(L, tau, steps, algorithm)
    n = 2*L*L
    idx = np.arange(2**n, dtype=np.int64)
    psi = np.zeros(2**n, complex)
    psi[sum(1 << k for k in occupied)] = 1
    bit = lambda q: (idx >> q) & 1
    for op in ops:
        if op.name == 'z':
            psi = psi*(1-2*bit(op.qubits[0]))
            continue
        flip, phase = 0, np.ones(2**n, complex)
        for q, p in zip(op.qubits, op.paulis):
            if p == 'Z':
                phase *= 1-2*bit(q)
            elif p == 'X':
                flip |= 1 << q
            else:  # Y|0> = i|1>, Y|1> = -i|0>: output bit 1 carries +i
                flip |= 1 << q
                phase *= np.where(bit(q) == 1, 1j, -1j)
        P = phase*psi[idx ^ flip]
        psi = np.cos(op.angle)*psi - 1j*np.sin(op.angle)*P
    modes = hc.mode_map(L)
    x = y = L//2
    both = (bit(modes[x, y, 0]) & bit(modes[x, y, 1])).astype(bool)
    return float((np.abs(psi[both])**2).sum())


def validate():
    rows = []
    for L, tau, steps in ((2, 1., 3), (3, 1., 2), (3, .7, 3)):
        for algorithm in ('termwise', 'plaquette'):
            sv = statevector_doublon(L, tau, steps, algorithm)
            ff = trotter_doublons(L, tau, steps, algorithm)[L//2, L//2]
            rows.append(dict(L=L, tau=tau, r=steps, algorithm=algorithm, statevector=sv, gaussian=ff, difference=abs(sv-ff)))
            assert abs(sv-ff) < 1e-10, rows[-1]
    for L in (2, 4):
        for algorithm in ('termwise', 'plaquette'):
            assert check_against_builder(L, 1., 3, algorithm), (L, algorithm)
    return rows


# ---------------------------------------------------------------- scan

def scan_case(args):
    algorithm, tau, L, steps = args
    t = L/4 if tau == 'L/4' else tau
    exact = exact_doublons(L, t)
    c = L//2
    rows = []
    for r in steps:
        err = trotter_doublons(L, t, r, algorithm) - exact
        rows.append(dict(algorithm=algorithm, family=tau if tau == 'L/4' else 'fixed', tau=t, L=L, r=r,
                         central_coefficient=abs(err[c, c])*r*r,
                         max_site_coefficient=float(np.abs(err).max())*r*r,
                         exact_central_doublon=float(exact[c, c])))
    return rows


def scan(sizes, taus, steps, workers):
    cases = [(a, t, L, steps) for a in ('termwise', 'plaquette') for t in taus for L in sizes
             if t != 'L/4' or L >= 4]
    cases.sort(key=lambda v: -v[2])
    with mp.get_context('fork').Pool(workers) as pool:
        parts = pool.map(scan_case, cases, chunksize=1)
    rows = sorted((v for part in parts for v in part), key=lambda v: (v['algorithm'], v['family'], v['tau'], v['L'], v['r']))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-size', type=int, default=20)
    parser.add_argument('--workers', type=int, default=64)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    checks = validate()
    rows = scan(list(range(2, args.max_size+1)), [1., 2., 4., 'L/4'], [16, 32], args.workers)
    (OUT/'coefficients.json').write_text(json.dumps(dict(validation=checks, rows=rows), indent=1)+'\n')


if __name__ == '__main__':
    main()
