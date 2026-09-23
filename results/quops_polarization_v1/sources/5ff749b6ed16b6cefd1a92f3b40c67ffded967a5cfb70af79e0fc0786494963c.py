"""Exact QUOPS process-polarization checks (standalone; no project imports).

The full PTM backend computes each concrete arbitrary-angle circuit exactly.
The ensemble backend integrates axes, angles, and oriented random matchings
analytically; it computes the SAME ensemble mean, not individual circuits.
All arithmetic is float64. No trajectories, state-fidelity proxy, or truncation.

Run with the qb_gsee interpreter and -B. Outputs stay in the v1 result directory.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import socket
import sys
import time

# Avoid implicit BLAS teams and writes outside the authorized output directory.
ROOT = Path(os.environ.get("QUOPS_ROOT", str(Path(__file__).resolve().parents[1])))
OUT = ROOT / "results/quops_polarization_v1"
os.environ.setdefault("MPLCONFIGDIR", str(OUT / "mplconfig"))
import numpy as np
from numba import njit

SEED = 26092303
PS = (0.001, 0.003, 0.01)
TARGET = 0.195
PAULIS = np.array([np.eye(2), [[0, 1], [1, 0]],
                   [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], dtype=complex)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def csv_write(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def npz_write(path, **arrays):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("wb") as f:
        np.savez_compressed(f, **arrays)
    temp.replace(path)


def local_cnot():
    # Independent, tiny Hilbert-space construction avoids phase-table mistakes.
    u = np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                  [0, 0, 0, 1], [0, 0, 1, 0]])
    basis = [np.kron(a, b) for a in PAULIS for b in PAULIS]
    perm, signs = [], []
    for a in basis:
        coeff = np.array([np.trace(b @ u @ a @ u.T).real / 4 for b in basis])
        k = int(np.argmax(abs(coeff)))
        assert abs(coeff[k]) == 1
        perm.append(k)
        signs.append(coeff[k])
    return np.array(perm), np.array(signs)


CNOT_MAP, CNOT_SIGN = local_cnot()


@lru_cache(maxsize=None)
def geometry(w):
    ids = np.arange(4 ** w, dtype=np.int64)
    digits = np.array([(ids // 4 ** q) % 4 for q in range(w)])
    weight = (digits != 0).sum(axis=0)
    counts = np.array([(digits == p).sum(axis=0) for p in (1, 2, 3)]).T
    _, groups = np.unique(counts, axis=0, return_inverse=True)
    sizes = np.bincount(groups)
    return ids, digits, weight, groups, sizes


def matching_map(w, order):
    ids, digits, weight, groups, sizes = geometry(w)
    perm = ids.copy()
    signs = np.ones(len(ids))
    active = np.zeros(len(ids), dtype=np.int64)
    for q, r in zip(order[:2*(w//2):2], order[1:2*(w//2):2]):
        old = 4 * digits[q] + digits[r]
        new = CNOT_MAP[old]
        perm += (new // 4 - digits[q]) * 4 ** int(q)
        perm += (new % 4 - digits[r]) * 4 ** int(r)
        signs *= CNOT_SIGN[old]
        active += old != 0
    idle = (digits[order[-1]] != 0).astype(np.int64) if w % 2 else np.zeros(len(ids), dtype=np.int64)
    return perm, signs, active, idle


def draw_circuit(w, depth, replica):
    # Seeds fixed by width and replica; p and idle comparisons share circuits.
    rng = np.random.default_rng(np.random.SeedSequence([SEED, w, replica]))
    axes = rng.integers(1, 4, size=(depth, w), dtype=np.int64)
    angles = rng.uniform(0, 2*np.pi, size=(depth, w))
    orders = np.array([rng.permutation(w) for _ in range(depth)])
    return axes, angles, orders


@njit
def rotate_overlap(m, q, axis, angle):
    """M <- R M R^T, where M = R_noisy R_ideal^T (unnormalized)."""
    n = len(m)
    stride = 4 ** q
    b, c = (axis % 3) + 1, ((axis + 1) % 3) + 1
    co, si = math.cos(2*angle), math.sin(2*angle)
    for base in range(0, n, 4*stride):
        for offset in range(stride):
            i, k = base + b*stride + offset, base + c*stride + offset
            for j in range(n):
                a, z = m[i, j], m[k, j]
                m[i, j], m[k, j] = co*a-si*z, si*a+co*z
    for i in range(n):
        for base in range(0, n, 4*stride):
            for offset in range(stride):
                j, k = base + b*stride + offset, base + c*stride + offset
                a, z = m[i, j], m[i, k]
                m[i, j], m[i, k] = co*a-si*z, si*a+co*z


@njit
def noisy_cnot_layer(m, scratch, perm, signs, factor):
    n = len(m)
    for i in range(n):
        k = perm[i]
        f = factor[i]*signs[i]
        for j in range(n):
            scratch[i, j] = f*signs[j]*m[k, perm[j]]


def exact_circuit(w, axes, angles, orders, p, idle=True):
    """All prefix depths. Full PTM overlap; O(16**w) memory and gate work."""
    n = 4 ** w
    m, scratch = np.eye(n), np.empty((n, n))
    weight = geometry(w)[2]
    single = (1 - 4*(p/10)/3) ** weight
    curve = np.empty(len(axes))
    for layer, (ax, theta, order) in enumerate(zip(axes, angles, orders)):
        for q in range(w):
            rotate_overlap(m, q, int(ax[q]), float(theta[q]))
        perm, signs, active, waiting = matching_map(w, order)
        # Disjoint CNOTs are self-inverse. Pair support is unchanged by CNOT.
        factor = single[perm] * (1 - 16*p/15) ** active
        if idle:
            factor *= (1 - 4*(p/10)/3) ** waiting
        noisy_cnot_layer(m, scratch, perm, signs, factor)
        m, scratch = scratch, m
        curve[layer] = (np.trace(m) - 1) / (n - 1)
    return curve


@njit
def average_rotation(v, w, axes):
    # axes=0 integrates the axis as well; otherwise integrate only its angle.
    n = len(v)
    for q in range(w):
        stride = 4 ** q
        for base in range(0, n, 4*stride):
            for offset in range(stride):
                i = base + offset
                if axes[q] == 0:
                    a, b, c = v[i+stride], v[i+2*stride], v[i+3*stride]
                    total = (a+b+c)/6
                    v[i+stride] = a/2 + total
                    v[i+2*stride] = b/2 + total
                    v[i+3*stride] = c/2 + total
                else:
                    b, c = (axes[q] % 3)+1, ((axes[q]+1) % 3)+1
                    mean = (v[i+b*stride] + v[i+c*stride])/2
                    v[i+b*stride] = mean
                    v[i+c*stride] = mean


@njit
def project_permutations(v, groups, sizes):
    sums = np.zeros(len(sizes))
    for i in range(len(v)):
        sums[groups[i]] += v[i]
    for i in range(len(v)):
        v[i] = sums[groups[i]]/sizes[groups[i]]


def exact_ensemble(w, p, depth, idle=True):
    """Analytically integrate all angles, axes, and uniform oriented matchings.

    Starting with M=I, the angle integral keeps E[M] diagonal. A random-axis
    rotation mixes its three nonidentity entries with diag=2/3, offdiag=1/6.
    With qubit-permutation-invariant input, a fixed matching followed by the
    permutation projection equals averaging all oriented matchings. The latter
    projection uses histograms (nX,nY,nZ), never a Pauli-weight approximation.
    """
    _, _, weight, groups, sizes = geometry(w)
    perm, _, active, waiting = matching_map(w, np.arange(w))
    factors = (1-4*(p/10)/3)**weight[perm] * (1-16*p/15)**active
    if idle:
        factors *= (1-4*(p/10)/3)**waiting
    v = np.ones(4**w)
    axes = np.zeros(w, dtype=np.int64)
    curve = np.empty(depth)
    for layer in range(depth):
        average_rotation(v, w, axes)
        v = v[perm]*factors
        project_permutations(v, groups, sizes)
        curve[layer] = (v.sum()-1)/(len(v)-1)
    return curve


def conditional_ensemble(w, axes, orders, p, idle=True):
    """Exact angle-integrated mean conditional on sampled axes and matchings."""
    weight = geometry(w)[2]
    single = (1-4*(p/10)/3)**weight
    v = np.ones(4**w)
    curve = np.empty(len(axes))
    for layer, (ax, order) in enumerate(zip(axes, orders)):
        average_rotation(v, w, ax)
        perm, _, active, waiting = matching_map(w, order)
        factor = single[perm]*(1-16*p/15)**active
        if idle:
            factor *= (1-4*(p/10)/3)**waiting
        v = v[perm]*factor
        curve[layer] = (v.sum()-1)/(len(v)-1)
    return curve


def prediction(w, p, depth, idle=True):
    depth = np.asarray(depth)
    n1, n2, ni = w*depth, (w//2)*depth, (w%2)*depth*int(idle)
    exposure = p*(n2+(n1+ni)/10)
    local = (1-4*(p/10)/3)**(n1+ni)*(1-16*p/15)**n2
    aw = 4**w/(4**w-1)
    embedded = (1-aw*(p/10))**(n1+ni)*(1-aw*p)**n2
    return exposure, local, embedded


def depth_for(w, p):
    # Start below the target conservatively; validate the actual endpoint later.
    return int(math.ceil(-math.log(0.18)/(p*(w//2 + w/10))))


def cirq_choi_curve(w, axes, angles, orders, p, idle=True):
    """Independent complex128 Choi-state propagation with Cirq channels."""
    import cirq
    d = 2**w
    qubits = cirq.LineQubit.range(2*w)
    vector = np.eye(d, dtype=complex).reshape(-1)/math.sqrt(d)
    rho = np.outer(vector, vector.conj())
    ideal = np.eye(d, dtype=complex)
    curve = []
    sim = cirq.DensityMatrixSimulator(dtype=np.complex128)
    for ax, theta, order in zip(axes, angles, orders):
        gates = [cirq.MatrixGate(math.cos(t)*PAULIS[0]-1j*math.sin(t)*PAULIS[a])(qubits[q])
                 for q, (a, t) in enumerate(zip(ax, theta))]
        pairs = [cirq.CNOT(qubits[int(q)], qubits[int(r)])
                 for q, r in zip(order[:2*(w//2):2], order[1:2*(w//2):2])]
        circuit = cirq.Circuit(cirq.Moment(gates),
                               cirq.Moment(cirq.depolarize(p/10)(q) for q in qubits[:w]),
                               cirq.Moment(pairs),
                               cirq.Moment(cirq.depolarize(p, n_qubits=2).on(*g.qubits) for g in pairs))
        if idle and w % 2:
            circuit.append(cirq.depolarize(p/10)(qubits[int(order[-1])]))
        rho = sim.simulate(circuit, initial_state=rho, qubit_order=qubits).final_density_matrix
        unitary = cirq.Circuit(cirq.Moment(gates), cirq.Moment(pairs)).unitary(qubit_order=qubits[:w])
        ideal = unitary @ ideal
        target = ideal.reshape(-1)/math.sqrt(d)
        f = np.vdot(target, rho @ target).real
        curve.append((d*d*f-1)/(d*d-1))
    return np.array(curve)


def validate():
    import cirq
    checks = []
    for w in (2, 3):
        ax, theta, order = draw_circuit(w, 5, 911)
        for p in (0., .01, .13):
            for idle in (False, True):
                a = exact_circuit(w, ax, theta, order, p, idle)
                b = cirq_choi_curve(w, ax, theta, order, p, idle)
                err = float(np.max(abs(a-b)))
                assert err < 2e-12, (w, p, idle, err)
                checks.append(dict(check="PTM_vs_Cirq_Choi", w=w, p=p, idle=idle, max_abs=err))
    # Exact quadrature: cos(2theta),sin(2theta) first/second moments are integrated
    # by theta in {0, pi/4, pi/2, 3pi/4}. Exhaust all 4^6 angle choices.
    w, depth, p = 2, 3, .13
    ax, _, order = draw_circuit(w, depth, 912)
    total = np.zeros(depth)
    for values in itertools.product(range(4), repeat=w*depth):
        theta = np.array(values).reshape(depth, w)*np.pi/4
        total += exact_circuit(w, ax, theta, order, p)
    quadrature = total/(4**(w*depth))
    conditional = conditional_ensemble(w, ax, order, p)
    err = float(np.max(abs(quadrature-conditional)))
    assert err < 3e-12, err
    checks.append(dict(check="angle_integral_vs_4096_exact_circuits", max_abs=err))
    # Explicit matching average, on a permutation-invariant input, equals the
    # orbit projection. All 3! and 4! ordered-qubit permutations are enumerated.
    for w in (3, 4):
        ids, _, weight, groups, sizes = geometry(w)
        rng = np.random.default_rng(SEED+w)
        v = rng.uniform(.1, 1, len(sizes))[groups]
        a = np.zeros(len(ids))
        for order in itertools.permutations(range(w)):
            perm, _, active, waiting = matching_map(w, np.array(order))
            f = (1-4*p/30)**weight[perm]*(1-16*p/15)**active*(1-4*p/30)**waiting
            a += v[perm]*f/math.factorial(w)
        perm, _, active, waiting = matching_map(w, np.arange(w))
        b = v[perm]*(1-4*p/30)**weight[perm]*(1-16*p/15)**active*(1-4*p/30)**waiting
        project_permutations(b, groups, sizes)
        err = float(np.max(abs(a-b)))
        assert err < 1e-13, err
        checks.append(dict(check="matching_integral_vs_all_permutations", w=w, max_abs=err))
    # Single noisy gate embedded in w qubits has process infidelity exactly p.
    for w in range(2, 9):
        weight = geometry(w)[2]
        digit = (np.arange(4**w) % 4 != 0)
        gamma = (((1-4*p/3)**digit).sum()-1)/(4**w-1)
        expected = 1 - 4**w*p/(4**w-1)
        assert abs(gamma-expected) < 1e-14
        checks.append(dict(check="embedded_single_gate_normalization", w=w, max_abs=float(abs(gamma-expected))))
    result = dict(status="passed", cirq_version=cirq.__version__, checks=checks,
                  source_sha256=sha(__file__))
    json_write(OUT/"validation.json", result)
    print(json.dumps(result), flush=True)


def run_ensemble():
    rows = []
    start = time.time()
    for w in range(2, 9):
        for p in PS:
            depth = depth_for(w, p)
            on = exact_ensemble(w, p, depth)
            off = exact_ensemble(w, p, depth, False) if w % 2 else on
            for idle, curve in ((True, on), (False, off)):
                end = int(np.flatnonzero(curve <= TARGET)[0])+1
                for layer, gamma in enumerate(curve[:end], 1):
                    exposure, local, embedded = prediction(w, p, layer, idle)
                    rows.append(dict(method="exact_full_ensemble_mean", w=w, p=p, idle=int(idle),
                                     D=layer, gamma=float(gamma), F_proc=float((1+(4**w-1)*gamma)/4**w),
                                     n_1q=w*layer, n_2q=(w//2)*layer,
                                     n_idle_10ns=(w%2)*layer*int(idle),
                                     exposure=float(exposure), product_local=float(local),
                                     product_embedded=float(embedded), ratio_local=float(gamma/local),
                                     ratio_embedded=float(gamma/embedded)))
            print(f"ensemble w={w} p={p:g} D={end} gamma={on[end-1]:.6g}", flush=True)
    csv_write(OUT/"ensemble.csv", rows)
    json_write(OUT/"ensemble_run.json", dict(status="complete", seconds=time.time()-start,
                host=socket.gethostname(), python=sys.executable, python_version=platform.python_version(),
                numpy=np.__version__, source_sha256=sha(__file__), rows=len(rows)))


def run_sample(w, replica, max_depth=None, method="exact_ptm"):
    start = time.time()
    source_hash = sha(__file__)
    depth = max_depth or max(depth_for(w, p) for p in PS)
    axes, angles, orders = draw_circuit(w, depth, replica)
    circuit_path = OUT/"circuits"/f"w{w}_r{replica:02d}_d{depth}.npz"
    if not circuit_path.exists():
        npz_write(circuit_path, axes=axes, angles=angles, orders=orders)
    records = []
    for p in PS:
        stop = min(depth, depth_for(w, p))
        for idle in (True, False) if w % 2 else (True,):
            if method == "exact_ptm":
                curve = exact_circuit(w, axes[:stop], angles[:stop], orders[:stop], p, idle)
            else:
                curve = conditional_ensemble(w, axes[:stop], orders[:stop], p, idle)
            for layer, gamma in enumerate(curve, 1):
                records.append(dict(method=method, w=w, replica=replica, p=p, idle=int(idle),
                                    D=layer, gamma=float(gamma)))
    stem = f"{method}_w{w}_r{replica:02d}_d{depth}"
    csv_write(OUT/"raw"/(stem+".csv"), records)
    receipt = dict(status="complete", method=method, w=w, replica=replica, max_depth=depth,
                   rows=len(records), seconds=time.time()-start, host=socket.gethostname(),
                   source_sha256=source_hash, circuit_sha256=sha(circuit_path),
                   raw_sha256=sha(OUT/"raw"/(stem+".csv")),
                   seed=[SEED, w, replica])
    json_write(OUT/"receipts"/(stem+".json"), receipt)
    print(json.dumps(receipt), flush=True)


def benchmark():
    rows = []
    for w, depth in ((2, 2), (4, 10), (5, 4), (6, 2)):
        ax, theta, orders = draw_circuit(w, depth, 999)
        start = time.time()
        curve = exact_circuit(w, ax, theta, orders, .01)
        elapsed = time.time()-start
        rows.append(dict(w=w, D=depth, seconds=elapsed, seconds_per_layer=elapsed/depth,
                         gamma=float(curve[-1]), matrix_bytes=2*8*16**w))
        print(rows[-1], flush=True)
    json_write(OUT/"benchmark.json", rows)


def campaign(workers, replicas, exact_max_width):
    """Disjoint atomic case outputs; restart only from verified receipts."""
    source_hash = sha(__file__)
    jobs = [(w, r, None, "exact_ptm") for w in range(2, exact_max_width+1)
            for r in range(replicas)]
    jobs += [(w, r, None, "exact_angle_mean") for w in range(2, 9)
             for r in range(replicas)]
    rng = np.random.default_rng(SEED)
    rng.shuffle(jobs)
    pending = []
    for w, r, depth, method in jobs:
        depth = max(depth_for(w, p) for p in PS)
        circuit_path = OUT/"circuits"/f"w{w}_r{r:02d}_d{depth}.npz"
        if not circuit_path.exists():
            axes, angles, orders = draw_circuit(w, depth, r)
            npz_write(circuit_path, axes=axes, angles=angles, orders=orders)
        stem = f"{method}_w{w}_r{r:02d}_d{depth}"
        receipt = OUT/"receipts"/(stem+".json")
        raw = OUT/"raw"/(stem+".csv")
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            if saved["status"] == "complete" and saved["raw_sha256"] == sha(raw):
                continue
        pending.append((w, r, None, method))
    start = time.time()
    print(f"campaign {len(jobs)} cases, {len(pending)} pending, {workers} workers", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run_sample, *job) for job in pending]
        for future in as_completed(futures):
            future.result()
    json_write(OUT/"campaign_run.json", dict(status="complete", expected_cases=len(jobs),
                workers=workers, replicas=replicas, exact_max_width=exact_max_width,
                seconds=time.time()-start, source_sha256=source_hash, host=socket.gethostname(),
                pbs_jobid=os.environ.get("PBS_JOBID")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["validate", "ensemble", "sample", "conditional", "benchmark", "campaign"])
    parser.add_argument("--width", type=int, default=4)
    parser.add_argument("--replica", type=int, default=0)
    parser.add_argument("--max-depth", type=int)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--replicas", type=int, default=20)
    parser.add_argument("--exact-max-width", type=int, default=6)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    source_path = OUT/"sources"/(sha(__file__)+".py")
    source_path.parent.mkdir(parents=True, exist_ok=True)
    if not source_path.exists():
        source_path.write_bytes(Path(__file__).read_bytes())
    if args.command == "validate":
        validate()
    elif args.command == "ensemble":
        run_ensemble()
    elif args.command == "benchmark":
        benchmark()
    elif args.command == "campaign":
        campaign(args.workers, args.replicas, args.exact_max_width)
    else:
        run_sample(args.width, args.replica, args.max_depth,
                   "exact_ptm" if args.command == "sample" else "exact_angle_mean")


if __name__ == "__main__":
    main()
