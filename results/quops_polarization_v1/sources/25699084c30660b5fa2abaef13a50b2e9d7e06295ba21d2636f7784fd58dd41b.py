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


def exact_two_layers(w, axes, angles, orders, p, idle=True):
    """Exact individual-circuit prefixes D<=2 without constructing a full PTM.

    M starts diagonal and remains diagonal after the first layer. The trace
    after the second layer needs only diag(R M R^T). This shortcut MUST NOT
    be extended to a third layer: discarded off-diagonals can then return.
    """
    assert len(axes) <= 2
    v = np.ones(4**w)
    weight = geometry(w)[2]
    curve = []
    for ax, theta, order in zip(axes, angles, orders):
        for q in range(w):
            stride = 4**q
            a = v.reshape(-1, 4, stride)
            b, c = int(ax[q] % 3)+1, int((ax[q]+1) % 3)+1
            co, si = math.cos(2*theta[q])**2, math.sin(2*theta[q])**2
            x, y = a[:, b, :].copy(), a[:, c, :].copy()
            a[:, b, :], a[:, c, :] = co*x+si*y, si*x+co*y
        perm, _, active, waiting = matching_map(w, order)
        factor = (1-4*p/30)**weight[perm]*(1-16*p/15)**active
        if idle:
            factor *= (1-4*p/30)**waiting
        v = v[perm]*factor
        curve.append((v.sum()-1)/(len(v)-1))
    return np.array(curve)


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
    for w in range(2, 6):
        ax, theta, order = draw_circuit(w, 2, 913)
        a = exact_circuit(w, ax, theta, order, .13)
        b = exact_two_layers(w, ax, theta, order, .13)
        err = float(np.max(abs(a-b)))
        assert err < 1e-13
        checks.append(dict(check="two_layer_shortcut_vs_full_PTM", w=w, max_abs=err))
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
    if not 2 <= exact_max_width <= 8 or workers < 1 or replicas < 1:
        raise ValueError("Require width 2..8 and positive workers/replicas")
    physical_memory = os.sysconf("SC_PAGE_SIZE")*os.sysconf("SC_PHYS_PAGES")
    if workers*2*8*16**exact_max_width > .7*physical_memory:
        raise ValueError("Dense PTM workers would exceed 70% of physical memory; reduce --workers")
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


def run_shallow():
    rows = []
    for w in (7, 8):
        for r in range(20):
            depth = max(depth_for(w, p) for p in PS)
            ax, theta, order = draw_circuit(w, depth, r)
            for p in PS:
                for idle in (True, False) if w % 2 else (True,):
                    curve = exact_two_layers(w, ax[:2], theta[:2], order[:2], p, idle)
                    for d, gamma in enumerate(curve, 1):
                        rows.append(dict(method="exact_two_layer_PTM", w=w, replica=r, p=p,
                                         idle=int(idle), D=d, gamma=float(gamma)))
    csv_write(OUT/"shallow_raw.csv", rows)


def analyze():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from collections import defaultdict
    import scipy.stats

    with (OUT/"ensemble.csv").open() as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k in r:
            if k != "method":
                r[k] = float(r[k])
        for k in ("w", "D", "idle"):
            r[k] = int(r[k])
    lookup = {(r["w"], r["p"], r["idle"], r["D"]): r for r in rows}
    grouped = defaultdict(list)
    for r in rows:
        grouped[r["w"], r["p"], r["idle"]].append(r)

    def fit_rate(selected):
        # Equal weight per (w,p,idle) cell, not per row. Dense low-p depth grids
        # must not receive ten times the weight of the high-p cells.
        count = defaultdict(int)
        for r in selected:
            count[r["w"], r["p"], r["idle"]] += 1
        xy = xx = 0.
        for r in selected:
            weight = 1/count[r["w"], r["p"], r["idle"]]
            xy += weight*r["exposure"]*(-math.log(r["gamma"]))
            xx += weight*r["exposure"]**2
        return xy/xx

    train = [r for r in rows if r["p"] in (0.001, 0.01) and r["idle"]]
    global_c = fit_rate(train)
    by_width = {w: fit_rate([r for r in train if r["w"] == w]) for w in range(2, 9)}
    for r in rows:
        r["fit_global"] = math.exp(-global_c*r["exposure"])
        r["fit_width"] = math.exp(-by_width[r["w"]]*r["exposure"])
        r["ratio_fit_global"] = r["gamma"]/r["fit_global"]
        r["ratio_fit_width"] = r["gamma"]/r["fit_width"]
    csv_write(OUT/"comparison.csv", rows)

    samples = defaultdict(list)
    receipts = []
    for path in sorted((OUT/"receipts").glob("*.json")):
        receipt = json.loads(path.read_text())
        if receipt["replica"] >= 20:
            continue
        raw = OUT/"raw"/(path.stem+".csv")
        assert sha(raw) == receipt["raw_sha256"]
        circuit = OUT/"circuits"/f'w{receipt["w"]}_r{receipt["replica"]:02d}_d{receipt["max_depth"]}.npz'
        assert sha(circuit) == receipt["circuit_sha256"]
        receipts.append(receipt)
        with raw.open() as f:
            for record in csv.DictReader(f):
                key = (record["method"], int(record["w"]), float(record["p"]),
                       int(record["idle"]), int(record["D"]))
                samples[key].append(float(record["gamma"]))
    with (OUT/"shallow_raw.csv").open() as f:
        for record in csv.DictReader(f):
            key = (record["method"], int(record["w"]), float(record["p"]),
                   int(record["idle"]), int(record["D"]))
            samples[key].append(float(record["gamma"]))
    summary = []
    for (method, w, p, idle, depth), values in sorted(samples.items()):
        if (w, p, idle, depth) not in lookup:
            continue
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=1)) if len(values) > 1 else 0.
        se = sd/math.sqrt(len(values))
        ensemble = lookup[w, p, idle, depth]["gamma"]
        half = float(scipy.stats.t.ppf(.975, len(values)-1))*se
        summary.append(dict(method=method, w=w, p=p, idle=idle, D=depth, n=len(values),
                            mean=mean, sd=sd, se=se, ci95_low=mean-half, ci95_high=mean+half,
                            exact_ensemble=ensemble, mean_minus_ensemble=mean-ensemble,
                            ratio_mean_to_ensemble=mean/ensemble))
    csv_write(OUT/"sample_summary.csv", summary)
    endpoints = []
    for w in range(2, 9):
        for p in PS:
            r = grouped[w, p, 1][-1]
            off = lookup[w, p, 0, r["D"]]["gamma"]
            endpoints.append(dict(w=w, p=p, D=r["D"], gamma=r["gamma"],
                                  ratio_local=r["ratio_local"], ratio_embedded=r["ratio_embedded"],
                                  idle_on_over_off=r["gamma"]/off))
    csv_write(OUT/"endpoints.csv", endpoints)
    held = [r for r in rows if r["p"] == .003 and r["idle"]]
    def max_error(key, selected):
        return max(abs(r[key]-1) for r in selected)
    fit = dict(model="gamma=exp(-c*exposure); exposure=p*(n2+(n1+n_idle_10ns)/10)",
               training_p=[.001, .01], heldout_p=.003, fit_idle=True,
               weighting="equal total weight per width/p cell; zero intercept; squared log residuals",
               global_c=global_c, c_by_width=by_width,
               heldout_max_relative_error_global=max_error("ratio_fit_global", held),
               heldout_max_relative_error_width=max_error("ratio_fit_width", held),
               local_product_max_relative_ratio=max(r["ratio_local"] for r in rows),
               embedded_product_max_relative_error=max_error("ratio_embedded", rows),
               sample_mean_max_relative_deviation=max(abs(r["ratio_mean_to_ensemble"]-1) for r in summary if r["method"] == "exact_ptm"),
               scope="Exact all-ensemble means w=2..8; 20 individual circuits w=2..6 at all depths; w=7,8 individual checks only D=1,2. Full-depth individual w=7,8 requirement not completed.")
    json_write(OUT/"fit.json", fit)
    counts = {method: {str(w): sum(r["method"] == method and r["w"] == w for r in receipts)
                       for w in range(2, 9)} for method in ("exact_ptm", "exact_angle_mean")}
    assert all(counts["exact_ptm"][str(w)] == 20 for w in range(2, 7)), counts
    assert all(counts["exact_angle_mean"][str(w)] == 20 for w in range(2, 9)), counts
    assert all(r["n"] == 20 for r in summary), "Incomplete sample cells"
    for method, widths in (("exact_ptm", range(2, 7)), ("exact_angle_mean", range(2, 9))):
        for w in widths:
            assert sorted(r["replica"] for r in receipts if r["w"] == w and r["method"] == method) == list(range(20))
    for r in receipts:
        expected_rows = sum(depth_for(r["w"], p) for p in PS)*(2 if r["w"] % 2 else 1)
        assert r["rows"] == expected_rows
        assert (OUT/"sources"/(r["source_sha256"]+".py")).exists()
    audit = dict(status="complete_for_reported_scope", original_brief_status="incomplete_individual_deep_w7_w8",
                 receipts=len(receipts), counts=counts, ensemble_cells=len(grouped),
                 all_integer_depths=True, target=TARGET, sample_cells=len(summary),
                 input_sources={str(path.relative_to(ROOT)): sha(path) for path in
                                (ROOT/"PLAN.md", ROOT/"refs/2609.12146v1.pdf", ROOT/"python/quops_study.py")})
    json_write(OUT/"audit.json", audit)

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.0), constrained_layout=True)
    colors = {w: plt.get_cmap("viridis")((w-2)/6) for w in range(2, 9)}
    for w in range(2, 9):
        cell = grouped[w, .003, 1]
        x = np.array([r["exposure"] for r in cell])
        y = np.array([r["gamma"] for r in cell])
        axes[0].plot(x, y, color=colors[w], label=f"w={w}")
        axes[1].plot(x, [r["ratio_local"] for r in cell], color=colors[w])
        axes[1].plot(x, [r["ratio_embedded"] for r in cell], color=colors[w], linestyle=":", linewidth=1)
        idle_ratio = [r["gamma"]/lookup[w, .003, 0, r["D"]]["gamma"] for r in cell]
        axes[2].plot(x, idle_ratio, color=colors[w])
        if w <= 6:
            sm = [r for r in summary if r["method"] == "exact_ptm" and r["w"] == w and r["p"] == .003 and r["idle"]]
            for r in sm[::max(1, len(sm)//7)]:
                ex = lookup[w, .003, 1, r["D"]]["exposure"]
                axes[0].errorbar(ex, r["mean"], yerr=2.093*r["se"], color=colors[w], fmt="o", ms=2.5, capsize=1)
    axes[0].set(ylabel="Mean process polarization", title="Exact mean; dots: 20 exact circuits")
    axes[0].legend(ncol=2, fontsize=7, frameon=False)
    axes[1].set(ylabel="Exact mean / gate product", title="Solid: local; dotted: embedded")
    axes[2].set(ylabel="Mean with idle / without idle", title="One 10 ns idle per odd-width layer")
    for ax in axes:
        ax.set_xlabel(r"Exposure $x=p[n_2+(n_1+n_{\rm idle})/10]$")
        ax.grid(alpha=.18)
    fig.savefig(OUT/"polarization.png", dpi=180)
    fig.savefig(OUT/"polarization.pdf")
    plt.close(fig)
    print(json.dumps(fit, indent=2), flush=True)


def replay():
    checks = []
    for w in range(2, 9):
        for r in range(20):
            depth = max(depth_for(w, p) for p in PS)
            path = OUT/"circuits"/f"w{w}_r{r:02d}_d{depth}.npz"
            ax, theta, order = draw_circuit(w, depth, r)
            with np.load(path) as data:
                assert np.array_equal(data["axes"], ax)
                assert np.array_equal(data["angles"], theta)
                assert np.array_equal(data["orders"], order)
    for w in (2, 3, 6):
        depth = max(depth_for(w, p) for p in PS)
        with np.load(OUT/"circuits"/f"w{w}_r00_d{depth}.npz") as data:
            stop = depth_for(w, .01)
            curve = exact_circuit(w, data["axes"][:stop], data["angles"][:stop], data["orders"][:stop], .01)
        path = OUT/"raw"/f"exact_ptm_w{w}_r00_d{depth}.csv"
        with path.open() as f:
            old = np.array([float(r["gamma"]) for r in csv.DictReader(f)
                            if float(r["p"]) == .01 and int(r["idle"]) == 1])
        error = float(np.max(abs(curve-old)))
        assert error < 1e-12
        checks.append(dict(w=w, p=.01, replica=0, max_abs=error, bitwise_equal=bool(np.array_equal(curve, old))))
    json_write(OUT/"replay.json", dict(status="passed", all_140_circuit_files_seed_reproduced=True,
                                      checks=checks, source_sha256=sha(__file__)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["validate", "ensemble", "sample", "conditional", "benchmark", "campaign", "shallow", "analyze", "replay"])
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
    elif args.command == "shallow":
        run_shallow()
    elif args.command == "analyze":
        analyze()
    elif args.command == "replay":
        replay()
    else:
        run_sample(args.width, args.replica, args.max_depth,
                   "exact_ptm" if args.command == "sample" else "exact_angle_mean")


if __name__ == "__main__":
    main()
