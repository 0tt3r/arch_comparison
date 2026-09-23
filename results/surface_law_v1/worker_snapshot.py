"""Single-process C6/C7 surface-memory calibration; all output is isolated.

Commands: validate, sample, fit, audit. See results/surface_law_v1/manifest.json.
No multiprocessing, sinter worker pool, or implicit BLAS threading is used.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import socket
import sys
import time

import numpy as np
import pymatching
import scipy
from scipy.optimize import minimize, brentq
from scipy.special import xlogy
from scipy.stats import beta, chi2
import stim

ROOT = Path('/lcrc/project/QC_NA_Sim/architecture_comparison')
OUT = ROOT / 'results/surface_law_v1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def circuit(distance, rounds, p, basis):
    """Compile Stim's hook-safe rotated schedule into dressed CZ macro-layers.

    Every CNOT is H(target)-CZ-H(target), each H has p/10 noise.
    The generator's ancilla H gates are also noisy. Basis changes are fused
    into interaction slots: no additional idle macro-layer or time is assigned.
    Four interaction slots, M slot, R slot; idle p per unused slot. Initial
    product prep and final data measurement carry flip p. MR is split M/R.
    """
    template = stim.Circuit.generated('surface_code:rotated_memory_' + basis,
                                      distance=distance, rounds=rounds)
    coordinates = template.get_final_qubit_coordinates()
    all_qubits = set(coordinates)
    data = {q for q, (x, y) in coordinates.items() if x % 2 and y % 2}
    out = stim.Circuit()

    def noise(gate, targets, probability):
        if p and targets:
            out.append(gate, targets, probability)

    for ins in template.flattened():
        name = ins.name
        ts = [t.value for t in ins.targets_copy() if t.is_qubit_target]
        if name == 'TICK':
            continue  # Ticks below mark precisely the six declared macro-slots.
        if name == 'CX':
            targets = ts[1::2]
            out.append('H', targets)
            noise('DEPOLARIZE1', targets, p / 10)
            out.append('CZ', ts)
            noise('DEPOLARIZE2', ts, p)
            out.append('H', targets)
            noise('DEPOLARIZE1', targets, p / 10)
            noise('DEPOLARIZE1', sorted(all_qubits - set(ts)), p)
            out.append('TICK')
        elif name == 'H':
            out.append(ins)
            noise('DEPOLARIZE1', ts, p / 10)
        elif name in ('R', 'RX'):
            out.append(ins)
            noise('Z_ERROR' if name == 'RX' else 'X_ERROR', ts, p)
        elif name == 'MR':
            noise('X_ERROR', ts, p)
            out.append('M', ts)
            noise('DEPOLARIZE1', sorted(data), p)
            out.append('TICK')
            out.append('R', ts)
            noise('X_ERROR', ts, p)
            noise('DEPOLARIZE1', sorted(data), p)
            out.append('TICK')
        elif name in ('M', 'MX'):
            noise('Z_ERROR' if name == 'MX' else 'X_ERROR', ts, p)
            out.append(ins)
        else:
            out.append(ins)
    return out


def interval(k, n, alpha=.05):
    return [0. if k == 0 else float(beta.ppf(alpha/2, k, n-k+1)),
            1. if k == n else float(beta.ppf(1-alpha/2, k+1, n-k))]


def per_round(q, rounds):
    return float(-.5 * np.expm1(np.log1p(-2 * min(q, .499999999)) / rounds))


def validate():
    rows = []
    for d in (3, 5, 7, 9):
        for basis in ('x', 'z'):
            c = circuit(d, 2*d, .001, basis)
            ideal = circuit(d, 2*d, 0, basis)
            det, obs = ideal.compile_detector_sampler(seed=17).sample(256, separate_observables=True)
            assert not det.any() and not obs.any()
            dem = c.detector_error_model(decompose_errors=True)
            decoder = pymatching.Matching.from_detector_error_model(dem)
            coords = c.get_final_qubit_coordinates()
            data = [q for q, (x, y) in coords.items() if x % 2 and y % 2]
            logical = [q for q in data if coords[q][1 if basis == 'x' else 0] == 1]
            control = stim.Circuit()
            inserted = False
            for ins in ideal:
                if ins.name == ('MX' if basis == 'x' else 'M') and set(t.value for t in ins.targets_copy()) == set(data):
                    control.append('Z_ERROR' if basis == 'x' else 'X_ERROR', logical, 1)
                    inserted = True
                control.append(ins)
            det, obs = control.compile_detector_sampler(seed=18).sample(256, separate_observables=True)
            assert inserted and not det.any() and obs.all()
            assert np.all(decoder.decode_batch(det) != obs)
            # Both representations of batch decoding must agree bit-for-bit.
            det, obs = c.compile_detector_sampler(seed=19).sample(2048, separate_observables=True)
            packed = np.packbits(det, axis=1, bitorder='little')
            assert np.array_equal(decoder.decode_batch(det), decoder.decode_batch(packed, bit_packed_shots=True))
            assert sum(ins.name == 'TICK' for ins in c) == 12*d
            assert sum(len(ins.targets_copy()) // 2 for ins in c if ins.name == 'CZ') == 2*d * (4*d*(d-1))
            distance_found = len(c.shortest_graphlike_error())
            assert distance_found == d, (d, basis, distance_found)
            rows.append(dict(distance=d, basis=basis, graphlike_distance=distance_found,
                             noiseless_shots=256, logical_flip_shots=256,
                             packed_decode_equal=True, macro_layers_per_round=6))
    result = dict(status='pass', rows=rows, stim=stim.__version__, pymatching=pymatching.__version__)
    atomic(OUT/'validation.json', result)
    print(json.dumps(result), flush=True)


def sample(args):
    folder = OUT / args.group
    folder.mkdir(parents=True, exist_ok=True)
    key = f'd{args.distance}_{args.basis}_p{args.p:g}_r{args.rounds}_s{args.seed}'
    path = folder / (key + '.json')
    contract = dict(distance=args.distance, basis=args.basis, p=args.p, rounds=args.rounds,
                    target_shots=args.shots, seed=args.seed, batch=args.batch,
                    program_sha256=sha(__file__), decoder='pymatching MWPM; correlations disabled')
    if path.exists():
        previous = json.loads(path.read_text())
        assert previous['contract'] == contract and previous['status'] == 'complete'
        print('Already complete:', key, flush=True)
        return
    lock = folder / (key + '.lock')
    with lock.open('x') as f:
        json.dump(dict(pid=os.getpid(), host=socket.gethostname(), job=os.environ.get('PBS_JOBID')), f)
    started = time.monotonic()
    try:
        c = circuit(args.distance, args.rounds, args.p, args.basis)
        cpath = folder / (key + '.stim')
        cpath.write_text(str(c) + '\n')
        decoder = pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
        sampler = c.compile_detector_sampler(seed=args.seed)
        shots = failures = 0
        batch_records = []
        while shots < args.shots:
            n = min(args.batch, args.shots-shots)
            det, obs = sampler.sample(n, separate_observables=True, bit_packed=True)
            predicted = decoder.decode_batch(det, bit_packed_shots=True, bit_packed_predictions=True)
            k = int(np.count_nonzero(np.any(predicted != obs, axis=1)))
            failures += k
            shots += n
            batch_records.append([n, k])
            elapsed = time.monotonic()-started
            if shots == args.shots or len(batch_records) % 10 == 0:
                result = dict(status='complete' if shots == args.shots else 'running', contract=contract,
                    shots=shots, failures=failures, interval95=interval(failures, shots),
                    finite_window_per_round=per_round(failures/shots, args.rounds),
                    finite_window_interval95=[per_round(q,args.rounds) for q in interval(failures,shots)],
                    batches=batch_records, elapsed_seconds=elapsed, shots_per_second=shots/elapsed,
                    circuit_sha256=sha(cpath), host=socket.gethostname(), pbs_job=os.environ.get('PBS_JOBID'),
                    stim=stim.__version__, pymatching=pymatching.__version__, numpy=np.__version__,
                    scipy=scipy.__version__, python=sys.version, machine=platform.machine(),
                    peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                atomic(path if shots == args.shots else path.with_suffix('.progress.json'), result)
                print(key, shots, failures, round(elapsed,2), flush=True)
        progress = path.with_suffix('.progress.json')
        if progress.exists():
            progress.unlink()
    finally:
        lock.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['validate', 'sample', 'fit', 'audit'])
    parser.add_argument('--distance', type=int, default=3)
    parser.add_argument('--basis', choices=['x','z'], default='x')
    parser.add_argument('--p', type=float, default=.001)
    parser.add_argument('--rounds', type=int, default=30)
    parser.add_argument('--shots', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=20260923)
    parser.add_argument('--batch', type=int, default=4096)
    parser.add_argument('--group', default='pilot')
    args = parser.parse_args()
    if args.command == 'validate':
        validate()
    elif args.command == 'sample':
        sample(args)
    elif args.command == 'fit':
        fit()
    else:
        audit()


if __name__ == '__main__':
    main()
