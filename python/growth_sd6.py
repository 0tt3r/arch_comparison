"""STAR patch-growth bounds under standard uniform circuit noise (SD6), results/growth_sd6_v1.

The growth campaign of results/frontier_v1 (star_growth_probe.growth_circuit) used one-qubit
and idle noise p/10 plus a 0.8p per-round reset allowance. SD6 uses p after every gate, p
per idle qubit per layer, and flips p on reset and measurement, with no time-based allowance.
The circuit here is growth_circuit with exactly those three noise edits, applied to its source
and asserted, so the schedule, detectors and observable are unchanged. The 24 cells, shot counts
and the summary format (pairs with pauli_total_upper) repeat frontier_growth.py and
exposure_diagnostics.growth_summary.

    python growth_sd6.py --workers 24          # on a PBS node; resumable
"""
import argparse
import hashlib
import inspect
import json
import os
import resource
import socket
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pymatching
import stim
from scipy.stats import beta

import star_growth_probe as probe
from model import log_gamma2_dephasing

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/growth_sd6_v1'


def _sd6_growth_circuit():
    source = inspect.getsource(probe.growth_circuit)
    edits = [("elif name=='H':circuit.append('DEPOLARIZE1',targets,p/10)",
              "elif name=='H':circuit.append('DEPOLARIZE1',targets,p)"),
             ("if idle:circuit.append('DEPOLARIZE1',idle,p/10)",
              "if idle:circuit.append('DEPOLARIZE1',idle,p)"),
             ("if p:circuit.append('DEPOLARIZE1',sorted(coordinates),.8*p)",
              "pass  # SD6: no time-based reset allowance")]
    for old, new in edits:
        assert source.count(old) == 1, old
        source = source.replace(old, new)
    namespace = dict(vars(probe))
    exec(compile(source.replace('def growth_circuit(', 'def growth_circuit_sd6(', 1), 'growth_sd6', 'exec'), namespace)
    return namespace['growth_circuit_sd6']


growth_circuit = _sd6_growth_circuit()


def cases():
    c = json.loads((ROOT/'configs/frontier_study.json').read_text())['growth']
    rows = []
    for di, df in c['pairs']:
        for p in c['p_values']:
            for basis in c['bases']:
                key = f'd{di}_d{df}_{basis}_p{p:g}'
                seed = c['seed'] + 7 + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
                rows.append(dict(key=key, di=di, df=df, p=p, basis=basis, shots=c['shots_per_case'], seed=seed))
    return rows


def execute(row):
    folder = OUT/'cases'
    folder.mkdir(parents=True, exist_ok=True)
    path = folder/(row['key']+'.json')
    if path.exists() and json.loads(path.read_text())['status'] == 'complete':
        return json.loads(path.read_text())
    start = time.monotonic()
    circuit = growth_circuit(row['di'], row['df'], row['p'], row['basis'])
    ideal = growth_circuit(row['di'], row['df'], 0, row['basis'])
    det, obs = ideal.compile_detector_sampler(seed=7).sample(64, separate_observables=True)
    assert not det.any() and not obs.any()
    decoder = pymatching.Matching.from_detector_error_model(circuit.detector_error_model(decompose_errors=True))
    sampler = circuit.compile_detector_sampler(seed=row['seed'])
    shots = failures = 0
    batch = max(128, min(5000, 20000000//circuit.num_detectors))
    while shots < row['shots']:
        n = min(batch, row['shots']-shots)
        det, obs = sampler.sample(n, separate_observables=True)
        failures += int(np.any(decoder.decode_batch(det) != obs, axis=1).sum())
        shots += n
    cells = len(cases())
    upper = float(beta.ppf(1-.05/cells, failures+1, shots-failures))
    text = str(circuit)+'\n'
    result = dict(status='complete', contract=row, shots=shots, failures=failures, rate=failures/shots,
                  simultaneous_upper95=upper, circuit_sha256=hashlib.sha256(text.encode()).hexdigest(),
                  detectors=circuit.num_detectors, elapsed_seconds=time.monotonic()-start,
                  peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, host=socket.gethostname(),
                  pbs_job=os.environ.get('PBS_JOBID'), stim=stim.__version__, pymatching=pymatching.__version__,
                  noise='SD6 uniform circuit noise',
                  scope='Ideal encoded input; decoded growth + target memory + readout. Not a full injected-state channel.')
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(result, indent=2)+'\n')
    tmp.replace(path)
    print(row['key'], failures, shots, round(result['elapsed_seconds']), flush=True)
    return result


def summary():
    rows = [json.loads((OUT/'cases'/(r['key']+'.json')).read_text()) for r in cases()]
    assert all(v['status'] == 'complete' for v in rows)
    pairs = []
    for di, df in sorted({(r['contract']['di'], r['contract']['df']) for r in rows}):
        for p in (.001, .0001):
            selected = [v for v in rows if (v['contract']['di'], v['contract']['df'], v['contract']['p']) == (di, df, p)]
            assert len(selected) == 2
            q = min(.5, sum(v['simultaneous_upper95'] for v in selected))
            pairs.append(dict(di=di, df=df, p=p, failures={v['contract']['basis']: v['failures'] for v in selected},
                              pauli_total_upper=q, chi_per_growth_upper=log_gamma2_dephasing(q),
                              scope='Union of X/Z upper bounds under SD6 noise; ideal encoded input'))
    (OUT/'growth_summary.json').write_text(json.dumps(dict(cases=rows, pairs=pairs,
        shots=sum(v['shots'] for v in rows), noise='SD6 uniform circuit noise',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2)+'\n')
    print('pairs', len(pairs))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=24)
    parser.add_argument('--summary', action='store_true')
    args = parser.parse_args()
    if not args.summary:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(execute, cases()))
    summary()


if __name__ == '__main__':
    main()
