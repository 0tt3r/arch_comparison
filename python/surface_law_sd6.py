"""Surface-code law under standard uniform circuit noise (SD6), results/surface_law_sd6_v1.

The first calibration (surface_law_calibration.py, results/surface_law_v1) used a CZ schedule
with one-qubit noise p/10. This version uses the standard depolarizing model of the
cultivation and factory literature: Stim's rotated memory circuit with native CNOTs and
depolarizing error p after every one- and two-qubit gate, p on every idle qubit in every
layer (including the measurement and reset layers), and flip probability p after resets
and before measurements. Sampling, decoding (PyMatching), fitting and the holdout reuse
surface_law_calibration.py unchanged; only the circuit and the output folder differ.

    python surface_law_sd6.py validate
    python surface_law_sd6.py run --workers 40        # training and holdout cells, resumable
    python surface_law_sd6.py fit
"""
import argparse
import json
import multiprocessing as mp
from types import SimpleNamespace

import stim

import surface_law_calibration as base

OUT = base.ROOT/'results/surface_law_sd6_v1'
base.OUT = OUT

ONE_QUBIT = {'H', 'S', 'S_DAG', 'X', 'Y', 'Z', 'SQRT_X', 'SQRT_X_DAG'}
TWO_QUBIT = {'CX', 'CZ', 'CY'}


def circuit(distance, rounds, p, basis):
    """SD6: noiseless Stim rotated memory circuit, then uniform noise p layer by layer."""
    template = stim.Circuit.generated('surface_code:rotated_memory_' + basis,
                                      distance=distance, rounds=rounds)
    qubits = sorted(template.get_final_qubit_coordinates())
    out = stim.Circuit()
    busy = set()

    def noise(gate, targets):
        if p and targets:
            out.append(gate, sorted(targets), p)

    def close_layer():
        noise('DEPOLARIZE1', [q for q in qubits if q not in busy])
        busy.clear()

    for ins in template.flattened():
        name = ins.name
        ts = [t.value for t in ins.targets_copy() if t.is_qubit_target]
        if name == 'TICK':
            close_layer()
            out.append('TICK')
        elif name in ('R', 'RX'):
            out.append(ins)
            noise('X_ERROR' if name == 'R' else 'Z_ERROR', ts)
            busy.update(ts)
        elif name in ('M', 'MX'):
            noise('X_ERROR' if name == 'M' else 'Z_ERROR', ts)
            out.append(ins)
            busy.update(ts)
        elif name == 'MR':
            noise('X_ERROR', ts)
            out.append(ins)
            noise('X_ERROR', ts)
            busy.update(ts)
        elif name in ONE_QUBIT:
            out.append(ins)
            noise('DEPOLARIZE1', ts)
            busy.update(ts)
        elif name in TWO_QUBIT:
            out.append(ins)
            noise('DEPOLARIZE2', ts)
            busy.update(ts)
        else:  # DETECTOR, OBSERVABLE_INCLUDE, QUBIT_COORDS, SHIFT_COORDS
            out.append(ins)
    return out


base.circuit = circuit
base.write_report = lambda result: None     # the v1 report text describes the p/10 model


def cases():
    """Same design as results/surface_law_v1 (training and holdout), fresh seeds."""
    old = json.loads((base.ROOT/'results/surface_law_v1/manifest.json').read_text())['cases']
    hold = json.loads((base.ROOT/'results/surface_law_v1/holdout_manifest.json').read_text())['cases']
    out = {}
    for group, rows, offset in (('cases', old, 0), ('holdout', hold, 500)):
        out[group] = []
        for i, c in enumerate(rows):
            seed = 202609260000 + offset + i
            key = f"d{c['distance']}_{c['basis']}_p{c['p']:g}_r{c['rounds']}_s{seed}"
            out[group].append({**{k: c[k] for k in ('p', 'distance', 'shots', 'rounds', 'basis', 'batch')},
                               'seed': seed, 'key': key, 'worker': 0})
    return out


def validate():
    """Noiseless circuits never fire a detector; every noisy DEM decomposes for matching."""
    rows = []
    for d in (3, 5):
        for basis in ('x', 'z'):
            clean = circuit(d, 3*d, 0., basis)
            det = clean.compile_detector_sampler(seed=1).sample(256)
            noisy = circuit(d, 3*d, 1e-3, basis)
            dem = noisy.detector_error_model(decompose_errors=True)
            rows.append(dict(distance=d, basis=basis, noiseless_detection_events=int(det.sum()),
                             detectors=noisy.num_detectors, error_mechanisms=dem.num_errors))
            assert det.sum() == 0
    result = dict(status='pass', model='SD6: p after every 1q/2q gate, p per idle qubit per layer, '
                  'flip p after reset and before measurement', rows=rows)
    base.atomic(OUT/'validation.json', result)
    print(json.dumps(result, indent=1))


def _sample(job):
    group, c = job
    base.sample(SimpleNamespace(group=group, **{k: c[k] for k in ('p', 'distance', 'shots', 'rounds',
                                                                    'basis', 'seed', 'batch')}))
    return c['key']


def run(workers):
    plan = cases()
    manifest = dict(created='2026-09-24', model='SD6 uniform circuit noise', cases=plan['cases'],
                    worker_sha256=base.sha(base.ROOT/'python/surface_law_calibration.py'),
                    sd6_module_sha256=base.sha(__file__))
    OUT.mkdir(parents=True, exist_ok=True)
    base.atomic(OUT/'manifest.json', manifest)
    base.atomic(OUT/'holdout_manifest.json', {**manifest, 'cases': plan['holdout']})
    jobs = [('cases', c) for c in plan['cases']] + [('holdout', c) for c in plan['holdout']]
    jobs = [j for j in jobs if not (OUT/j[0]/(j[1]['key']+'.json')).exists()]
    jobs.sort(key=lambda j: -j[1]['shots']*j[1]['rounds']*j[1]['distance']**2)
    for lock in OUT.rglob('*.lock'):
        lock.unlink()            # stale locks from an interrupted job; completed cells are kept
    with mp.get_context('fork').Pool(workers) as pool:
        for key in pool.imap_unordered(_sample, jobs):
            print('done', key, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('validate', 'run', 'fit'))
    parser.add_argument('--workers', type=int, default=40)
    args = parser.parse_args()
    if args.command == 'validate':
        validate()
    elif args.command == 'run':
        run(args.workers)
    else:
        base.fit()


if __name__ == '__main__':
    main()
