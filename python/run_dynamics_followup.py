"""Direct checks of frozen allocations and a bounded, restartable time sweep."""
import argparse
import json
import os

from budget_optimization import algorithm_for
from large_dynamics import ROOT, atomic_json, solve
from small_dynamics import rectangle


def settings():
    return json.loads((ROOT/'configs/dynamics_followup.json').read_text())


def headline_rows(config=None):
    config = settings() if config is None else config
    rows = json.loads((ROOT/config['selected_input']).read_text())
    return [v for v in rows if v['layout'] in config['headline_layouts'][v['arm']]]


def selected_cases(rows):
    return sorted({(algorithm_for(v['arm']), v['r']) for v in rows})


def tail_spread(records, reference, site):
    a, b = [abs(v['result']['doublons'][site]-reference)*v['contract']['r']**2
            for v in records[-2:]]
    return abs(a-b)/max(a, b, 1e-10)


def run(task):
    if not os.getenv('PBS_JOBID'):
        raise RuntimeError('Run full 4x4 cases on a PBS compute node')
    config = settings(); out = ROOT/config['output']
    out.mkdir(parents=True, exist_ok=True)
    threads = config['threads']
    site = rectangle(4, 4).index(tuple(config['site']))
    if task == 'selected':
        cases = []
        for algorithm, r in selected_cases(headline_rows(config)):
            name = f'4x4_tau1_{algorithm}_r{r}_t{threads}.json'
            prior = ROOT/'results/dynamics_4x4/tau1'/name
            # Reuse only completed existing records; never add to the old study.
            folder = prior.parent if prior.exists() else out/'tau1'
            solve(4, 4, 1., algorithm, r, folder, threads)
            cases.append(str((folder/name).relative_to(ROOT)))
        atomic_json(out/'selected_completion.json', {'task': task, 'cases': cases,
            'row_count': len(headline_rows(config)), 'status': 'complete'})
        return
    tau = float(task.removeprefix('tau'))
    if tau not in config['times']:
        raise ValueError(f'Unconfigured time: {tau}')
    folder = out/f'tau{tau:g}'
    exact = solve(4, 4, tau, 'exact', 0, folder, threads)
    reference = exact['result']['doublons'][site]
    convergence = []
    for algorithm in ('plaquette', 'termwise'):
        records = [solve(4, 4, tau, algorithm, int(k*tau), folder, threads)
                   for k in config['initial_steps_per_tau']]
        spread = tail_spread(records, reference, site)
        while spread > config['maximum_relative_tail_spread']:
            r = 2*records[-1]['contract']['r']
            if r > config['maximum_steps_per_tau']*tau:
                raise RuntimeError(f'{algorithm}, tau={tau}: coefficient has not converged')
            records.append(solve(4, 4, tau, algorithm, r, folder, threads))
            spread = tail_spread(records, reference, site)
        convergence.append({'algorithm': algorithm, 'steps': [v['contract']['r'] for v in records],
                            'relative_tail_coefficient_spread': spread})
    atomic_json(out/f'{task}_completion.json', {'task': task, 'tau': tau,
        'rows': convergence, 'status': 'complete', 'site': config['site'],
        'scope': 'Tail convergence on this finite OBC lattice; no bulk error bound'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task', choices=['selected', 'tau2', 'tau4'])
    run(parser.parse_args().task)
