"""Audit saved step allocations, then isolate finite-patch capacity sensitivity."""
import argparse
import csv
import json
from math import ceil, exp, isclose, log, log10

from budget_optimization import algorithm_for, coefficient, optimize, settings as budget_settings
from delivery_sensitivity import choices, key, settings as delivery_settings
from discussion_model import assumptions, capacity
from large_dynamics import ROOT, sha
from run_dynamics_followup import headline_rows, settings
from small_dynamics import rectangle

OUT = ROOT/'results/dynamics_followup_v1'


def read(path):
    return json.loads(path.read_text())


def write(name, value):
    (OUT/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def measurements():
    paths = sorted((ROOT/'results/dynamics_4x4/tau1').glob('*.json'))
    paths += sorted((ROOT/'results/dynamics_4x4/rectangle_tau1').glob('*.json'))
    paths += sorted(OUT.glob('tau*/*.json'))
    records = [(p, read(p)) for p in paths]
    binary_hash = read(ROOT/'results/dynamics_4x4/build.json')['sha256']['build/hubbard_sector']
    seen = set(); rows = []
    for path, value in records:
        c = value['contract']; data = value['result']
        identity = (c['nx'], c['ny'], c['tau'], c['algorithm'], c['r'])
        assert identity not in seen, f'Duplicate scientific case: {identity}'
        seen.add(identity)
        assert value['status'] == 'complete'
        assert c['source_sha256'] == sha(ROOT/'cpp/hubbard_sector.cpp')
        assert c['driver_sha256'] == sha(ROOT/'python/large_dynamics.py')
        assert c['binary_sha256'] == binary_hash
        assert value['coefficient_sha256'] == sha(path.with_suffix('.coefficients.txt'))
        assert abs(data['norm_squared']-1) < 1e-8
        if c['algorithm'] == 'exact':
            assert value['chebyshev_state_tail_bound'] < 1e-12
            continue
        refs = [(p, v) for p, v in records if v['contract']['algorithm'] == 'exact'
                and all(v['contract'][k] == c[k] for k in ('nx', 'ny', 'tau'))]
        assert len(refs) == 1
        refpath, reference = refs[0]
        site = [2, 2] if c['nx'] == c['ny'] == 4 else [1, 1]
        i = rectangle(c['nx'], c['ny']).index(tuple(site))
        D = reference['result']['doublons'][i]
        error = data['doublons'][i]-D
        rows.append(dict(patch=f"{c['nx']}x{c['ny']}", tau=c['tau'], algorithm=c['algorithm'],
            r=c['r'], site=site, exact_doublon=D, doublon=data['doublons'][i],
            signed_error=error, coefficient=abs(error)*c['r']**2,
            source=str(path.relative_to(ROOT)), exact_source=str(refpath.relative_to(ROOT))))
    return rows


def saved_allocation_checks(rows):
    """Compare with the original reserve, before changing any allocation."""
    allowance = settings()['numerical_allowance']; checks = []
    for plan in headline_rows():
        match = [v for v in rows if v['patch'] == '4x4' and v['tau'] == 1
                 and v['algorithm'] == algorithm_for(plan['arm']) and v['r'] == plan['r']]
        assert len(match) == 1
        measured = match[0]; error = abs(measured['signed_error'])+allowance
        reserve = plan['absolute_error']*plan['budget_fractions']['trotter']
        checks.append({k: plan[k] for k in ('arm', 'layout', 'absolute_error', 'r', 'log10_Q_min')})
        checks[-1].update(original_trotter_reserve=reserve, measured_error=abs(measured['signed_error']),
            numerical_allowance=allowance, error_to_reserve=error/reserve,
            passed=error <= reserve, measurement_source=measured['source'])
    return checks


def completion_checks(rows):
    expected = {(tau, algorithm) for tau in (1, 2, 4) for algorithm in ('plaquette', 'termwise')}
    actual = {(v['tau'], v['algorithm']) for v in rows if v['patch'] == '4x4'}
    assert actual == expected, 'Missing time/ordering in the finite-lattice census'
    selected = read(OUT/'selected_completion.json')
    assert selected['status'] == 'complete' and selected['row_count'] == 16
    expected_sources = {v['measurement_source'] for v in saved_allocation_checks(rows)}
    assert len(selected['cases']) == len(set(selected['cases'])) == 9
    assert set(selected['cases']) == expected_sources
    config = settings()
    for tau in config['times']:
        completed = read(OUT/f'tau{tau}_completion.json')
        assert completed['status'] == 'complete' and completed['site'] == config['site']
        assert {v['algorithm'] for v in completed['rows']} == {'plaquette', 'termwise'}
        for v in completed['rows']:
            group = sorted([m for m in rows if m['patch'] == '4x4' and m['tau'] == tau
                            and m['algorithm'] == v['algorithm']], key=lambda m: m['r'])
            assert v['steps'] == [m['r'] for m in group]
            assert {int(k*tau) for k in config['initial_steps_per_tau']} <= set(v['steps'])
            assert max(v['steps']) <= config['maximum_steps_per_tau']*tau
            a, b = [m['coefficient'] for m in group[-2:]]
            spread = abs(a-b)/max(a, b, 1e-10)
            assert isclose(spread, v['relative_tail_coefficient_spread'], abs_tol=1e-12)
            assert spread <= config['maximum_relative_tail_spread']


def layout_changes(arm, layout):
    if arm == 'nisq':
        return {} if layout == 'effective' else {'nisq_model': layout}
    return dict(choices(arm, delivery_settings()))[layout]


def direct_costs(rows):
    options = budget_settings(); config = delivery_settings()
    base = {**assumptions(), 'factory_choices': config['factory_choices']}
    chosen = []; allowance = settings()['numerical_allowance']
    for plan in headline_rows():
        group = [v for v in rows if v['patch'] == '4x4' and v['tau'] == 1
                 and v['algorithm'] == algorithm_for(plan['arm'])]
        candidates = []
        for v in group:
            error = abs(v['signed_error'])+allowance
            row = optimize(4, 1, plan['arm'], {**base, **layout_changes(plan['arm'], plan['layout'])},
                plan['absolute_error'], {**options, 'explicit_step_candidates': [v['r']]},
                C_override=error*v['r']**2)
            if row['log10_Q_min'] is None:
                continue
            row.update(layout=plan['layout'], measured_trotter_error=abs(v['signed_error']),
                       numerical_allowance=allowance, measurement_source=v['source'],
                       coefficient_patch='4x4', coefficient_policy='direct measured error',
                       original_r=plan['r'], original_log10_Q_min=plan['log10_Q_min'])
            candidates.append(row)
        same_r = next(v for v in candidates if v['r'] == plan['r'])
        for menu, eligible in [('saved_r', [same_r]), ('all_measured_r', candidates)]:
            chosen.append({**min(eligible, key=key), 'step_menu': menu})
    return chosen


def patch_coefficients(rows):
    """Retain pointwise converged tails separately from planning envelopes."""
    result = []
    for patch, tau, algorithm in sorted({(v['patch'], v['tau'], v['algorithm']) for v in rows}):
        group = sorted([v for v in rows if (v['patch'], v['tau'], v['algorithm']) ==
                        (patch, tau, algorithm)], key=lambda v: v['r'])
        a, b = [v['coefficient'] for v in group[-2:]]
        spread = abs(a-b)/max(a, b, 1e-10)
        assert spread <= settings()['maximum_relative_tail_spread']
        result.append(dict(patch=patch, tau=tau, algorithm=algorithm, coefficient=max(a, b),
            finest_coefficient=b, tail_spread=spread, largest_measured_coefficient=max(v['coefficient'] for v in group),
            exact_doublon=group[0]['exact_doublon'], steps=[v['r'] for v in group]))
    return result


def extended_coefficient(tau, algorithm, rows, margin=1.1):
    if tau not in read(ROOT/'configs/patch_sensitivity.json')['times']:
        raise ValueError('No extended-patch calibration at this time')
    # All finite-r measurements enter, not just the leading-order tail.
    extra = [v['coefficient'] for v in rows if v['tau'] <= tau and v['algorithm'] == algorithm]
    rectangle_rows = read(ROOT/'results/discussion_v3/rectangle_checks/calibration.json')['rows']
    extra.extend(abs(m['coefficient']) for v in rectangle_rows
                 if v['tau'] <= tau and v['algorithm'] == algorithm for m in v['measurements'])
    return max(coefficient(tau, algorithm), margin*max(extra))


def sensitivity(rows):
    config = read(ROOT/'configs/patch_sensitivity.json')
    base = {**assumptions(), 'factory_choices': delivery_settings()['factory_choices']}
    scenarios = [('nisq', 'effective', {}), ('nisq', 'pec_two_qubit', {'nisq_model': 'pec_two_qubit'}),
                 ('nisq', 'compiled_pec', {'nisq_model': 'compiled_pec'}),
                 ('star', 'reference', {}), ('hybrid', 'reference', {}), ('ft', 'reference', {})]
    costs = []
    for eps in config['precisions']:
        for tau in config['times']:
            for L in config['sizes']:
                for arm, scenario, changes in scenarios:
                    for policy in ('3x3_reference', 'added_patches'):
                        C = None
                        if policy == 'added_patches' and L != 2:
                            C = extended_coefficient(tau, algorithm_for(arm), rows, config['coefficient_margin'])
                        v = optimize(L, tau, arm, {**base, **changes}, eps, C_override=C)
                        v.update(scenario=scenario, patch_policy=policy)
                        if C is not None:
                            v.update(coefficient_patch='envelope of 3x3, 4x3, 4x4',
                                     coefficient_policy='running maximum over sampled times and step counts')
                        costs.append(v)
            print('patch sensitivity', eps, tau, flush=True)
    return costs


def capacity_rows(costs):
    config = read(ROOT/'configs/patch_sensitivity.json'); result = []
    for eps in config['precisions']:
        for tau in config['times']:
            for policy in ('3x3_reference', 'added_patches'):
                for Q in config['qubit_budgets']:
                    row = dict(absolute_error=eps, tau=tau, patch_policy=policy, Q=Q)
                    for arm, scenario in [('nisq', 'effective'), ('nisq', 'pec_two_qubit'),
                                          ('nisq', 'compiled_pec'), ('star', 'reference'),
                                          ('hybrid', 'reference'), ('ft', 'reference')]:
                        group = [v for v in costs if v['absolute_error'] == eps and v['tau'] == tau
                                 and v['patch_policy'] == policy and v['arm'] == arm and v['scenario'] == scenario]
                        row[scenario if arm == 'nisq' else arm] = capacity(group, Q)
                    result.append(row)
    return result


def plot(rows, coefficients, costs, output_dir=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullFormatter
    import numpy as np
    from plot_capacity import COLORS
    output_dir = ROOT/'figures' if output_dir is None else output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    small = read(ROOT/'results/discussion_v3/calibration.json')['rows']
    rectangles = read(ROOT/'results/discussion_v3/rectangle_checks/calibration.json')['rows']
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6), sharey=True)
    for ax, algorithm in zip(axes, ('plaquette', 'termwise')):
        for patch, color in [('3x3', '#888888'), ('4x3', '#bb772d'), ('4x4', '#14667e')]:
            group = [v for v in (small if patch == '3x3' else rectangles+coefficients)
                     if v['patch'] == patch and v['algorithm'] == algorithm and v['tau'] in (1, 2, 4)]
            group = sorted(group, key=lambda v: v['tau'])
            ax.plot([v['tau'] for v in group], [v['coefficient'] for v in group], 'o-', color=color,
                    label=patch.replace('x', ' × '))
        ax.plot([1, 2, 4], [extended_coefficient(t, algorithm, rows) for t in (1, 2, 4)],
                '--', color='#222222', label='Added-patch planning envelope')
        ax.set(yscale='log', xlabel='Evolution time τ', title=algorithm.capitalize(), xticks=[1, 2, 4])
        ax.set_ylim(.05, 8); ax.set_yticks([.1, .3, 1, 3], ['0.1', '0.3', '1', '3'])
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel('Finite-patch coefficient |Dᵣ − Dexact| r²')
    fig.suptitle('Finite-patch corrections depend on the Trotter ordering', x=.08, ha='left')
    fig.text(.5, .01, 'Solid points: converged tail pairs. Dashed: all sampled steps/times with 10% margin; not a bulk bound.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .05, 1, .94))
    for ext in ('pdf', 'png'): fig.savefig(output_dir/f'patch_time_extension.{ext}', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.2), sharex=True, sharey=True)
    for axrow, eps in zip(axes, (.001, .01)):
        for ax, tau in zip(axrow, (1, 2, 4)):
            for arm, scenario, label in [('nisq', 'effective', 'NISQ effective'),
                                         ('nisq', 'pec_two_qubit', 'NISQ PEC (2q)'),
                                         ('nisq', 'compiled_pec', 'NISQ compiled PEC'),
                                         ('star', 'reference', 'STAR'), ('hybrid', 'reference', 'Hybrid'),
                                         ('ft', 'reference', 'FTQC')]:
                for policy, ls in [('3x3_reference', ':'), ('added_patches', '-')]:
                    group = [v for v in costs if v['arm'] == arm and v['scenario'] == scenario
                             and v['absolute_error'] == eps and v['tau'] == tau and v['patch_policy'] == policy]
                    Q = sorted(set(np.geomspace(30, 1e9, 180)) | {10**v['log10_Q_min'] for v in group
                              if v['log10_Q_min'] is not None and np.log10(30) <= v['log10_Q_min'] <= 9})
                    L = np.array([capacity(group, q) for q in Q], float); L[L == 0] = np.nan
                    color = {'pec_two_qubit': '#aaaaaa', 'compiled_pec': '#444444'}.get(scenario, COLORS[arm])
                    ax.step(Q, L, where='post', color=color,
                            ls=ls, lw=1.7, label=label if policy == 'added_patches' else None)
            ax.set(xscale='log', xlim=(30, 1e9), ylim=(0, 50), title=f'τ = {tau}, ε = {eps:g}')
            ax.set_xticks([1e2, 1e4, 1e6, 1e8]); ax.set_yticks([2, 8, 16, 24, 32, 40, 48]); ax.grid(alpha=.2)
    axes[0, 0].legend(frameon=False, fontsize=8, loc='upper left')
    for ax in axes[-1]: ax.set_xlabel('Available physical qubits, Q')
    for ax in axes[:, 0]: ax.set_ylabel('Largest modeled side L (L × L sites)')
    fig.legend([Line2D([], [], color='#555555', ls=ls) for ls in (':', '-')],
               ['3 × 3 planning reference', 'Envelope including 4 × 3 and 4 × 4'],
               loc='lower center', bbox_to_anchor=(.5, .057), ncol=2, frameon=False, fontsize=9)
    fig.suptitle('Capacity sensitivity to the finite-patch accuracy model', x=.07, ha='left')
    fig.text(.5, .01, 'One week · p = 10⁻³ · same reference delivery and source menu · budgets reoptimized · sparse L grid to 48\n'
             'L = 2 retains its own calibration. Larger-L transfer and NISQ residual-bias affordability remain assumptions.',
             ha='center', fontsize=9, color='#555555')
    fig.tight_layout(rect=(0, .115, 1, .95))
    for ext in ('pdf', 'png'): fig.savefig(output_dir/f'patch_capacity.{ext}', dpi=180)
    plt.close(fig)


def source_paths():
    paths = [ROOT/name for name in ['python/dynamics_followup_analysis.py', 'python/run_dynamics_followup.py',
        'python/large_dynamics.py', 'cpp/hubbard_sector.cpp', 'python/budget_optimization.py',
        'python/discussion_model.py', 'python/delivery_sensitivity.py', 'configs/dynamics_followup.json',
        'configs/patch_sensitivity.json', 'configs/discussion_model.json', 'configs/assumptions.json',
        'configs/budget_optimization.json', 'configs/delivery_sensitivity.json',
        'pbs/dynamics_followup.pbs', 'results/dynamics_followup_v1/inputs.sha256',
        'results/dynamics_4x4/interpolated.json', 'results/dynamics_4x4/build.json',
        'results/discussion_v3/calibration.json', 'results/discussion_v3/finite_step_guard.json',
        'results/discussion_v3/rectangle_checks/calibration.json']]
    for directory in [ROOT/'results/dynamics_4x4/tau1', ROOT/'results/dynamics_4x4/rectangle_tau1']+list(OUT.glob('tau*')):
        if directory.is_dir(): paths.extend(p for p in directory.iterdir() if p.is_file())
    paths.extend(OUT.glob('*_completion.json'))
    return sorted(set(paths))


def summary(rows, checks, direct, coefficients, costs):
    text = ['# Direct allocation checks and longer-time 4 × 4 dynamics', '',
        f"All {len(checks)} original headline allocations pass their saved Trotter reserves, before reoptimization.",
        'The reference precision remains epsilon=0.001 (1% of D*=0.1), with 0.01 as the loose sensitivity.',
        'This checks finite-lattice Trotter error only. NISQ mitigation bias and experimental-ordering transfer remain assumptions.', '',
        '| Error | Arm | Layout | Saved r | Error / saved reserve | Saved Q | Q with measured-error budget at saved r | Best measured r | Best measured Q |',
        '|---:|---|---|---:|---:|---:|---:|---:|---:|']
    def area(v):
        q = v['log10_Q_min']; return f'{10**q:.5g}' if q < 9 else f'10^{q:.2f}'
    for v in checks:
        pair = [x for x in direct if all(x[k] == v[k] for k in ('arm', 'layout', 'absolute_error'))]
        same = next(x for x in pair if x['step_menu'] == 'saved_r')
        best = next(x for x in pair if x['step_menu'] == 'all_measured_r')
        text.append(f"| {v['absolute_error']} | {v['arm']} | {v['layout']} | {v['r']} | {v['error_to_reserve']:.4f} | {area(v)} | {area(same)} | {best['r']} | {area(best)} |")
    text += ['', 'The original Q column is validated with its unchanged allocation. Subsequent columns reoptimize the allocation using measured error plus 1e-10. The last column searches only the measured step menu; it is not a global optimum.', '',
        '| Patch | Time | Ordering | Exact central doublon | Tail coefficient | Tail spread | Largest sampled coefficient |',
        '|---|---:|---|---:|---:|---:|---:|']
    for v in coefficients:
        text.append(f"| {v['patch']} | {v['tau']:g} | {v['algorithm']} | {v['exact_doublon']:.10f} | {v['coefficient']:.6g} | {v['tail_spread']:.3%} | {v['largest_measured_coefficient']:.6g} |")
    text += ['', '## Capacity sensitivity', '',
        'The new figure holds delivery and source choices fixed while changing the accuracy coefficient. Each curve reoptimizes error allocations and steps. The added-patch envelope retains the old 3x3 running maximum and adds every measured 4x3/4x4 coefficient through the requested time, with a 10% margin. L=2 keeps its own calibration.',
        'This envelope is a planning sensitivity, not a bulk bound, a continuously sampled time maximum, or a validation at every selected r. The main discussion-v3 curves are preserved, and no 4x4 inference at tau=8 is made.', '',
        '[Capacity figure](../../figures/patch_capacity.pdf) · [Patch/time figure](../../figures/patch_time_extension.pdf)', '',
        '| Time | Ordering | Old planning coefficient | Added-patch coefficient |', '|---:|---|---:|---:|']
    for tau in (1, 2, 4):
        for algorithm in ('plaquette', 'termwise'):
            text.append(f'| {tau} | {algorithm} | {coefficient(tau, algorithm):.6g} | {extended_coefficient(tau, algorithm, rows):.6g} |')
    (OUT/'summary.md').write_text('\n'.join(text)+'\n')


def run():
    for task in ('selected', 'tau2', 'tau4'):
        assert read(OUT/f'{task}_completion.json')['status'] == 'complete'
    rows = measurements(); completion_checks(rows); checks = saved_allocation_checks(rows)
    write('measurements.json', rows); write('allocation_checks.json', checks)
    assert all(v['passed'] for v in checks), 'A saved allocation failed: inspect before reoptimization'
    direct = direct_costs(rows); coefficients = patch_coefficients(rows); costs = sensitivity(rows)
    write('direct_costs.json', direct); write('patch_coefficients.json', coefficients); write('costs.json', costs)
    write('assumptions.json', {'hardware_and_delivery': {**assumptions(),
          'factory_choices': delivery_settings()['factory_choices']},
          'budget_search': budget_settings(), 'dynamics': settings(),
          'patch_sensitivity': read(ROOT/'configs/patch_sensitivity.json'),
          'step_selection': 'Measured coefficients and budget optimization override the legacy fixed-allocation step defaults'})
    write('capacity.json', capacity_rows(costs))
    with (OUT/'capacity.csv').open('w', newline='') as f:
        data = capacity_rows(costs); writer = csv.DictWriter(f, fieldnames=list(data[0]))
        writer.writeheader(); writer.writerows(data)
    with (OUT/'direct_budgets.csv').open('w', newline='') as f:
        columns = ['absolute_error', 'arm', 'layout', 'step_menu', 'r', 'original_r',
                   'log10_Q_min', 'original_log10_Q_min', 'measured_trotter_error']
        data = [{**{k: v[k] for k in columns}, **v['budget_fractions']} for v in direct]
        writer = csv.DictWriter(f, fieldnames=list(data[0])); writer.writeheader(); writer.writerows(data)
    plot(rows, coefficients, costs); summary(rows, checks, direct, coefficients, costs)
    write('validation.json', {'status': 'PASS: measured step allocations, converged finite-patch coefficients, and capacity accounting',
        'allocation_checks': len(checks), 'measurements': len(rows), 'cost_rows': len(costs),
        'scope': settings()['scope'], 'sha256': {str(p.relative_to(ROOT)): sha(p) for p in source_paths()}})
    check(require_jobs=False)


def check(require_jobs=True):
    report = read(OUT/'validation.json')
    for path, value in report['sha256'].items(): assert sha(ROOT/path) == value, path
    rows = measurements(); completion_checks(rows); checks = saved_allocation_checks(rows)
    assert rows == read(OUT/'measurements.json')
    assert checks == read(OUT/'allocation_checks.json') and len(checks) == 16
    assert all(v['passed'] for v in checks)
    assert patch_coefficients(rows) == read(OUT/'patch_coefficients.json')
    direct = read(OUT/'direct_costs.json'); costs = read(OUT/'costs.json')
    base = {**assumptions(), 'factory_choices': delivery_settings()['factory_choices']}
    from discussion_model import estimate
    for v in direct+costs:
        if v['log10_Q_min'] is None: continue
        b = v['budget_fractions']; assert isclose(sum(b.values()), 1, abs_tol=1e-12)
        assert min(b.values()) >= 0 and b['stat'] > 0
        changes = layout_changes(v['arm'], v.get('layout', v.get('scenario')))
        c = {**base, **changes, 'stat_fraction': b['stat'], 'logical_fraction': b['logical'],
             'magic_fraction': b['magic'], 'synthesis_fraction': b['synthesis'],
             'forced_cleanup': v.get('cleanup_levels', 0)}
        replay = estimate(v['L'], v['tau'], v['arm'], c, v['absolute_error'], steps=v['r'])
        assert isclose(replay['log10_Q_min'], v['log10_Q_min'], abs_tol=1e-12)
        eps = v['absolute_error']; seconds = v['seconds_per_attempt']
        shots = ceil((log(2/c['delta'])/2 if v['shot_rule'] == 'hoeffding' else
                      c['gaussian_z']**2*c['signal_reference']*(1-c['signal_reference']))/(eps*b['stat'])**2)
        assert isclose(log10(shots)+v['chi']/log(10), v['log10_total_shots'], abs_tol=1e-10)
        copies = max(0, log(shots)+v['chi']+log(seconds/c['week_seconds']))
        if copies < 30: copies = log(ceil(exp(copies)))
        assert isclose((log(v['qubits_per_copy'])+copies)/log(10), v['log10_Q_min'], abs_tol=1e-10)
        if v['arm'] != 'nisq':
            d = v['distance']; f = v.get('factories', 0)
            tiles = (v['qubits_per_copy']-f*c['source_qubits'])/(2*d*d)
            faults = log(tiles*ceil(seconds/c['cycle_seconds'])*c['surface_A'])
            faults += (d+1)/2*log(c['p']/c['surface_threshold'])
            assert faults+v['chi']/2 <= log(eps*b['logical'])+1e-12
        if v['arm'] in ('ft', 'hybrid'):
            assert log(v['T_count']*v['magic_error'])+v['chi']/2 <= log(eps*b['magic'])+1e-12
            assert seconds >= max(v['supply_seconds'], v['consumption_seconds'])
        assert isclose(v['absolute_error']*b['trotter'], v['coefficient']/v['r']**2, abs_tol=1e-12)
        if 'measured_trotter_error' in v:
            assert isclose(v['coefficient']/v['r']**2, v['measured_trotter_error']+v['numerical_allowance'], abs_tol=1e-12)
        assert v.get('factories', 0) < max(delivery_settings()['factory_choices'])
    assert capacity_rows(costs) == read(OUT/'capacity.json')
    for algorithm in ('plaquette', 'termwise'):
        C = [extended_coefficient(tau, algorithm, rows) for tau in (1, 2, 4)]
        assert C == sorted(C)
        assert all(a >= coefficient(tau, algorithm) for a, tau in zip(C, (1, 2, 4)))
    if require_jobs:
        jobs = read(OUT/'jobs.json'); submissions = read(OUT/'submissions.json')
        assert {v['jobid'] for v in jobs} == {v['jobid'] for v in submissions}
        assert all(v['job_state'] == 'F' and v['Exit_status'] == 0 for v in jobs)
        job_for = {v['task']: v['jobid'] for v in submissions}
        for path in OUT.glob('tau*/*.json'):
            task = 'selected' if path.parent.name == 'tau1' else path.parent.name
            assert read(path)['pbs_jobid'] == job_for[task]
    for line in (OUT/'inputs.sha256').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        if name == 'build/hubbard_sector':
            assert digest == read(ROOT/'results/dynamics_4x4/build.json')['sha256'][name]
        else:
            assert sha(ROOT/name) == digest, name
    print(report['status'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    if parser.parse_args().check: check()
    else: run()
