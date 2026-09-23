"""Replay published ledgers and compare experimental observables on their own tasks."""
import argparse
import csv
import hashlib
import json
from math import ceil, exp, floor, log10, sqrt

import numpy as np
from scipy.optimize import brentq

from project_paths import ROOT

OUT = ROOT/'results/published_benchmarks_v1'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(name, data):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/name).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def csv_rows(name, rows):
    fields = list(dict.fromkeys(k for v in rows for k in v))
    with (OUT/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def qcels_branches(r, J, N=5):
    """All nonzero times in Appendix C, n=1,...,N, with integer Trotter steps."""
    return [(n/N*2.**(j-J), ceil(r*n/N*2.**(j-J)))
            for j in range(1, J+1) for n in range(1, N+1)]


def qcels_time(row, cycle, ns=100, quadratures=2, growth=0., beta=0., control=True):
    L, U = row['L'], row['U']
    terms = 8*L*(L-1)+L*L
    norm = 4*L*(L-1)+U*L*L/4
    Tmax = .06/(2*row['epsilon']/3)
    clock_seconds = (row['di']+row['df']+6)*cycle
    total = 0.
    for fraction, steps in qcels_branches(row['r'], row['J']):
        clocks = row['clocks']*steps + (18*steps+16 if control else 0)
        # Eq. 11: every rotation has P_L = beta*|theta| + 2 I_growth.
        chi = 4*beta*norm*Tmax*fraction/2 + 8*growth*(2*terms*steps)
        total += clocks*clock_seconds*exp(chi)
    return ns*quadratures*total


def topqad_replay(cfg):
    assert (cfg['N'], cfg['Ns'], cfg['delta'], cfg['quadratures']) == (5, 100, .06, 2), \
        'This replay implements the published Table IV sampling conventions'
    star = []
    for v in cfg['star']:
        L, U, df = v['L'], v['U'], v['df']
        base = qcels_time(v, cfg['cycle_seconds'])
        beta = cfg['hardware'][v['hardware']]['beta_d9']
        observed_multiplier = v['seconds_published']/v['pec_free_seconds_published']
        angle_multiplier = qcels_time(v, cfg['cycle_seconds'], beta=beta)/base
        # This is an explicitly inverse diagnostic, not a fit used to claim
        # independent reproduction of the paper's unknown growth parameters.
        growth = None
        if observed_multiplier >= angle_multiplier:
            growth = brentq(lambda g: qcels_time(v, cfg['cycle_seconds'], growth=g, beta=beta)/base-observed_multiplier,
                            0., 1e-4, xtol=1e-18)
        theta_total = (4*L*(L-1)+U*L*L/4)*(.06/(2*v['epsilon']/3))/2
        single = v['r']*v['clocks']*(v['di']+df+6)*cfg['cycle_seconds']
        row = dict(v, qubits_replayed=(4*L*L+1)*2*df*df,
            rotations_replayed=2*(8*L*(L-1)+L*L)*v['r'],
            r_from_rounded_w=.06/(2*(2*v['epsilon']/3))*sqrt(v['w']/(v['epsilon']/3)),
            deepest_total_angle=theta_total, single_seconds_replayed=single,
            pec_free_seconds_replayed=base, published_variance_multiplier=observed_multiplier,
            angle_only_multiplier=angle_multiplier, implied_growth_per_state=growth,
            clock_ratio=single/v['single_seconds_published'],
            pec_free_ratio=base/v['pec_free_seconds_published'],
            qubit_ratio=(4*L*L+1)*2*df*df/v['qubits_published'],
            cycle_only_240ns_seconds=base*240e-9/cfg['cycle_seconds'])
        flags = []
        if abs(row['clock_ratio']-1) > .05:
            flags.append('Published clock and single-shot columns disagree')
        if abs(row['pec_free_ratio']-1) > .05:
            flags.append('PEC-free replay differs by more than 5 percent')
        if growth is None:
            flags.append('Published PEC multiplier is below the stated small-angle contribution')
        row['flags'] = '; '.join(flags) or 'none'
        star.append(row)
    # Calibrate one growth parameter on U=4,L=4 for each hardware, then carry
    # it unchanged to the other seven rows. This is more informative than
    # claiming success from inverting every row independently.
    for v in star:
        anchor = next(a for a in star if (a['U'],a['L'],a['hardware'])==(4,4,v['hardware']))
        g = anchor['implied_growth_per_state']
        beta = cfg['hardware'][v['hardware']]['beta_d9']
        seconds = qcels_time(v,cfg['cycle_seconds'],growth=g,beta=beta)
        v.update(anchor_growth_per_state=g,anchor_transfer_seconds=seconds,
                 anchor_transfer_ratio=seconds/v['seconds_published'],
                 used_to_fit_growth=(v['U'],v['L'])==(4,4))
    ft = []
    for v in cfg['ft']:
        seconds = v['T_count']*v['df']*cfg['cycle_seconds']
        ft.append(dict(v, one_port_seconds=seconds, runtime_ratio=seconds/v['seconds_published'],
                       effective_T_clocks=v['seconds_published']/seconds,
                       algorithm_data_qubits_only=2*v['df']**2*v['n_logical']))
    return star, ft


def fixed_workload_sensitivities(star,cfg):
    """One changed input at a time, without calling this a hardware migration."""
    reference=next(v for v in star if (v['U'],v['L'],v['hardware'])==(4,4,'target'))
    beta=cfg['hardware']['target']['beta_d9'];g=reference['anchor_growth_per_state']
    growth_data=read(ROOT/'results/frontier_v1/growth_summary.json')['pairs']
    g_upper=next(v['pauli_total_upper'] for v in growth_data if (v['di'],v['df'],v['p'])==(9,13,.001))
    rows=[]
    for name,cycle,b,growth in [('Published anchor',350e-9,beta,g),
          ('Clock only: 240 ns',240e-9,beta,g),('Angle coefficient only: doubled',350e-9,2*beta,g),
          ('Growth only: our d9-to-d13 upper bound',350e-9,beta,g_upper)]:
        seconds=qcels_time(reference,cycle,growth=growth,beta=b)
        rows.append(dict(scenario=name,seconds=seconds,relative_to_anchor=seconds/reference['anchor_transfer_seconds'],
                         cycle_seconds=cycle,beta=b,growth_per_state=growth))
    return rows


def nearest_rows(rows, target_times):
    """Use only matching saved times; never smooth or interpolate the benchmark."""
    found = []
    for t in target_times:
        matches = [v for v in rows if abs(v['tau']-t) < 1e-8]
        if len(matches) == 1:
            found.append(matches[0])
    return found


def phasecraft_summary(data):
    raw_audit, comparisons = [], []
    observables = data['observables']
    for U in (0, 4):
        shots = [v for v in data['shots'] if v['U'] == U]
        times = [v['tau'] for v in shots]
        published_raw = [v for v in observables if v['U'] == U and v['observable'] == 'mean_doublon' and v['group'] == 'raw']
        raw = nearest_rows(published_raw, times)
        assert len(raw) == len(shots), 'Raw data should have one value per shot time'
        differences = [abs(a['mean_doublon']-b['value']) for a,b in zip(shots, raw)]
        raw_audit.append(dict(U=U, time_points=len(shots), total_shots=sum(v['shots'] for v in shots),
            shots_per_point=sorted(set(v['shots'] for v in shots)),
            max_raw_mean_difference=max(differences),
            mean_number_acceptance=float(np.mean([v['number_acceptance'] for v in shots])),
            median_variance_per_shot=float(np.median([v['variance_per_shot'] for v in shots])),
            median_local_variance=float(np.median([v['mean_local_bernoulli_variance'] for v in shots]))))
        references = [('FLO', None), ('Exact', None)] if U == 0 else [('TDVP chi=1024', 1024), ('TDVP chi=2048', 2048)]
        for observable in ('mean_doublon', 'local_doublon'):
          sites = sorted(set(v['site'] for v in observables if v['U']==U and v['observable']==observable),key=str)
          for ref_group, bond in references:
           for group in ('Raw', 'TFLO', 'GPR', 'Mitigated'):
            diffs, errors = [], []
            for site in sites:
                selected = [v for v in observables if v['U']==U and v['observable']==observable and v['site']==site]
                refs = nearest_rows([v for v in selected if v['method']==ref_group],times)
                method = nearest_rows([v for v in selected if v['method']==group],times)
                reference = {round(v['tau'],8):v for v in refs}
                for v in method:
                    if round(v['tau'],8) in reference:
                        diffs.append(v['value']-reference[round(v['tau'],8)]['value'])
                        if v['reported_se'] is not None: errors.append(v['reported_se'])
            if diffs:
                assert len(diffs) == len(times)*(28 if observable=='local_doublon' else 1), \
                    'Comparison must retain every experimental time and site'
                comparisons.append(dict(U=U,observable=observable,group=group,reference=ref_group,bond_dimension=bond,
                    points=len(diffs),rmse=float(np.sqrt(np.mean(np.square(diffs)))),
                    max_absolute_difference=float(np.max(np.abs(diffs))),
                    median_reported_se=float(np.median(errors)) if errors else None,
                    meaning='Training consistency, not held-out accuracy' if U==0 else 'Finite-bond reference agreement, not residual bias'))
    return raw_audit, comparisons


def hashes():
    files = ['python/published_benchmarks.py','python/phasecraft_data.py',
        'configs/published_benchmarks.json','refs/2603.13093v1.pdf','refs/2510.26300v2.pdf',
        'results/published_benchmarks_v1/phasecraft_reduced.json',
        'results/published_benchmarks_v1/data_availability.json','results/frontier_v1/growth_summary.json']
    return {name:sha(ROOT/name) for name in files}


def run(check=False):
    cfg = read(ROOT/'configs/published_benchmarks.json')
    data = read(OUT/'phasecraft_reduced.json')
    assert data['reducer_sha256'] == sha(ROOT/'python/phasecraft_data.py')
    star, ft = topqad_replay(cfg)
    raw, comparisons = phasecraft_summary(data)
    result = dict(topqad_star=star,topqad_ft=ft,phasecraft_raw=raw,phasecraft_comparisons=comparisons,
                  fixed_workload_sensitivities=fixed_workload_sensitivities(star,cfg))
    validation = dict(status='PASS: published-input ledger replay and reduced experimental data checks',
        star_rows=len(star),ft_rows=len(ft),experimental_shot_batches=len(data['shots']),
        observable_rows=len(data['observables']),comparison_rows=len(comparisons),
        flagged_star_rows=[dict(U=v['U'],L=v['L'],hardware=v['hardware'],flags=v['flags']) for v in star if v['flags']!='none'],
        max_star_area_relative_difference=max(abs(v['qubit_ratio']-1) for v in star),
        max_ft_time_optimal_relative_difference=max(abs(v['runtime_ratio']-1) for v in ft if v['mode']=='time'),
        scope='Table inputs are not an independently reproduced injection schedule. Implied growth is diagnostic. Experimental processing is supplied, not refitted; U=4 references are not exact truth.')
    assert max(v['max_raw_mean_difference'] for v in raw) < 1e-12
    # Table IV reports two significant digits; 1.0e5 represents a 5000-wide
    # half-bin, so a blanket three-percent tolerance would reject rounding.
    assert all(abs(v['qubits_replayed']-v['qubits_published']) <=
               .5*10**(floor(log10(v['qubits_published']))-1) for v in star)
    assert validation['max_ft_time_optimal_relative_difference'] < .03
    assert len(star)==16 and len(ft)==32 and len(raw)==2 and len(comparisons)>0
    if check:
        assert read(OUT/'inputs.json') == hashes()
        assert read(OUT/'benchmarks.json') == result
        assert read(OUT/'validation.json') == validation
        print(validation['status'])
    else:
        save('inputs.json',hashes());save('benchmarks.json',result);save('validation.json',validation)
        csv_rows('topqad_star.csv',star);csv_rows('topqad_ft.csv',ft)
        csv_rows('phasecraft_comparisons.csv',comparisons)
        print(json.dumps(validation,indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    run(parser.parse_args().check)
