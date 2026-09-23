"""Read-only checks for docs/RAW_RUS_NOTE.md; no result/cache files are written.

Run with qb_gsee Python. Default output is the compact numerical summary;
--json includes every angle, input hashes, and the explicit failure controls.
The synthetic controls are sensitivity examples, not measured gadget channels.
"""
import sys
sys.dont_write_bytecode = True

import argparse
import hashlib
import json
from math import expm1, fsum, log, log1p, pi
from pathlib import Path

import numpy as np

from model import workload
from star_channel_study import (growth_error, preparation, raw_rus_details,
                                rus_channel, wrap)

ROOT = Path(__file__).resolve().parents[1]
CAP = 48


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def clear_caches():
    preparation.cache_clear()
    raw_rus_details.cache_clear()


def branch_sum(cs, success=0.5, conjugate_failure=True):
    """Condition on completion before len(cs); only scalar-channel toy models."""
    reached, failed, total = 1.0, 1.0 + 0j, 0j
    for c in cs:
        total += reached * success * failed * c
        failed *= c.conjugate() if conjugate_failure else c
        reached *= 1 - success
    return total / (1 - reached), reached


def grid():
    angles = sorted(set(np.geomspace(1e-3, pi/8, 25).tolist() + [.01, pi/16]))
    rows = []
    cap_change = 0.0
    for d in (5, 9, 25):
        for p in (1e-3, 1e-4):
            for theta in angles:
                prep = preparation(theta, d, p)
                raw = rus_channel(theta, d, p, 'raw', floor=0, cap=CAP)
                pcec = rus_channel(theta, d, p, 'noisy_pcec', floor=0, cap=CAP)
                deeper = rus_channel(theta, d, p, 'raw', floor=0, cap=64)
                z = complex(*raw['multiplier'])
                change = abs(z - complex(*deeper['multiplier']))
                cap_change = max(cap_change, change)
                independent, tail = branch_sum([
                    preparation(wrap(theta * 2**m), d, p)['multiplier']
                    for m in range(CAP)])
                check(abs(z - independent) < 2e-13, 'raw branch sum disagrees')
                check(change < 2e-13, 'cap 48 is not converged numerically')
                for result in (raw, pcec):
                    check(0 <= result['q_dephasing'] < 0.5, 'invalid dephasing')
                    check(abs(complex(*result['multiplier'])) <= 1 + 1e-13,
                          'noncontractive multiplier')
                # This is the implemented *single-trial* ideal-PCEC identity.
                c = prep['multiplier']
                ideal_q = (1 - abs(c)**2)/2
                q = (1 - c.real)/2
                check(abs(ideal_q - (2*q - 2*q*q - c.imag*c.imag/2)) < 1e-15,
                      'ideal-PCEC algebra')
                rows.append(dict(theta=theta, d=d, p=p,
                    preparation_acceptance=prep['acceptance'],
                    detector_adjusted_acceptance=prep['acceptance']*(1-p)**(d*d+6*d),
                    raw_q=raw['q_dephasing'], pcec_q=pcec['q_dephasing'],
                    raw_imag=raw['coherent_component'],
                    pcec_imag=pcec['coherent_component'],
                    raw_chi=-2*log(1-2*raw['q_dephasing']),
                    pcec_chi=-2*log(1-2*pcec['q_dephasing']),
                    raw_tail=raw['uncorrected_rus_tail'],
                    pcec_tail=pcec['uncorrected_rus_tail']))
                clear_caches()
    summary = []
    for d in (5, 9, 25):
        for p in (1e-3, 1e-4):
            selected = [v for v in rows if (v['d'], v['p']) == (d, p)]
            span = lambda key: [min(v[key] for v in selected), max(v[key] for v in selected)]
            ratios = [v['pcec_chi']/v['raw_chi'] for v in selected]
            summary.append(dict(d=d, p=p, raw_q_range=span('raw_q'),
                pcec_q_range=span('pcec_q'),
                max_abs_raw_imag=max(abs(v['raw_imag']) for v in selected),
                max_abs_pcec_imag=max(abs(v['pcec_imag']) for v in selected),
                chi_ratio_range=[min(ratios), max(ratios)],
                acceptance_range=span('preparation_acceptance')))
    return rows, summary, cap_change


def controls():
    theta, d, p = .01, 9, .001
    cs = [preparation(wrap(theta*2**m), d, p)['multiplier'] for m in range(64)]
    symmetric, _ = branch_sum(cs)
    biased, _ = branch_sum(cs, .501)
    linear_biased = fsum((1-.501)**m * (2*.501-1) * c.imag
                        for m, c in enumerate(cs))
    check(abs(biased.imag-linear_biased) < 1e-8,
          'biased-success first-order sensitivity')
    caps = []
    for k in (4, 8, 16, 24, 32, 48):
        z, tail = branch_sum(cs[:k])
        first_order = tail/(1-tail)*fsum(c.imag for c in cs[:k])
        # The coefficient of each b_m is exactly the omitted terminal mass.
        check(all(abs(2**(-m-1)-fsum(2**(-n-1) for n in range(m+1,k))-tail)
                  < 1e-16 for m in range(k)), 'finite-cap cancellation identity')
        check(abs(z-symmetric) <= 2*tail + 3e-14, 'conditional tail bound')
        caps.append(dict(cap=k, tail=tail, imag=z.imag,
                         first_order_imag=first_order,
                         change_from_cap64=abs(z-symmetric)))

    pulse_error = 1e-4
    c = complex(np.exp(-2j*pulse_error))
    conjugate = .5*c/(1-.5*c.conjugate())
    same_sign = .5*c/(1-.5*c)
    check(abs(conjugate.imag) < 1e-10 and abs(same_sign.imag) > 3.9e-4,
          'systematic sign control')
    finite_angle = .2
    rare = complex(np.exp(-2j*finite_angle))
    correlated = 1-p+p*(.5*rare/(1-.5*rare.conjugate()))
    single = 1-p+p*rare
    independent = .5*single/(1-.5*single.conjugate())
    check(abs(correlated.imag) > 100*abs(independent.imag),
          'rare correlated-fault counterexample')
    acceptances = []
    for seed in (5, 9, 25):
        A = preparation(theta, seed, p)['acceptance']
        # Sensitivity only: an independent omitted fault occurs with probability
        # u, always passes selection, and maximally polarizes the accepted state.
        # Without that fault, selection succeeds with the modeled probability A.
        u = p
        acceptances.append(dict(d=seed, theta=theta, p=p, acceptance=A,
            mean_preparations=1/A,
            illustrative_omitted_fault_u=u,
            upper_success_imbalance_if_omitted_fault=u/(2*(u+(1-u)*A))))
    for seed in (5, 9, 25):
        z = rus_channel(.013, seed, 0., 'raw', cap=CAP)
        check(abs(complex(*z['multiplier'])-1) < 2e-13, 'zero-noise control')
    clear_caches()
    return dict(biased_success=dict(theta=theta,d=d,p=p,success=.501,
                    balanced_imag=symmetric.imag,biased_imag=biased.imag,
                    first_order_biased_imag=linear_biased),
                caps=caps, acceptance_examples=acceptances,
                systematic=dict(pulse_error=pulse_error,
                    conjugate_imag=conjugate.imag,same_sign_imag=same_sign.imag),
                correlated=dict(p=p,finite_angle=finite_angle,
                    correlated_imag=correlated.imag,independent_imag=independent.imag))


def workload_examples():
    """Fixed r: reproduce the review's two exposures, without reoptimization."""
    pairs = json.loads((ROOT/'results/frontier_v1/growth_summary.json').read_text())['pairs']
    g = next(v['pauli_total_upper'] for v in pairs
             if (v['di'],v['df'],v['p']) == (9,13,.001))
    anchor = g/(growth_error(9,.001,1e-10,13)/1e-10)
    rows = []
    for L,r in ((4,10),(16,8)):
        w = workload(L,1,r,'termwise',4)
        for mode in ('raw','noisy_pcec'):
            chi = coherent_log = 0.
            for angle,count in w.dynamic_angles:
                v = rus_channel(angle,9,.001,mode,floor=anchor,data_distance=13,cap=CAP)
                real = 1-2*v['q_dephasing']
                chi -= 2*count*log(real)
                coherent_log += count*log1p(abs(v['coherent_component'])/real)
            rows.append(dict(L=L,r=r,mode=mode,chi=chi,
                coherent_bias_bound=.5*expm1(coherent_log),growth_upper=g))
            clear_caches()
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', action='store_true', help='include all angle rows')
    args = parser.parse_args()
    rows, summary, change = grid()
    result = dict(grid_points=len(rows),cap=CAP,correction_cap=32,
        angles='25 geometric points in [0.001, pi/8], plus 0.01 and pi/16',
        growth_in_grid=0,summary=summary,max_cap48_to64_change=change,
        controls=controls(),workloads=workload_examples(),
        source_sha256={str(path):hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
            for path in map(Path,('python/star_channel_study.py','python/model.py',
                'python/raw_rus_check.py','refs/PhysRevX.15.021057.pdf',
                'results/frontier_v1/growth_summary.json'))},
        validation='All numerical and algebraic checks passed; restricted-model evidence only.')
    if args.json:
        result['rows'] = rows
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
