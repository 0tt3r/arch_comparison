"""Accuracy/cost boundaries and conditional component composition controls."""
import json
from math import exp, log, sqrt

import numpy as np
import pytest

from integrated_comparison import (OUT, component_exposure, configuration,
    encoded_candidate, growth_probability, nisq_attempt, offset_limit,
    planned_runtime, rotation_channel, threshold)
from model import workload


def test_accuracy_floor_cannot_be_bought_away_with_shots():
    assert threshold(10**20, .01-.03) is None
    assert offset_limit(100,.01,100) is None
    assert offset_limit(10**20,.01,100) == pytest.approx(.01,abs=1e-8)


def test_break_even_boundary_matches_independent_gaussian_estimator():
    # A Bernoulli signal with affine slope s has variance s^2 D(1-D)/M.
    M,eps,b,trotter = 200000,.01,.002,.001
    F = threshold(M,eps-b-trotter)
    slope = sqrt(F)
    variance = slope*slope*.1*.9/M
    assert b+trotter+1.96*sqrt(variance) == pytest.approx(eps)
    assert offset_limit(M,eps-trotter,F) == pytest.approx(b)
    assert threshold(10*M,eps-b-trotter) == pytest.approx(10*F)


def test_training_and_postselection_cost_conventions():
    raw,training,acceptance = 100000,10000,.2
    # Rejection has to appear once: per-attempt variance is V/acceptance.
    F_allowed = threshold(raw-training,.005)
    accepted_bound = (raw-training)*acceptance*.005**2/(1.96**2*.09)
    assert F_allowed*acceptance == pytest.approx(accepted_bound)
    assert threshold(0,.01) is None


def test_growth_does_not_invent_unmeasured_distances_or_zero_error():
    assert growth_probability(5,.001,'component') > 5e-4
    assert 1e-5 < growth_probability(9,.001,'component') < 1.3e-5
    assert growth_probability(25,.001,'component') == 0
    assert growth_probability(9,.001,'zero_growth') == 0
    with pytest.raises(ValueError):
        growth_probability(5,.001,'zero_growth')
    with pytest.raises(StopIteration):
        growth_probability(7,.001,'component')


def test_noisy_channel_inverse_has_independent_matrix_norm_control():
    v = rotation_channel(.01,9,.001,'component')
    q = v['q_dephasing']
    assert q > rotation_channel(.01,9,.001,'zero_growth')['q_dephasing']
    w = workload(2,1,100,'termwise')
    c = component_exposure(w,.001,9,'component')
    # The x/y block of each channel divided by its real multiplier leaves
    # I + a rotation-generator perturbation. Its product residual is bounded
    # by the product of (1+|imag/real|) used in the campaign ledger.
    product = np.eye(2)
    for theta,n in w.dynamic_angles:
        ch = rotation_channel(theta,9,.001,'component')
        real,imag = ch['multiplier']
        matrix = np.array([[real,-imag],[imag,real]])/real
        product = np.linalg.matrix_power(matrix,n) @ product
    assert np.linalg.norm(product-np.eye(2),2)/2 <= c['coherent_bias_bound'] + 2e-12
    assert c['serial_analog_seconds'] > 0 and c['rus_tail_bound'] >= 0


def test_replica_floor_and_minimum_attempt_time():
    row = dict(qubits_per_copy=100,seconds_per_attempt=2.,log10_total_shots=3.,log10_Q_min=2.)
    assert planned_runtime(row,99) is None
    assert planned_runtime(row,199) == pytest.approx(2000)
    assert planned_runtime(row,200) == pytest.approx(1000)
    assert planned_runtime(row,10**9) == 2.


def test_lower_distance_anchor_uses_its_own_growth_and_service():
    smaller = rotation_channel(.01,9,.0001,'component',13)
    larger = rotation_channel(.01,9,.0001,'component',25)
    assert smaller['data_distance']==13 and larger['data_distance']==25
    assert smaller['growth_probability']==growth_probability(9,.0001,'component',13)
    assert smaller['mean_seconds'] < larger['mean_seconds']


def test_analog_candidate_includes_growth_receiver_and_real_service():
    w = workload(2,1,16,'termwise')
    c = configuration(.001)
    component = dict(component_exposure(w,.001,9,'component'),seed=9)
    b = dict(stat=.5,logical=.1)
    row = encoded_candidate(w,'star',c,.01,component,b,0,'compact')
    assert row['distance'] >= 25
    assert row['analog_seconds'] >= component['serial_analog_seconds']/c['parallel_lanes']
    assert row['seconds_per_attempt'] >= row['analog_seconds']
    assert exp(row['log_logical_fault_bound']+row['chi']/2) <= .01*b['logical']


def test_saved_integration_and_requirement_contracts_when_available():
    if not (OUT/'validation.json').exists():
        pytest.skip('Resource sweep not generated yet')
    from integrated_comparison import load_channels, validate
    load_channels()
    rows = json.loads((OUT/'costs.json').read_text())
    thresholds = json.loads((OUT/'break_even.json').read_text())
    assert validate(rows,thresholds)['status'].startswith('PASS:')
