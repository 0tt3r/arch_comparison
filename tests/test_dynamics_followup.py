"""Prevent reallocation and finite-patch extrapolation from masquerading as checks."""
import pytest

import dynamics_followup_analysis as analysis
from run_dynamics_followup import headline_rows, selected_cases, tail_spread


def test_frozen_headlines_share_nine_unique_circuits():
    rows = headline_rows()
    assert len(rows) == 16
    assert len(selected_cases(rows)) == 9
    assert ('termwise', 33) in selected_cases(rows)
    assert ('plaquette', 21) in selected_cases(rows)


def test_audit_uses_original_reserve_and_does_not_reallocate(monkeypatch):
    plan = dict(arm='nisq', layout='effective', absolute_error=.001, r=20,
                log10_Q_min=2, budget_fractions={'trotter': .3, 'stat': .6, 'residual_bias': .1})
    monkeypatch.setattr(analysis, 'headline_rows', lambda: [plan])
    measurement = dict(patch='4x4', tau=1, algorithm='plaquette', r=20,
                       signed_error=-.00031, source='synthetic')
    check = analysis.saved_allocation_checks([measurement])[0]
    assert not check['passed']  # Would fit by stealing a small amount from statistics.
    assert check['original_trotter_reserve'] == pytest.approx(.0003)
    assert plan['budget_fractions']['stat'] == .6
    measurement['signed_error'] = .00029
    assert analysis.saved_allocation_checks([measurement])[0]['passed']


def test_selected_case_missing_or_duplicated_is_rejected(monkeypatch):
    plan = dict(arm='star', r=33)
    monkeypatch.setattr(analysis, 'headline_rows', lambda: [plan])
    with pytest.raises(AssertionError): analysis.saved_allocation_checks([])
    row = dict(patch='4x4', tau=1, algorithm='termwise', r=33)
    with pytest.raises(AssertionError): analysis.saved_allocation_checks([row, row])


def test_completion_requires_both_orderings_at_each_time():
    incomplete = [dict(patch='4x4', tau=tau, algorithm='plaquette') for tau in (1, 2, 4)]
    with pytest.raises(AssertionError, match='Missing time/ordering'): analysis.completion_checks(incomplete)


def test_added_patch_envelope_keeps_past_peaks_and_excludes_future_data():
    rows = [dict(tau=1, algorithm='plaquette', coefficient=10),
            dict(tau=2, algorithm='plaquette', coefficient=1),
            dict(tau=4, algorithm='plaquette', coefficient=1000),
            dict(tau=1, algorithm='termwise', coefficient=2000)]
    assert analysis.extended_coefficient(1, 'plaquette', rows) == pytest.approx(11)
    assert analysis.extended_coefficient(2, 'plaquette', rows) == pytest.approx(11)
    assert analysis.extended_coefficient(4, 'plaquette', rows) == pytest.approx(1100)
    with pytest.raises(ValueError, match='No extended-patch'): analysis.extended_coefficient(8, 'plaquette', rows)


def test_tail_convergence_does_not_accept_a_zero_crossing():
    records = [dict(contract={'r': r}, result={'doublons': [D]})
               for r, D in [(16, .10001), (32, .1)]]
    assert tail_spread(records, .1, 0) == 1
    records[-1]['result']['doublons'][0] = .1000025
    assert tail_spread(records, .1, 0) < 1e-10


def test_measured_budget_adds_numerical_allowance_before_optimization():
    from budget_optimization import optimize, settings
    from discussion_model import assumptions
    error = .00059431702957663
    row = optimize(4, 1, 'nisq', assumptions(), .001,
                   {**settings(), 'explicit_step_candidates': [21]},
                   C_override=(error+1e-10)*21**2)
    assert row['budget_fractions']['trotter']*.001 == pytest.approx(error+1e-10, abs=1e-15)
    assert row['budget_fractions']['residual_bias'] == .1
