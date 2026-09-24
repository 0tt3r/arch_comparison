"""Factory re-pricing, RUS sign reconciliation, and the strip solver copy."""
import json

import pytest

import factory_calibration as fc
import rus_reconciliation as rr
from project_paths import ROOT


def test_vendored_factory_model_reproduces_published_row_and_uses_calibrated_law():
    ref = fc.evaluate(('reference', 'one', 1e-4, (7, 3, 3)))
    cal = fc.evaluate(('calibrated', 'one', 1e-4, (7, 3, 3)))
    assert ref['error'] == pytest.approx(4.4e-8, rel=.02) and ref['qubits'] == 810
    assert cal['error'] > ref['error'] and cal['qubits'] == ref['qubits']


def test_factory_menu_is_pareto_by_decade():
    menu = json.loads((fc.OUT/'menu.json').read_text())['litinski_factories']
    for rows in menu.values():
        assert all(0 < r[1] < 1e-6 and r[2] > 0 and r[3] > 0 for r in rows)


@pytest.mark.parametrize('k', [2, 3, 4, 6])
def test_conjugate_failures_cancel_coherence_and_halve_pcec(k):
    for theta in (1e-3, .05, .3):
        x, y, pcec = rr.averaged(theta, k, 1e-3, 'conjugate')
        assert abs(y) < 1e-12*x and x == pytest.approx(pcec/2, rel=1e-9)
        xs, ys, _ = rr.averaged(theta, k, 1e-3, 'same_sign')
        assert abs(ys) > 1e3*abs(y)


def test_strip_solver_differs_from_pinned_solver_only_in_site_limit():
    a = (ROOT/'cpp/hubbard_sector.cpp').read_text().splitlines()
    b = (ROOT/'cpp/hubbard_strip.cpp').read_text().splitlines()
    changed = [line for line in b if line not in a]
    assert len(changed) == 3 and any('n>18' in line for line in changed)
