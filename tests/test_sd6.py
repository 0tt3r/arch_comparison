"""Uniform-noise (SD6) additions: the physical PEC reference and the round plumbing."""
import pytest

import optimization_round as o
from model import log_gamma2_depolarizing


def test_pec_uniform_adds_one_qubit_exposure_at_p():
    c = {**o.sweep_configuration(1e-3), 'nisq_model': 'pec_uniform', 'nisq_one_qubit_ratio': .8}
    two = o._ESTIMATE(6, 1.5, 'nisq', {**c, 'nisq_model': 'pec_two_qubit'}, .01, steps=12)
    uni = o.estimate_uniform(6, 1.5, 'nisq', c, .01, steps=12)
    extra = .8*two['two_qubit_gates']*log_gamma2_depolarizing(1e-3, 1)
    assert uni['chi'] == pytest.approx(two['chi'] + extra, rel=1e-12)
    assert uni['one_qubit_gates'] == pytest.approx(.8*two['two_qubit_gates'])


def test_rounds_select_growth_source_and_restore_it():
    o.use_round('v5')
    try:
        assert o.v1.growth_pairs is not o._GROWTH_PAIRS
        assert o.v1.growth_probability(5, 1e-3, 'component', 25) > 0
    finally:
        o.use_round('v1')
    assert o.v1.growth_pairs is o._GROWTH_PAIRS
