"""First-order coherent cancellation in the branch-averaged RUS channel (docs/RAW_RUS_NOTE.md)."""
import numpy as np
import pytest

from raw_rus_check import branch_sum


def channels(scale, count=40, seed=3):
    rng = np.random.default_rng(seed)
    a, b = rng.uniform(0, 1, count), rng.uniform(-1, 1, count)
    return [1 - scale*x + 1j*scale*y for x, y in zip(a, b)], a


def test_coherence_cancels_at_first_order_for_any_sequence():
    for scale in (1e-3, 1e-4):
        cs, a = channels(scale)
        C, _ = branch_sum(cs)
        # Im C is second order; dephasing keeps the first-order weights 2^-m.
        assert abs(C.imag) < 5*scale**2
        assert 1 - C.real == pytest.approx(scale*sum(2.**-m*x for m, x in enumerate(a)), rel=1e-2)


def test_cancellation_needs_conjugate_failure_and_fair_success():
    cs, _ = channels(1e-3)
    assert abs(branch_sum(cs, conjugate_failure=False)[0].imag) > 1e-4
    assert abs(branch_sum(cs, success=.55)[0].imag) > 10*abs(branch_sum(cs)[0].imag)
