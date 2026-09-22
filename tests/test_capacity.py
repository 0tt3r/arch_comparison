"""Boundary checks for discrete inversion of the resource estimates."""
import math
import pytest
from plot_capacity import maximum_tested_side


def test_capacity_includes_threshold_and_preserves_unknown_region():
    points=[{'L':2,'log10_Q_min':2.},{'L':4,'log10_Q_min':4.}]
    assert math.isnan(maximum_tested_side(points,99))
    assert maximum_tested_side(points,100)==2
    assert maximum_tested_side(points,9999)==2
    assert maximum_tested_side(points,10000)==4
    assert maximum_tested_side(points,1e12)==4  # largest tested, not a physical limit


def test_capacity_takes_envelope_for_nonmonotone_costs():
    points=[{'L':4,'log10_Q_min':5.},{'L':2,'log10_Q_min':2.},{'L':6,'log10_Q_min':4.}]
    assert maximum_tested_side(points,1e4)==6


def test_capacity_rejects_nonpositive_hardware_budget():
    with pytest.raises(ValueError):maximum_tested_side([],0)
