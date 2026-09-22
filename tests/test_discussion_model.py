"""Small regression checks for the discussion estimates, not calibration."""
import math

from discussion_model import assumptions,capacity,estimate


def test_willow_anchor_and_integer_device_capacity():
    c=assumptions()
    row=estimate(6,.75,'nisq',c,.001)  # Three Trotter steps, the published point.
    assert row['two_qubit_gates']==4372
    assert row['depth']==294
    assert math.isclose(math.exp(row['chi']),100)
    assert row['qubits_per_copy']==81
    assert capacity([row],80)==0 and capacity([row],81)==6


def test_star_angle_law_is_step_independent_without_floor():
    c={**assumptions(),'star_floor':0}
    a=estimate(4,1,'star',c,.001)
    b=estimate(4,1,'star',{**c,'steps_per_tau':8},.001)
    assert math.isclose(a['chi'],.384)  # Theta = (5 L^2 - 4 L) tau = 64.
    assert a['chi']==b['chi']
    assert b['seconds_per_attempt']>a['seconds_per_attempt']


def test_ft_cleanup_and_actual_logical_budget():
    c=assumptions()
    strict=estimate(6,1,'ft',c,.001)
    loose=estimate(6,1,'ft',c,.01)
    assert strict['cleanup_levels']==1 and loose['cleanup_levels']==0
    for row,epsilon in [(strict,.001),(loose,.01)]:
        assert row['T_count']*row['magic_error']<=epsilon*c['magic_fraction']
        protected_qubits=row['qubits_per_copy']-row['factories']*c['source_qubits']
        tiles=protected_qubits/(2*row['distance']**2)
        cycles=math.ceil(row['seconds_per_attempt']/c['cycle_seconds'])
        faults=tiles*cycles*c['surface_A']*(c['p']/c['surface_threshold'])**((row['distance']+1)/2)
        assert faults<=epsilon*c['logical_fraction']
