"""Small regression checks for the discussion estimates, not calibration."""
import math

from discussion_model import assumptions,capacity,estimate,fit_evidence,step_count


def test_willow_anchor_and_integer_device_capacity():
    c=assumptions()
    row=estimate(6,.75,'nisq',c,.001,steps=3)  # Published point, separate from the fit scenario.
    assert row['two_qubit_gates']==4372
    assert row['depth']==294
    assert math.isclose(math.exp(row['chi']),100)
    assert row['qubits_per_copy']==81
    assert capacity([row],80)==0 and capacity([row],81)==6


def test_star_angle_law_is_step_independent_without_floor():
    c={**assumptions(),'star_floor':0}
    a=estimate(4,1,'star',c,.001,steps=4)
    b=estimate(4,1,'star',c,.001,steps=8)
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


def test_trotter_budget_conventions_and_fit_provenance():
    c=assumptions()
    assert [step_count(1,e,c) for e in (.001,.01)]==[30,10]
    split={**c,'trotter_budget_fraction':.125}
    assert [step_count(1,e,split) for e in (.001,.01)]==[42,14]
    fit=fit_evidence(c)
    assert len(fit['points'])==6
    assert all(.217<v['coefficient']<.224 for v in fit['points'])
    assert round(fit['median_coefficient'],2)==c['trotter_coefficient']
    assert step_count(4,.001,c)==238


def test_compiled_pec_reproduces_saved_four_step_headline():
    c={**assumptions(),'nisq_model':'compiled_pec'}
    row=estimate(4,1,'nisq',c,.001,steps=4)
    assert row['shot_rule']=='hoeffding'
    assert math.isclose(row['chi'],20.017057134229688)
    assert round(10**row['log10_Q_min'])==2132604


def test_hybrid_supplies_only_fixed_states_and_budgets_amplification():
    c=assumptions()
    row=estimate(4,1,'hybrid',c,.001,steps=4)
    assert row['T_count']==row['fixed_T']==576
    assert row['analog_rotations']==288
    assert math.isclose(row['chi'],.27552)
    assert row['T_count']*row['magic_error']*math.exp(row['chi']/2)<=.001*c['magic_fraction']
