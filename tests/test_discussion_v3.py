"""Physics-ordering checks and independent error-budget accounting."""
import math

import numpy as np
import pytest

from budget_optimization import coefficient,optimize,settings
from discussion_model import assumptions,estimate
from hubbard_circuits import build_workload
from small_dynamics import HubbardSector,rectangle
from time_calibration import SplitSector
from test_refinement_v3 import evolve_qubit_state,measured_doublon


@pytest.mark.parametrize('algorithm',['plaquette','termwise'])
def test_factorized_sector_matches_the_actual_qubit_ordering(algorithm):
    model=SplitSector(2,2)
    ops,observed=build_workload(2,.83,3,algorithm)
    qubit=evolve_qubit_state(ops,8)
    sector=model.trotter(.83,3,algorithm)
    assert abs(model.doublon(sector,(1,1))-measured_doublon(qubit,8,observed))<2e-12
    assert abs(np.linalg.norm(sector)-1)<2e-12


def test_factorized_odd_patch_matches_independent_sparse_plaquette():
    fast=SplitSector(3,3);sparse=HubbardSector(rectangle(3,3))
    assert np.allclose(fast.trotter(.8,12),sparse.trotter(.8,12),atol=2e-12)


def test_time_lookup_preserves_ordering_and_rejects_extrapolation():
    assert coefficient(1,'termwise')<coefficient(1,'plaquette')
    assert coefficient(8,'termwise')>coefficient(8,'plaquette')
    with pytest.raises(ValueError,match='No measured coefficient'):coefficient(9,'plaquette')


def test_finite_step_envelope_covers_loose_precision_correction():
    model=SplitSector(3,3)
    exact=model.doublon(model.exact(4),(1,1))
    error=abs(model.doublon(model.trotter(4,35,'termwise'),(1,1))-exact)
    # This point exceeded the asymptotic coefficient plus 10 percent.
    assert error>1.1*1.9568859593050547/35**2
    assert error<=coefficient(4,'termwise')/35**2


@pytest.mark.parametrize('arm',['nisq','star','hybrid','ft'])
def test_optimized_budget_closes_and_individual_faults_fit(arm):
    c=assumptions();eps=.001;row=optimize(4,1,arm,c,eps);b=row['budget_fractions']
    assert math.isclose(sum(b.values()),1,abs_tol=1e-12)
    assert min(b.values())>=0 and b['stat']>0
    assert row['assumed_trotter_error']<=b['trotter']*eps*(1+1e-12)
    if arm!='nisq':
        qprotected=row['qubits_per_copy']-row.get('factories',0)*c['source_qubits']
        tiles=qprotected/(2*row['distance']**2)
        faults=tiles*math.ceil(row['seconds_per_attempt']/c['cycle_seconds'])*c['surface_A']*(c['p']/c['surface_threshold'])**((row['distance']+1)/2)
        assert math.exp(row['chi']/2)*faults<=eps*b['logical']*(1+1e-12)
    if arm in ('ft','hybrid'):
        assert math.exp(row['chi']/2)*row['T_count']*row['magic_error']<=eps*b['magic']*(1+1e-12)
    if arm=='ft':
        assert row['cleanup_levels']==0
        assert row['log10_Q_min']<estimate(4,1,'ft',c,eps)['log10_Q_min']


def test_nisq_optimizer_reserves_unmeasured_bias():
    row=optimize(4,1,'nisq',assumptions(),.001)
    assert row['budget_fractions']['residual_bias']==.1
    assert row['budget_fractions']['magic']==row['budget_fractions']['synthesis']==0


def test_two_by_two_resource_point_uses_its_own_calibration():
    row=optimize(2,8,'star',assumptions(),.001)
    assert row['coefficient_patch']=='2x2'
    expected=coefficient(8,'termwise',{**settings(),'coefficient_patch':'2x2'})
    assert row['coefficient']==expected
