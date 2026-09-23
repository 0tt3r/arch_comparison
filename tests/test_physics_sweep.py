"""Independent invariants for source scenarios and repeated-angle compilation."""
import cmath
import math

import numpy as np
import pytest

from budget_optimization import optimize
from discussion_model import assumptions, estimate
from local_ordering import OrderedSector
from model import log_gamma2_dephasing, workload
from physics_sweep import configuration, planning_coefficient
from resource_extensions import (analog_exposure, compact_physical_counts,
                                 hamming_weight_counts, plaquette_rotation_stages,
                                 source_error, synthesis_plan)


def test_unanchored_source_requires_an_explicit_scenario():
    c = {**assumptions(), 'p':1e-4}
    with pytest.raises(ValueError, match='anchored'):source_error(c)
    assert source_error(configuration(1e-4,'held')) == 2e-9
    assert source_error(configuration(1e-4,'scaled')) == pytest.approx(2e-12)
    for policy in ('held_at_anchor','power3'):
        with pytest.raises(ValueError, match='anchor'):source_error({**c,'p':.002,'source_p_policy':policy})


def test_count_and_angle_injection_have_distinct_step_scaling():
    c={**assumptions(),'star_floor':0}
    a,b=[workload(4,1,r,'termwise') for r in (8,16)]
    assert analog_exposure(a,c)==analog_exposure(b,c)
    c['star_law']='original_count'
    x,y=analog_exposure(a,c)[0],analog_exposure(b,c)[0]
    assert y/x==pytest.approx(b.rotations/a.rotations)
    assert x==pytest.approx(a.rotations*log_gamma2_dephasing(4*c['p']/15))


def test_count_law_optimizer_keeps_the_exposure_ledger_consistent():
    c={**configuration(.001),'star_law':'original_count'}
    v=optimize(4,1,'star',c,.01,C_override=planning_coefficient(4,1,'termwise'))
    assert v['chi_angle']==v['chi_floor']==0
    assert v['chi_count']==pytest.approx(v['chi'])
    assert v['chi_count']==pytest.approx(v['rotations']*log_gamma2_dephasing(4*c['p']/15))


@pytest.mark.parametrize('size',[1,2,3,4,7,8,16])
def test_hamming_weight_phases_match_all_weight_sectors(size):
    # Independent action of the weight-register rotations, up to global phase.
    theta=.173;bits=size.bit_length()
    phase_offset=cmath.exp(-1j*theta*(size-(2**bits-1)))
    for weight in range(size+1):
        ideal=cmath.exp(-1j*theta*(size-2*weight))
        actual=phase_offset
        for bit in range(bits):
            actual*=cmath.exp(-1j*theta*2**bit*(1-2*((weight>>bit)&1)))
        assert abs(actual-ideal)<2e-14
    # Number of carries in a binary adder tree, independently summed by level.
    assert sum(size//2**k for k in range(1,bits+1))==size-size.bit_count()


@pytest.mark.parametrize('L',[2,4,6,8])
def test_hwp_does_not_pool_rotations_from_different_stages(L):
    w=workload(L,1,7,'plaquette')
    assert sum(width*repeat for width,repeat in plaquette_rotation_stages(w))==w.rotations
    identity=hamming_weight_counts(w,1)
    assert identity['synthesized_rotations']==w.rotations
    assert identity['hwp_toffolis']==identity['hwp_workspace']==0
    c=assumptions()
    assert synthesis_plan(w,.001,c)['T_count']==synthesis_plan(w,.001,{**c,'rotation_synthesis':'hamming_weight','hwp_batch_size':1})['T_count']


def test_hwp_workspace_and_synthesis_error_are_both_charged():
    c={**configuration(.001),'rotation_synthesis':'hamming_weight','hwp_batch_size':16}
    v=optimize(4,1,'ft',c,.001,C_override=planning_coefficient(4,1,'plaquette'))
    assert v['hwp_workspace']==15
    assert v['data_tiles']>=math.ceil(1.5*(32+6+15))+24
    w=workload(4,1,v['r'],'plaquette');b=v['budget_fractions']
    plan=synthesis_plan(w,.001,{**c,'synthesis_fraction':b['synthesis']})
    assert plan['T_count']==v['T_count']
    assert v['T_count']==w.fixed_T+4*v['hwp_toffolis']+v['synthesized_rotations']*v['T_per_rotation']
    eta=.001*b['synthesis']/(2*v['synthesized_rotations'])
    assert 2*v['synthesized_rotations']*eta==pytest.approx(.001*b['synthesis'])


def test_compact_counts_include_preparation_and_feasible_busy_time():
    c=configuration(.001);v=compact_physical_counts(4,11,c)
    # 24 bonds, four Paulis/bond, two Strang half steps.
    hopping=2*11*4*24;onsite=2*11*(16-1)+1
    assert v['two_qubit_gates']==4*hopping+2*onsite+8*16
    assert v['qubits_per_copy']>=3*16+2
    assert v['depth']*(3*16+2)>=2*v['two_qubit_gates']
    assert compact_physical_counts(4,11,{**c,'local_routing_multiplier':2})['two_qubit_gates']>v['two_qubit_gates']


def test_colored_order_covers_bonds_and_remains_second_order():
    model=OrderedSector(2,2)
    edges=[e for group in model.colored_edges() for e in group]
    assert set(edges)==set(model.edges) and len(edges)==len(model.edges)
    exact=model.doublon(model.exact(.7),(1,1))
    errors=[]
    for r in (16,32):
        psi=model.trotter(.7,r,'colored')
        assert abs(np.linalg.norm(psi)-1)<1e-12
        errors.append(abs(model.doublon(psi,(1,1))-exact))
    assert 3.8<errors[0]/errors[1]<4.2


def test_hwp_batch_one_replays_the_frozen_resource_model():
    c=configuration(.001)
    a=estimate(4,1,'ft',c,.01,steps=11)
    b=estimate(4,1,'ft',{**c,'rotation_synthesis':'hamming_weight','hwp_batch_size':1},.01,steps=11)
    for key in ('T_count','distance','seconds_per_attempt','log10_Q_min'):
        assert a[key]==b[key]


def test_frozen_source_archive_cannot_mask_changed_scientific_inputs(tmp_path,monkeypatch):
    import hashlib
    import zipfile
    import frozen_baseline
    monkeypatch.setattr(frozen_baseline,'ROOT',tmp_path)
    original=b'original source';digest=hashlib.sha256(original).hexdigest()
    (tmp_path/'python').mkdir();(tmp_path/'python/model.py').write_bytes(b'new source')
    folder=tmp_path/'results/physics_sweep_v1';folder.mkdir(parents=True)
    (tmp_path/'results/data.json').write_bytes(b'changed data')
    with zipfile.ZipFile(folder/'baseline_sources.zip','w') as archive:
        archive.writestr('python/model.py',original)
        archive.writestr('results/data.json',original)
    frozen_baseline.verify_sources({'python/model.py':digest})
    with pytest.raises(AssertionError):
        frozen_baseline.verify_sources({'python/model.py':'0'*64})
    with pytest.raises(AssertionError):
        frozen_baseline.verify_sources({'results/data.json':digest})
