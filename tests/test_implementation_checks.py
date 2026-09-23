"""Independent channel and arithmetic checks for the bounded implementation study."""
import gzip
import itertools
import json
from math import cos, sin, atan, tan

import numpy as np
import pytest

from hwp_stage import and_gates, nearest_neighbor_route, schedule_patch_gates, stage_operations, weight_network
from nisq_estimator_study import depolarize, native_circuit, observables, probabilities
from project_paths import ROOT
from small_dynamics import HubbardSector, rectangle
from star_channel_study import preparation, corrected_trial, growth_error, raw_rus_details, rus_channel


def gate_vector(psi,name,qubits,angle=0.):
    psi=psi.copy();indices=np.arange(len(psi))
    if name in ('cx','cz','conditional_cz'):
        a,b=qubits;mask=(indices>>a)&1
        if name=='cx':return psi[indices^(mask<<b)]
        psi[mask.astype(bool)&((indices>>b)&1).astype(bool)]*=-1
        return psi
    q=qubits[0];zero=indices[(indices&(1<<q))==0];one=zero^(1<<q)
    if name=='h':matrix=np.array([[1,1],[1,-1]])/np.sqrt(2)
    elif name=='rz':matrix=np.diag(np.exp(np.array([-1,1])*.5j*angle))
    else:matrix=np.diag([1,{'t':np.exp(.25j*np.pi),'tdg':np.exp(-.25j*np.pi),'s':1j,'sdg':-1j}[name]])
    psi[zero],psi[one]=matrix@np.array([psi[zero],psi[one]])
    return psi


def measure_reset(psi,q,outcome):
    indices=np.arange(len(psi));zero=indices[(indices&(1<<q))==0];one=zero^(1<<q)
    out=np.zeros_like(psi);out[zero]=(psi[zero]+(-1)**outcome*psi[one])/np.sqrt(2)
    return out


def test_temporary_and_is_exact_on_the_entire_clean_target_subspace():
    for x in range(4):
        psi=np.eye(8,dtype=complex)[:,x]
        for op in and_gates(0,1,2):psi=gate_vector(psi,op.name,op.qubits)
        expected=np.eye(8)[:,x+(4 if x==3 else 0)]
        assert np.max(abs(psi-expected))<1e-14


@pytest.mark.parametrize('size',[2,3,4,7,8,16])
def test_adder_network_all_binary_inputs(size):
    forward,bits,reverse=weight_network(list(range(size)),list(range(size,2*size)))
    for value in range(2**size):
        data=[(value>>k)&1 for k in range(2*size)]
        for op in forward:
            if op.name=='and':a,b,t=op.qubits;data[t]^=data[a]&data[b]
            else:a,b=op.qubits;data[b]^=data[a]
        assert sum(data[q]*2**k for k,q in bits)==value.bit_count()
        for op in reverse:
            if op.name=='unand':a,b,t=op.qubits;data[t]^=data[a]&data[b]
            else:a,b=op.qubits;data[b]^=data[a]
        assert data==[(value>>k)&1 for k in range(2*size)]


@pytest.mark.parametrize('size',[2,3,4])
def test_hwp_stage_preserves_arbitrary_superpositions_in_every_measurement_branch(size):
    theta=.137;ops,width=stage_operations(size,size,theta)
    rng=np.random.default_rng(917+size);data=rng.normal(size=4**size)+1j*rng.normal(size=4**size);data/=np.linalg.norm(data)
    psi=np.zeros(2**width,complex);psi[:len(data)]=data;branches=[psi]
    for op in ops:
        if op.name=='unand':
            a,b,t=op.qubits;new=[]
            for state in branches:
                for outcome in (0,1):
                    branch=measure_reset(state,t,outcome)
                    new.append(gate_vector(branch,'cz',(a,b)) if outcome else branch)
            branches=new
        else:branches=[gate_vector(v,op.name,op.qubits,op.angle) for v in branches]
    ideal=psi.copy()
    for x in range(len(data)):
        ideal[x]*=np.exp(-1j*theta*sum(1-2*(((x>>(2*i))^(x>>(2*i+1)))&1) for i in range(size)))
    assert sum(np.linalg.norm(v)**2 for v in branches)==pytest.approx(1,abs=1e-13)
    assert sum(abs(np.vdot(ideal,v))**2 for v in branches)==pytest.approx(1,abs=1e-13)


def test_native_circuit_keeps_training_noise_locations():
    a,_=native_circuit(.7,1,0);b,_=native_circuit(.7,1,4)
    assert [(g.name,g.qubits) for g in a]==[(g.name,g.qubits) for g in b]
    assert any(g.role=='onsite' and g.angle==0 for g in a)
    p,mapping=native_circuit(.7,1,4);probs,_=probabilities(p,0)
    exact=HubbardSector(rectangle(2,2));reference=exact.doublon(exact.trotter(.7,1),(1,1))
    assert observables(probs,mapping)['raw']['mean']==pytest.approx(reference,abs=2e-13)


@pytest.mark.parametrize('qubits',[(0,),(1,),(0,2)])
def test_depolarizing_channel_against_explicit_pauli_mixture(qubits):
    rng=np.random.default_rng(251);psi=rng.normal(size=8)+1j*rng.normal(size=8);psi/=np.linalg.norm(psi)
    rho=np.outer(psi,psi.conj());actual=rho.copy();p=.037
    masks=np.array([sum(1<<q for k,q in enumerate(qubits) if (s>>k)&1) for s in range(2**len(qubits))])
    depolarize(actual,masks,p,np.empty_like(actual))
    paulis=[np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    target=(1-p)*rho
    for labels in itertools.product(range(4),repeat=len(qubits)):
        if not any(labels):continue
        mapping=dict(zip(qubits,labels));op=np.array([[1.]])
        for q in reversed(range(3)):op=np.kron(op,paulis[mapping.get(q,0)])
        target+=p/(4**len(qubits)-1)*(op@rho@op.conj().T)
    assert np.max(abs(target-actual))<2e-15


@pytest.mark.parametrize('distance',[3,5,9,25])
def test_projection_and_pcec_limits(distance):
    theta=.024;k=(distance+1)//2;alpha=atan(tan(theta)**(1/k))
    ideal=preparation(theta,distance,0)
    assert ideal['acceptance']==pytest.approx(cos(alpha)**(2*k)+sin(alpha)**(2*k))
    assert ideal['multiplier']==pytest.approx(1)
    noisy=preparation(theta,distance,.001)
    corrected,*_=corrected_trial(theta,distance,.001,'ideal_pcec')
    assert corrected==pytest.approx(abs(noisy['multiplier'])**2,abs=2e-15)
    assert abs(corrected.imag)<1e-15


def test_rus_resolves_dyadic_angles_and_charges_angle_dependent_corrections():
    theta=np.pi/32
    c,rounds,attempts,tail=raw_rus_details(theta,5,0)
    assert c==pytest.approx(1) and tail==0 and attempts==pytest.approx(1.75)
    a=rus_channel(.013,5,.001,'noisy_pcec',floor=0)
    b=rus_channel(.013,5,.001,'noisy_pcec',floor=1e-5)
    assert a['mean_rounds']>rus_channel(.013,5,.001,'raw')['mean_rounds']
    assert a['extra_correction_attempts']>0 and b['q_dephasing']>a['q_dephasing']
    assert abs(a['coherent_component'])>0 # imperfect PCEC is not silently twirled


@pytest.mark.parametrize('distance',[3,5,9,15])
def test_growth_scenario_improves_with_p_at_every_fixed_distance(distance):
    a=growth_error(distance,.001,1e-5,25)
    b=growth_error(distance,.0001,1e-5,25)
    assert b/a==pytest.approx(.1**((distance+1)/2))


def test_routing_is_nearest_neighbor_and_feedforward_is_scheduled():
    ops,width=stage_operations(4,4,.13);gates,layout=nearest_neighbor_route(ops,width)
    for g in gates:
        if len(g['qubits'])==2:
            a,b=g['qubits'];s=layout['side']
            assert abs(a//s-b//s)+abs(a%s-b%s)==1
    angles={g['angle'].hex():dict(T=7,clifford=3) for g in gates if g['name']=='rz'}
    schedule=schedule_patch_gates(gates,layout['width'],angles,2)
    ready=[0]*layout['width']
    for event in schedule['events']:
        assert event['start']>=max(ready[q] for q in event['qubits'])
        if event['dependency'] is not None:
            assert event['start']>=schedule['events'][event['dependency']]['end']+1
        for q in event['qubits']:ready[q]=event['end']
    assert schedule['required_buffer_upper_estimate']<=2
    assert schedule['T_consumed']==4*(4-(4).bit_count())+3*7


def test_bounded_buffer_pauses_production_during_long_clifford_work():
    gates=[dict(name='t',qubits=[0],angle=0.,dependency=None)]
    gates += [dict(name='h',qubits=[0],angle=0.,dependency=None)]*30
    gates += [dict(name='t',qubits=[0],angle=0.,dependency=None)]*10
    a=schedule_patch_gates(gates,1,{},1,buffer_policy='bounded')
    b=schedule_patch_gates(gates,1,{},1,buffer_policy='continuous')
    assert a['logical_clocks']>b['logical_clocks'] and a['required_buffer_upper_estimate']==1
    assert b['required_buffer_upper_estimate']>1


@pytest.mark.parametrize('distance',[3,5,9])
def test_projected_channel_against_explicit_product_state_and_fault_patterns(distance):
    theta=.047;p=.013;k=(distance+1)//2;alpha=atan(tan(theta)**(1/k))
    # Independent repetition-subspace projection of product rotations, with
    # the known Clifford phase correction on the logical |1> component.
    rotation=np.array([[cos(alpha),1j*sin(alpha)],[1j*sin(alpha),cos(alpha)]])
    state=np.array([1.+0j])
    for _ in range(k):state=np.kron(state,rotation[:,0])
    density=np.zeros((2,2),complex);probabilities=[p/15]*(k-1)+[p/30]
    for pattern in range(2**k):
        probability=np.prod([q if (pattern>>i)&1 else 1-q for i,q in enumerate(probabilities)])
        selected=np.array([state[pattern],state[(2**k-1)^pattern]*1j**(1-k)])
        density+=probability*np.outer(selected,selected.conj())
    acceptance=np.trace(density).real;density/=acceptance
    c=np.exp(2j*theta)*(density[0,0]-density[1,1]-2j*density[1,0].imag)
    expected=preparation(theta,distance,p)
    assert acceptance==pytest.approx(expected['acceptance'],abs=2e-15)
    assert c==pytest.approx(expected['multiplier'],abs=2e-15)


def test_routed_hwp_stage_matches_unrouted_target_on_all_measurement_branches():
    sites=3;theta=.197;ops,width=stage_operations(sites,3,theta)
    gates,layout=nearest_neighbor_route(ops,width);n=layout['width']
    rng=np.random.default_rng(735);data=rng.normal(size=64)+1j*rng.normal(size=64);data/=np.linalg.norm(data)
    state=np.zeros(2**n,complex);state[:64]=data;branches=[(state,{})]
    for i,g in enumerate(gates):
        new=[]
        for psi,outcomes in branches:
            if g['name']=='mx_reset':
                for outcome in (0,1):new.append((measure_reset(psi,g['qubits'][0],outcome),outcomes|{i:outcome}))
            elif g['name']=='conditional_cz':
                new.append((gate_vector(psi,'cz',g['qubits']) if outcomes[g['dependency']] else psi,outcomes))
            else:new.append((gate_vector(psi,g['name'],g['qubits'],g['angle']),outcomes))
        branches=new
    ideal=np.zeros(2**n,complex)
    for value in range(64):
        location=sum(((value>>q)&1)<<layout['final_positions'][q] for q in range(6))
        ideal[location]=data[value]*np.exp(-1j*theta*sum(1-2*(((value>>(2*j))^(value>>(2*j+1)))&1) for j in range(3)))
    assert sum(abs(np.vdot(ideal,psi))**2 for psi,_ in branches)==pytest.approx(1,abs=2e-13)
