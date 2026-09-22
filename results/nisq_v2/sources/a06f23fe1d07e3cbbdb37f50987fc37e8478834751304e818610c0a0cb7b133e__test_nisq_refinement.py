"""Small gate/channel identities, not a Hubbard dynamics benchmark."""
from math import pi
import json
from pathlib import Path

import numpy as np
import pytest

from hubbard_circuits import Operation, build_workload, fermionic_hadamard, parity_rotation
from nisq_compilation import (Block, causal_cone, cudaq_operations, fermionic_blocks,
                              fused_route, gaussian_rzz, grid_parity_route, schedule)
from nisq_observables import observable_exposure, parity_acceptance
from test_hubbard_circuits import operator


def equal_up_to_phase(a,b):
    phase=np.angle(np.trace(b.conj().T@a))
    return np.allclose(a*np.exp(-1j*phase),b,atol=2e-12)


@pytest.mark.parametrize('seed',range(8))
def test_gaussian_native_decomposition(seed):
    rng=np.random.default_rng(seed)
    matrix,_=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
    desired=np.eye(4,dtype=complex)
    desired[1:3,1:3]=matrix
    desired[3,3]=np.linalg.det(matrix)
    gates=gaussian_rzz((0,1),matrix)
    assert sum(g.name=='rzz' for g in gates)<=2
    assert equal_up_to_phase(operator(cudaq_operations(gates),2),desired)


def fock_permutation(position):
    n=len(position);result=np.zeros((2**n,2**n))
    for basis in range(2**n):
        occupied=[i for i in range(n) if basis&(1<<(n-1-i))]
        out=[position[i] for i in occupied]
        sign=(-1)**sum(out[i]>out[j] for i in range(len(out)) for j in range(i+1,len(out)))
        row=sum(1<<(n-1-i) for i in out)
        result[row,basis]=sign
    return result


def test_fused_route_preserves_full_four_mode_unitary():
    ops=[]
    for i,j,theta in [(0,3,.13),(1,2,-.23),(0,2,.41),(1,3,-.14)]:
        ops.extend([parity_rotation(i,j,'XX',theta),parity_rotation(i,j,'YY',theta)])
        ops.extend(fermionic_hadamard(j,i))
        ops.append(Operation('rotate',(0,1),.17,'ZZ','dynamic'))
    _,blocks=fermionic_blocks(ops)
    gates,position,info=fused_route(blocks,list(range(4)),[])
    assert info['pair_fusions']>0
    desired=fock_permutation(position)@operator(ops,4)
    assert equal_up_to_phase(operator(cudaq_operations(gates),4),desired)


def test_grid_tree_preserves_nonlocal_pauli_rotations():
    ops=[Operation('rotate',(0,1,2,3),.17,'XZZY','dynamic'),
         Operation('rotate',(0,1,2,3),-.21,'YZZX','dynamic')]
    gates,info=grid_parity_route(ops,4)
    assert info['off_path_entanglers_before_local_simplification']>0
    assert equal_up_to_phase(operator(cudaq_operations(gates),4),operator(ops,4))
    schedule(gates,4)


def test_fermionic_cone_with_occupied_spectator():
    # The occupied spectator lies inside the full JW string. Compressing modes
    # preserves diagonal local occupation probabilities, not arbitrary qubit amplitudes.
    blocks=[Block('hop',(0,2),.23),Block('hop',(1,3),.47)]
    selected,modes=causal_cone(blocks,[0])
    assert selected==blocks[:1] and modes==[0,2]
    full_ops=[parity_rotation(0,2,p,.23) for p in ('XX','YY')]
    full_ops +=[parity_rotation(1,3,p,.47) for p in ('XX','YY')]
    full=operator(full_ops,4)[:,12]  # |1100>, including occupied mode 1
    gates,position,_=fused_route(selected,modes,[0,1])
    reduced=operator(cudaq_operations(gates),2)[:,0]
    full_p=sum(abs(full[i])**2 for i in range(16) if i&8)
    reduced_p=sum(abs(reduced[i])**2 for i in range(4) if i&(1<<(1-position[0])))
    assert np.isclose(full_p,reduced_p)


@pytest.mark.parametrize('just_in_time',[False,True])
def test_parity_acceptance_against_independent_density_channel(just_in_time):
    # One same-spin hopping block, with the actual native gate/basis noise.
    gates,position,_=fused_route([Block('hop',(0,2),.217)],[0,2],[0])
    gates=(gates[:1]+[Operation('h',(0,)),Operation('h',(0,)),Operation('rzz',(0,1),.17)]
           +gates[1:]+[Operation('h',(0,)),Operation('h',(0,))])
    c=json.loads(Path('assumptions.json').read_text());c['p']=.017
    sched,events=schedule(gates,2)
    post=parity_acceptance(events,sched,position,c)
    got=post['acceptance']
    from boundary_timing import boundary_gaps, refine
    head,tail=boundary_gaps(gates,sched,events)
    if just_in_time:
        record={'L':2,'tau':1.,'steps':1,'algorithm':'plaquette','configuration':c,
                'schedule':sched,'mode_to_wire':position,'observable_wires':[0],
                'observable_exposure':observable_exposure(events,sched,[0],c),'postselection':post}
        result=refine(record,gates)
        assert result['removed_idle_qubit_layers']>0
        got=result['postselection_acceptance']
    rho=np.diag([1.,0,0,0]).astype(complex)
    I=np.eye(2);X=np.array([[0,1],[1,0]]);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1,-1])
    def channel(rho,qubits,p):
        from itertools import product
        terms=[]
        for choices in product((I,X,Y,Z),repeat=len(qubits)):
            factors=[I,I]
            for q,pauli in zip(qubits,choices):factors[q]=pauli
            terms.append(np.kron(*factors))
        return (1-p)*rho+p/(len(terms)-1)*sum(P@rho@P.conj().T for P in terms[1:])
    def flip(rho,q,p):
        P=np.kron(X,I) if q==0 else np.kron(I,X)
        return (1-p)*rho+p*P@rho@P
    for q in range(2):rho=flip(rho,q,c['p'])
    by_layer={}
    for start,gate in events:by_layer.setdefault(start,[]).append(gate)
    from nisq_observables import idle_probability
    for layer in range(sched['depth_layers']):
        busy=set()
        for gate in by_layer[layer]:
            U=operator(cudaq_operations([gate]),2);rho=U@rho@U.conj().T
            rho=channel(rho,gate.qubits,c['p'] if len(gate.qubits)==2 else c['p']/10)
            busy.update(gate.qubits)
        for q in set(range(2))-busy:
            if not just_in_time or head[q]<=layer<sched['depth_layers']-tail[q]:
                rho=channel(rho,(q,),idle_probability(c))
    for q in range(2):rho=flip(rho,q,c['p'])
    expected=(rho[1,1]+rho[2,2]).real  # one up-spin particle; no down-spin wires
    assert np.isclose(got,expected,atol=2e-12)


@pytest.mark.parametrize('algorithm',['plaquette','termwise'])
def test_native_workload_parities_remain_pauli(algorithm):
    c=json.loads(Path('assumptions.json').read_text())
    ops,_=build_workload(4,1,4,algorithm)
    prep,blocks=fermionic_blocks(ops)
    gates,position,_=fused_route(blocks,list(range(32)),prep)
    sched,events=schedule(gates,32)
    result=parity_acceptance(events,sched,position,c)
    assert .25<=result['acceptance']<=1


def test_magic_injection_error_changes_factory_without_changing_physical_error():
    from compiled_model import encoded
    c=json.loads(Path('assumptions.json').read_text())
    record=json.loads(Path('results/cudaq_v1/plaquette_L4_tau1_r4/result.json').read_text())
    before=encoded(record,c,'hybrid',injection_error=c['p'])
    after=encoded(record,c,'hybrid',injection_error=46*c['p']/15)
    assert before['factory_level']==1 and after['factory_level']==2
    assert before['log10_gamma2']==after['log10_gamma2']
    assert c['p']==.001
    loose={**c,'relative_error':.1}
    assert encoded(record,loose,injection_error=c['p'])['factory_level']==1
    assert encoded(record,loose,injection_error=46*c['p']/15)['factory_level']==2
