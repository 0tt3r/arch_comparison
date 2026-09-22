"""Independent checks of parallel fermionic routing and sector dynamics."""

from project_paths import ROOT, project_file
import numpy as np
import pytest
from scipy.linalg import expm

from hubbard_circuits import Operation, build_workload, fermionic_hadamard, parity_rotation
from nisq_compilation import Block, cudaq_operations, fermionic_blocks, schedule
from parallel_nisq import Stage, observable_stages, parallel_route, plaquette_stages
from small_dynamics import HubbardSector, rectangle
from test_hubbard_circuits import annihilation, operator, pauli_matrix, PAULI
from test_nisq_refinement import equal_up_to_phase, fock_permutation


@pytest.mark.parametrize('L,r',[(2,1),(4,4),(6,3)])
def test_parallel_stages_are_original_ordered_workload(L,r):
    _,expected=fermionic_blocks(build_workload(L,1.,r)[0])
    assert expected==[b for s in plaquette_stages(L,1.,r) for group in s.groups for b in group]


@pytest.mark.parametrize('order',[(0,1,2,3),(0,2,1,3)])
@pytest.mark.parametrize('inside',[False,True])
def test_parallel_regrouping_unitary_with_fermionic_signs(order,inside):
    stages=[Stage('first',[[Block('hop',(0,3),.23)],[Block('hop',(1,2),-.19)]]),
            Stage('second',[[Block('density',(0,1),.27)],[Block('density',(2,3),-.31)]]),
            Stage('third',[[Block('hadamard',(3,1))],[Block('hop',(0,2),.17)]])]
    ops=[]
    for stage in stages:
        for group in stage.groups:
            for b in group:
                if b.kind=='hop':ops.extend(parity_rotation(*b.modes,p,b.theta) for p in ('XX','YY'))
                elif b.kind=='density':ops.append(Operation('rotate',b.modes,b.theta,'ZZ'))
                else:ops.extend(fermionic_hadamard(*b.modes))
    gates,final,_=parallel_route(stages,order,[],optimize_inside=inside)
    initial={m:q for q,m in enumerate(order)}
    desired=fock_permutation(final)@operator(ops,4)@fock_permutation(initial).T
    assert equal_up_to_phase(operator(cudaq_operations(gates),4),desired)
    schedule(gates,4)


def evolve_qubit_state(ops,n):
    psi=np.zeros(2**n,dtype=complex);psi[0]=1
    for op in ops:
        if op.name=='rotate':
            P=pauli_matrix(n,op.qubits,op.paulis)
            psi=np.cos(op.angle)*psi-1j*np.sin(op.angle)*(P@psi)
        else:
            # Only a single-qubit gate here; reuse the independently tested
            # dense helper as a tensor factor, not the compiler's decomposition.
            local=operator([Operation(op.name,(0,),op.angle)],1)
            q=op.qubits[0];shape=psi.reshape([2]*n)
            shape=np.moveaxis(shape,q,0);shape=np.tensordot(local,shape,axes=(1,0))
            psi=np.moveaxis(shape,0,q).reshape(-1)
    return psi


def measured_doublon(psi,n,wires):
    return sum(abs(psi[i])**2 for i in range(2**n) if all(i&(1<<(n-1-q)) for q in wires))


def test_observable_pruning_and_parallel_compile_against_original_qubit_circuit():
    ops,observed=build_workload(2,1.,4)
    prep,_=fermionic_blocks(ops)
    stages,modes=observable_stages(plaquette_stages(2,1.,4),observed)
    order=sorted(modes,key=lambda m:(m%2,m))
    native,final,_=parallel_route(stages,order,prep,optimize_inside=True)
    desired=evolve_qubit_state(ops,8)
    actual=evolve_qubit_state(cudaq_operations(native),8)
    assert np.isclose(measured_doublon(desired,8,observed),
                      measured_doublon(actual,8,[final[m] for m in observed]),atol=2e-12)


def test_sector_trotter_matches_original_eight_qubit_workload():
    model=HubbardSector(rectangle(2,2))
    ops,observed=build_workload(2,1.,4)
    psi=evolve_qubit_state(ops,8)
    assert np.isclose(model.doublon(model.trotter(1.,4),(1,1)),
                      measured_doublon(psi,8,observed),atol=2e-12)


def test_exact_sector_dynamics_against_full_qubit_hamiltonian():
    sites=rectangle(2,2);modes={(x,y,s):2*i+s for i,(x,y) in enumerate(sites) for s in (0,1)}
    a=[annihilation(8,i) for i in range(8)]
    H=np.zeros((256,256),dtype=complex)
    for x,y in sites:
        H+=pauli_matrix(8,(modes[x,y,0],modes[x,y,1]),'ZZ')
        for xx,yy in ((x+1,y),(x,y+1)):
            if (xx,yy) in sites:
                for spin in (0,1):
                    i,j=modes[x,y,spin],modes[xx,yy,spin]
                    H-=a[i].conj().T@a[j]+a[j].conj().T@a[i]
    ops,observed=build_workload(2,.37,1)
    initial=evolve_qubit_state([op for op in ops if op.name=='x'],8)
    psi=expm(-.37j*H)@initial
    model=HubbardSector(sites)
    assert np.isclose(model.doublon(model.exact(.37),(1,1)),measured_doublon(psi,8,observed),atol=2e-12)
    assert np.isclose(np.linalg.norm(model.exact(.37)),1.)


def test_second_order_asymptotic_local_error():
    model=HubbardSector(rectangle(2,2));exact=model.doublon(model.exact(1),(1,1))
    errors=[abs(model.doublon(model.trotter(1,r),(1,1))-exact) for r in (8,16,32)]
    assert 3.9<errors[0]/errors[1]<4.1
    assert 3.9<errors[1]/errors[2]<4.1


def test_star_instrument_covariance_for_all_two_qubit_pauli_faults():
    from itertools import product
    from star_parallel import I,Z,injection_kraus
    theta=.137
    for where in ('before','after'):
        for fault in product(range(4),repeat=2):
            for outcome in (0,1):
                K=injection_kraus(theta,outcome,**{where:fault})
                assert np.allclose(K.conj().T@K,I/2)
                assert np.allclose(K@Z,Z@K) or np.allclose(K@Z,-Z@K)
    success=injection_kraus(theta,0)
    assert np.allclose(success,np.diag(np.exp(-1j*theta*np.array([1,-1])))/np.sqrt(2))
    corrected=injection_kraus(2*theta,0)@injection_kraus(theta,1)
    assert np.allclose(corrected,success/np.sqrt(2))


def test_parallel_rus_tail_cost_and_dependency_support():
    from star_parallel import expected_max_attempts,rotation_batches
    assert np.isclose(expected_max_attempts(1),2.)
    assert 2<expected_max_attempts(4)<8
    ops=[Operation('rotate',(0,1),.13,'ZZ'),Operation('rotate',(2,3),.13,'ZZ'),
         Operation('rotate',(1,2),.13,'ZZ')]
    assert rotation_batches(ops,4,2)==([2,1],0)
    assert rotation_batches(ops,4,1)==([1,1,1],0)


def test_cultivation_checks_state_budget_and_handoff_error():
    import json
    from pathlib import Path
    from cultivation_model import cultivated_encoded
    c=json.loads(project_file('assumptions.json').read_text())
    source=json.loads(project_file('cultivation_assumptions.json').read_text())
    record=json.loads((ROOT / 'results/cudaq_v1/plaquette_L4_tau1_r4/result.json').read_text())
    ideal=cultivated_encoded(record,c,'ft',source)
    handoff=cultivated_encoded(record,c,'ft',source,handoff_error=1e-8)
    assert ideal['post_cultivation_distillation_levels']==0
    assert handoff['post_cultivation_distillation_levels']==1
    assert handoff['magic_output_error']<ideal['magic_output_error']
    assert cultivated_encoded(record,c,'hybrid',source,handoff_error=1e-8)['post_cultivation_distillation_levels']==0


@pytest.mark.parametrize('qubits',[(1,),(0,2)])
def test_density_depolarizing_channel_against_explicit_pauli_mixture(qubits):
    from itertools import product
    from noisy_small_observable import DensityNoise
    rng=np.random.default_rng(9127)
    matrix=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8))
    rho=matrix@matrix.conj().T;rho/=np.trace(rho)
    p=.017;expected=(1-p)*rho
    for labels in product('IXYZ',repeat=len(qubits)):
        if set(labels)=={'I'}:continue
        P=pauli_matrix(3,qubits,''.join(labels))
        expected+=p/(4**len(qubits)-1)*(P@rho@P)
    assert np.allclose(DensityNoise(3).depolarize(rho,qubits,p),expected,atol=2e-15)


def test_density_native_gates_against_independent_dense_operators():
    from noisy_small_observable import DensityNoise
    rng=np.random.default_rng(123)
    psi=rng.normal(size=8)+1j*rng.normal(size=8);psi/=np.linalg.norm(psi)
    rho=np.outer(psi,psi.conj())
    gates=[Operation('x',(1,)),Operation('h',(0,)),Operation('rx',(2,),.27),
           Operation('rz',(0,),-.38),Operation('rzz',(0,2),.17)]
    for gate in gates:
        U=operator(cudaq_operations([gate]),3)
        assert np.allclose(DensityNoise(3).gate(rho,gate),U@rho@U.conj().T,atol=2e-15)
