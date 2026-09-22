"""Local gate identities and compilation contracts, not dynamics calibration."""
from math import isclose, pi

import numpy as np
from scipy.linalg import expm
import pytest

from hubbard_circuits import (angle_summary, build_workload, fermionic_hadamard,
                              native_path, native_fermionic_route, parity_rotation,
                              plaquette_evolution, schedule_native, Operation)
from model import workload

PAULI = {"I": np.eye(2), "X": np.array([[0,1],[1,0]]),
         "Y": np.array([[0,-1j],[1j,0]]), "Z": np.diag([1,-1])}


def tensor(matrices):
    result = np.array([[1.]])
    for matrix in matrices:
        result = np.kron(result, matrix)
    return result


def pauli_matrix(n, qubits, letters):
    factors = [PAULI["I"] for _ in range(n)]
    for q, p in zip(qubits, letters):
        factors[q] = PAULI[p]
    return tensor(factors)


def operator(ops, n):
    result = np.eye(2**n, dtype=complex)
    for op in ops:
        if op.name == "rotate":
            p = pauli_matrix(n, op.qubits, op.paulis)
            gate = np.cos(op.angle)*np.eye(2**n) - 1j*np.sin(op.angle)*p
        elif op.name == "cx":
            a,b = op.qubits
            za, xb = pauli_matrix(n,(a,),"Z"), pauli_matrix(n,(b,),"X")
            gate = (np.eye(2**n)+za)/2 + (np.eye(2**n)-za) @ xb / 2
        else:
            q = op.qubits[0]
            local = {"h": (PAULI["X"]+PAULI["Z"])/np.sqrt(2),
                     "s": np.diag([1,1j]), "sdg": np.diag([1,-1j]),
                     "z": PAULI["Z"], "x": PAULI["X"]}.get(op.name)
            if op.name == "rz":
                local = np.diag(np.exp(-.5j*op.angle*np.array([1,-1])))
            factors = [PAULI["I"] for _ in range(n)]
            factors[q] = local
            gate = tensor(factors)
        result = gate @ result
    return result


def annihilation(n, i):
    factors = [PAULI["Z"] if j < i else PAULI["I"] for j in range(n)]
    factors[i] = np.array([[0,1],[0,0]])
    return tensor(factors)


@pytest.mark.parametrize("i,j", [(0,1),(0,3),(3,1)])
def test_fermionic_hadamard_on_full_fock_space(i,j):
    n = 4
    f = operator(fermionic_hadamard(i,j), n)
    a,b = annihilation(n,i), annihilation(n,j)
    assert np.allclose(f @ a @ f.conj().T, (a+b)/np.sqrt(2))
    assert np.allclose(f @ b @ f.conj().T, (a-b)/np.sqrt(2))
    assert np.allclose(f @ f, np.eye(2**n))


def test_plaquette_gate_against_independent_hopping_exponential():
    n, time = 4, .173
    modes = (0,2,3,1)
    a = [annihilation(n,i) for i in range(n)]
    h = np.zeros((2**n,2**n), dtype=complex)
    for i,j in zip(modes, modes[1:]+modes[:1]):
        h -= a[i].conj().T @ a[j] + a[j].conj().T @ a[i]
    ops = plaquette_evolution(modes,time)
    assert np.allclose(operator(ops,n), expm(-1j*time*h), atol=1e-12)
    assert np.allclose(operator(native_path(ops),n), expm(-1j*time*h), atol=1e-12)
    routed,permutation,swaps=native_fermionic_route(ops,n,restore_order=True)
    assert permutation==list(range(n))
    assert swaps>0
    assert np.allclose(operator(routed,n),expm(-1j*time*h),atol=1e-12)


def test_greedy_fermionic_route_preserves_hopping_and_density():
    ops=[]
    for i,j,theta in [(0,3,.123),(1,2,-.17),(0,2,.07),(1,3,-.13)]:
        order=('XX','YY') if theta>0 else ('YY','XX')
        ops.extend(parity_rotation(i,j,p,theta) for p in order)
        # Density interaction uses physical site modes, which need not stay
        # adjacent after earlier fermionic mode permutations.
        ops.append(Operation('rotate',(0,1),.23,'ZZ','dynamic'))
    routed,permutation,_=native_fermionic_route(ops,4,restore_order=True)
    assert permutation==list(range(4))
    assert np.allclose(operator(routed,4),operator(ops,4),atol=1e-12)


@pytest.mark.parametrize("L,steps", [(2,1),(2,4),(4,1),(4,4),(6,3)])
@pytest.mark.parametrize("algorithm", ["plaquette", "termwise"])
def test_common_workload_matches_analytical_contract(L,steps,algorithm):
    ops, central = build_workload(L,1,steps,algorithm)
    counted = angle_summary(ops)
    expected = workload(L,1,steps,algorithm)
    assert counted["dynamic_rotations"] == expected.rotations
    assert counted["fixed_T"] == expected.fixed_T
    assert isclose(counted["theta_dynamic"], expected.theta, rel_tol=1e-12)
    assert central[1] == central[0]+1
    physical = schedule_native(native_path(ops), 2*L*L)
    assert physical["two_qubit_gates"] > 0
    assert physical["idle_qubit_layers"] >= 0
    assert physical["physical_copy_qubits"] >= 2*L*L
    routed,permutation,_=native_fermionic_route(ops,2*L*L)
    assert sorted(permutation)==list(range(2*L*L))
    assert schedule_native(routed,2*L*L)['idle_qubit_layers']>=0
