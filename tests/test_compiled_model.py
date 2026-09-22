"""Checks of the review's logical-fault claim and compiled accounting."""

from project_paths import ROOT, project_file
import json
import re
from pathlib import Path

import numpy as np
import pytest

from compiled_model import compiled_physical, encoded


def test_branch_independent_fault_requires_a_commutation_condition():
    # Bloch channels: Z dephasing, a stochastic X fault, and a stochastic H
    # fault. Both faults are independent of the sampled PEC correction.
    s,q=.2,.01
    noise=np.diag([s,s,1.])
    pauli=(1-q)*np.eye(3)+q*np.diag([1.,-1.,-1.])
    hadamard=(1-q)*np.eye(3)+q*np.array([[0.,0.,1.],[0.,-1.,0.],[1.,0.,0.]])
    assert np.allclose(np.linalg.solve(noise,pauli@noise),pauli)
    assert not np.allclose(np.linalg.solve(noise,hadamard@noise),hadamard)
    # On |0>, D=(I+X)/2: PEC^{-1} F E has bias q/(2s)>q.
    mitigated=np.linalg.solve(noise,hadamard@noise)@np.array([0.,0.,1.])
    assert np.isclose(mitigated[0]/2,q/(2*s))
    assert mitigated[0]/2 > q


def test_compiled_costs_and_idle_sensitivity():
    root=ROOT
    c=json.loads((project_file('assumptions.json')).read_text())
    record=json.loads((root/'results/cudaq_v1/plaquette_L4_tau1_r4/result.json').read_text())
    time=compiled_physical(record,c)
    layer=compiled_physical(record,c,'per_layer')
    assert time['chi_idle'] < layer['chi_idle']/9.9
    assert np.isclose(time['chi_two_qubit'],layer['chi_two_qubit'])
    assert encoded(record,c)['T_count'] == record['synthesis']['T_total']
    assert encoded(record,c,'hybrid')['T_count'] == record['angles']['fixed_T']


@pytest.mark.parametrize('algorithm',['plaquette','termwise'])
@pytest.mark.parametrize('L',[2,4,6,8])
def test_actual_synthesized_words_meet_operator_precision(algorithm,L):
    """Independently reconstruct each emitted one-qubit Clifford+T word."""
    directory=ROOT/'results/cudaq_v1'/f'{algorithm}_L{L}_tau1_r4'
    record=json.loads((directory/'result.json').read_text())
    tolerance=record['synthesis']['operator_precision_per_rotation']
    gates={'h':np.array([[1,1],[1,-1]])/np.sqrt(2),'s':np.diag([1,1j]),
           't':np.diag([1,np.exp(1j*np.pi/4)]),'x':np.array([[0,1],[1,0]])}
    for i,angle in enumerate(record['synthesis']['angle_classes']):
        unitary=np.eye(2,dtype=complex)
        mlir=(directory/f'angle_{i}_ft.mlir').read_text()
        for name in re.findall(r'qlx.apply #qlx.action<(\w+)>',mlir):
            unitary=gates[name]@unitary
        target=np.diag(np.exp(-1j*angle['theta']*np.array([1,-1])))
        phase=np.angle(np.trace(target.conj().T@unitary))
        error=np.linalg.norm(unitary*np.exp(-1j*phase)-target,2)
        assert error<=tolerance+1e-13
