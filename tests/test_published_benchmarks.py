"""Independent controls for the published-input replay and experimental reduction."""
import io
import json
import pickle

import numpy as np
import pytest

from phasecraft_data import DataOnlyReader, shot_observables
from published_benchmarks import OUT, qcels_branches, qcels_time, read, topqad_replay
from project_paths import ROOT


def test_qcels_time_grid_sum_and_integer_rounding():
    r,J,N=5*2**10,10,5
    branches=qcels_branches(r,J,N)
    assert len(branches)==J*N
    assert max(t for t,_ in branches)==1
    assert sum(steps for _,steps in branches)==pytest.approx(r*(N+1)*(1-2**(-J)))
    assert all(steps>=37*t and steps-37*t<1 for t,steps in qcels_branches(37,J,N))


def test_two_quadratures_and_control_operations_are_counted():
    row=dict(L=4,U=4,epsilon=.08,r=5*2**10,J=10,di=9,df=19,clocks=314)
    cycle=350e-9;branch_steps=row['r']*6*(1-2**-10)
    one=100*branch_steps*314*34*cycle
    assert qcels_time(row,cycle,quadratures=1,control=False)==pytest.approx(one)
    assert qcels_time(row,cycle,control=False)==pytest.approx(2*one)
    overhead=200*(18*branch_steps+16*50)*34*cycle
    assert qcels_time(row,cycle)==pytest.approx(2*one+overhead)


def test_growth_increases_cost_and_clock_change_is_linear():
    row=read(ROOT/'configs/published_benchmarks.json')['star'][0]
    t=qcels_time(row,350e-9,beta=.001)
    assert qcels_time(row,350e-9,beta=.001,growth=1e-5)>t
    assert qcels_time(row,240e-9,beta=.001)/t==pytest.approx(240/350)


def test_published_inconsistencies_are_not_silently_repaired():
    star,ft=topqad_replay(read(ROOT/'configs/published_benchmarks.json'))
    odd=next(v for v in star if (v['U'],v['L'],v['hardware'])==(8,8,'desired'))
    assert odd['clocks']==750
    assert odd['single_seconds_replayed']==pytest.approx(.4704)
    assert odd['single_seconds_published']==.297
    below=next(v for v in star if (v['U'],v['L'],v['hardware'])==(8,6,'desired'))
    assert below['implied_growth_per_state'] is None
    assert below['published_variance_multiplier']<below['angle_only_multiplier']
    # Area/time tradeoffs are preserved, not reduced to one capacity ranking.
    same=[v for v in ft if (v['U'],v['L'],v['hardware'])==(4,8,'target')]
    space=next(v for v in same if v['mode']=='space');fast=next(v for v in same if v['mode']=='time')
    assert space['qubits_published']<fast['qubits_published']
    assert space['seconds_published']>fast['seconds_published']


def test_fitted_growth_transfers_without_refitting_other_rows():
    star,_=topqad_replay(read(ROOT/'configs/published_benchmarks.json'))
    for hardware in ('target','desired'):
        group=[v for v in star if v['hardware']==hardware]
        assert len(set(v['anchor_growth_per_state'] for v in group))==1
        assert sum(v['used_to_fit_growth'] for v in group)==1


def test_doublon_reconstruction_respects_qubit_labels_and_correlations():
    labels=[f'({x},{y},{spin})' for spin in ('up','down') for x in range(7) for y in range(4)]
    bits=np.zeros((2,56),dtype=np.uint8)
    bits[0,:14]=1;bits[0,28:42]=1
    bits[1,14:28]=1;bits[1,28:42]=1
    D,accepted=shot_observables(bits,labels)
    assert D.mean(axis=1).tolist()==[.5,0]
    assert accepted.tolist()==[True,True]
    perm=np.random.default_rng(9).permutation(56)
    shuffled,_=shot_observables(bits[:,perm],[labels[i] for i in perm])
    np.testing.assert_array_equal(shuffled,D)


def test_external_dataframe_reader_rejects_executable_globals():
    with pytest.raises(ValueError,match='Unsupported external serialized type'):
        DataOnlyReader(io.BytesIO(b'cos\nsystem\n.')).load()
    value=np.array([['Raw',(1,2)]],dtype=object)
    recovered=DataOnlyReader(io.BytesIO(pickle.dumps(value,protocol=5))).load()
    assert recovered[0,0]=='Raw' and recovered[0,1]==(1,2)


def test_experimental_census_does_not_promote_interacting_references_to_truth():
    data=read(OUT/'phasecraft_reduced.json')
    assert len(data['shots'])==40
    assert sum(v['shots'] for v in data['shots'])==6400
    assert all(v['qubits']==56 for v in data['shots'])
    assert not any(v['U']==4 and v['method'] in ('Exact','FLO') for v in data['observables'])
    # The archive group named 'exact' also contains approximate Majorana
    # propagation. Selecting by group alone would silently mix references.
    assert any(v['U']==4 and v['group']=='exact' and v['method']=='Majorana Propagation' for v in data['observables'])
