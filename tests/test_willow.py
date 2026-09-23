"""Controls for archive integrity, observables, and randomized-circuit sampling."""
import io
import pickle
import zipfile

import h5py
import numpy as np
import pytest

from small_dynamics import HubbardSector,rectangle
from willow_archive import archive_index,selected_members,unpack_member
from willow_benchmark import combine_paulis,design_forecast,free_dynamics
from willow_data import grouped_statistics,initial_state_issues,nonempty_columns,shot_values


def test_selected_zip_member_integrity_and_corruption():
    buffer=io.BytesIO()
    with zipfile.ZipFile(buffer,'w',compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr('example.txt',b'Pinned scientific data\n'*30)
    buffer.seek(0);entry=archive_index(buffer)[0]
    assert unpack_member(buffer,entry)==b'Pinned scientific data\n'*30
    with pytest.raises(ValueError,match='CRC'):
        unpack_member(buffer,dict(entry,crc=entry['crc']^1))


def test_archive_selection_requires_every_declared_case():
    cfg=dict(sizes=[4],fluxes=['0'],interactions=[0],steps=[3],initial_state='center_hole')
    prefix='programmable_simulation_fermi_hubbard/'
    entries=[dict(name=prefix+'initial_states.json')]+[dict(name=prefix+
        'center_hole/4x4/mag_field_0/U_0/trotter_layers_3/'+f) for f in ('shots.h5','exp_vals.h5')]
    assert selected_members(cfg,entries)==entries
    with pytest.raises(ValueError,match='complete declared workload'):
        selected_members(cfg,entries[:-1])


def test_empty_pandas_frame_is_missing_data_not_exact_truth(tmp_path):
    with h5py.File(tmp_path/'empty.h5','w') as h:
        group=h.create_group('exact');axis=group.create_dataset('axis1',data=[0.])
        axis.attrs['shape']=np.bytes_(pickle.dumps((0,),protocol=0))
        assert nonempty_columns(group) is None


def test_bad_initial_state_metadata_is_flagged():
    state=dict(ups=[[0,0]],downs=[[0,5]],holons=[[2,4]],doublons=[],singlet_pairs=[[[1,0],[1,1]]])
    issues=initial_state_issues({'center_hole':{'4x4':state}},4)
    assert issues['out_of_bounds']==[(0,5),(2,4)]
    assert not issues['expected_center_hole_product_state']


def test_pauli_combination_matches_direct_occupation_statistics():
    rng=np.random.default_rng(832)
    bits=rng.integers(0,2,size=(157,2,2,2))
    z=1-2*bits
    result=combine_paulis(z.mean(0),(z[:,:,:,0]*z[:,:,:,1]).mean(0))
    D=bits[:,:,:,0]*bits[:,:,:,1]
    np.testing.assert_allclose(result['local_d'],D.mean(0).ravel(),atol=1e-15)
    expected=((bits[:,:,:,0]-bits[:,:,:,1])*np.array([[1,-1],[-1,1]])).sum((1,2))/8
    assert result['staggered_m'][0]==pytest.approx(expected.mean())


def test_covariance_agnostic_error_bound_covers_correlated_inputs():
    z=np.zeros((2,2,2));zz=np.zeros((2,2));zs=np.full_like(z,.02);zzs=np.full_like(zz,.03)
    result=combine_paulis(z,zz,zs,zzs)
    # Perfect anticorrelation of the ZZ fluctuation with both Z fluctuations
    # attains the upper bound; independent-error propagation would miss it.
    assert result['se']['local_d']==pytest.approx([(.02+.02+.03)/4]*4)


def test_randomized_circuit_variance_is_not_iid_shot_variance():
    groups=[np.full((100,1),v) for v in (.1,.2,.3,.4)]
    result=grouped_statistics(groups)
    assert result['cluster_se'][0]==pytest.approx(np.std([.1,.2,.3,.4],ddof=1)/2)
    between=result['between_mean_variance'][0];within=result['within_mean_variance'][0]
    assert design_forecast(between,within,4,shot_multiplier=100)==pytest.approx(result['cluster_se'][0])
    assert design_forecast(between,within,4,twirl_multiplier=100)==pytest.approx(result['cluster_se'][0]/10)


def test_postselection_does_not_reweight_randomized_circuits_by_acceptance():
    groups=[np.full((n,1),v) for n,v in [(10,.1),(100,.9)]]
    assert grouped_statistics(groups)['mean'][0]==pytest.approx(.5)
    assert np.concatenate(groups).mean()!=pytest.approx(.5)


def test_shot_mapping_ignores_ancillas_and_obeys_permutation():
    labels=[f'({x},{y},{s})' for x in range(2) for y in range(2) for s in ('up','down')]+['ancilla','ancilla']
    bits=np.random.default_rng(32).integers(0,2,size=(40,10),dtype=np.uint8)
    values,up,down=shot_values(bits,labels,2)
    bits[:,-2:]=1-bits[:,-2:]
    np.testing.assert_array_equal(shot_values(bits,labels,2)[0],values)
    perm=np.random.default_rng(33).permutation(10)
    shifted=shot_values(bits[:,perm],[labels[i] for i in perm],2)
    for actual,expected in zip(shifted,(values,up,down)):np.testing.assert_array_equal(actual,expected)


def test_free_dynamics_matches_independent_many_body_sector_solver():
    model=HubbardSector(rectangle(2,2),u=0)
    occupations=[dict(site=[x,y],spin=s,n=int(((x+y)%2==0)==(s=='up')))
                 for x in range(2) for y in range(2) for s in ('up','down')]
    result=free_dynamics(2,0.,occupations,[0.,.3,1.2])
    for t,observables in result.items():
        psi=model.exact(t)
        expected=[model.doublon(psi,(x,y)) for x in range(2) for y in range(2)]
        np.testing.assert_allclose(observables['local_d'],expected,atol=2e-14)
        assert sum((1-np.array(observables['local_z']))/2)==pytest.approx(4.)


def test_free_flux_reversal_leaves_product_density_invariant():
    occupations=[dict(site=[x,y],spin=s,n=int(((x+y)%2==0)==(s=='up')))
                 for x in range(3) for y in range(3) for s in ('up','down')]
    a=free_dynamics(3,.73,occupations,[.43])[.43]
    b=free_dynamics(3,-.73,occupations,[.43])[.43]
    np.testing.assert_allclose(a['local_z'],b['local_z'],atol=1e-14)
