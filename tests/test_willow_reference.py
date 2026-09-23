"""Independent Fock-space controls for the experimental circuit and center-hole sector."""
import json
import numpy as np
import pytest
from scipy.linalg import expm

from small_dynamics import rectangle
from willow_reference import hopping_order, free_circuit
from willow_benchmark import free_dynamics
from willow_sector import compile_solver, masks, solve
from willow_reference_analysis import decompose, observables


@pytest.fixture(scope='module')
def binary(tmp_path_factory):
    path=tmp_path_factory.mktemp('willow_solver')/'sector'
    compile_solver(path)
    return path


def independent_fock():
    n=4;sites=rectangle(2,2);ann=[]
    for i in range(2*n):
        a=np.zeros((1<<(2*n),)*2)
        for state in range(1<<(2*n)):
            if (state>>i)&1:a[state^(1<<i),state]=(-1)**((state&((1<<i)-1)).bit_count())
        ann.append(a)
    number=[a.T@a for a in ann]
    V=sum(4*(number[i]-.5*np.eye(256))@(number[i+n]-.5*np.eye(256)) for i in range(n))
    hops=[]
    for i,j in hopping_order(2,2):
        for s in (0,n):hops.append(-(ann[i+s].T@ann[j+s]+ann[j+s].T@ann[i+s]))
    iu,idown=masks(2)
    initial=np.zeros(256,complex);initial[iu+(idown<<n)]=1
    up=[m for m in range(1<<n) if m.bit_count()==iu.bit_count()]
    down=[m for m in range(1<<n) if m.bit_count()==idown.bit_count()]
    indices=[u+(d<<n) for u in up for d in down]
    return initial,V,hops,indices,number


@pytest.mark.parametrize('mode',['exact','circuit','trajectory'])
def test_product_sector_state_and_observables_against_full_fock(mode,binary,tmp_path):
    tau=.23;r=3 if mode!='exact' else 0
    state=tmp_path/'state.bin'
    records=solve(2,tau,mode,r,tmp_path,threads=2,binary=binary,state_path=state)
    actual=np.fromfile(state,np.complex128)
    initial,V,hops,indices,number=independent_fock()
    if mode=='circuit':
        expected=initial
        for _ in range(r):
            for h in hops:expected=expm(-1j*tau/(2*r)*h)@expected
            expected=expm(-1j*tau/r*V)@expected
            for h in hops[::-1]:expected=expm(-1j*tau/(2*r)*h)@expected
    else:expected=expm(-1j*tau*(r if mode=='trajectory' else 1)*(V+sum(hops)))@initial
    np.testing.assert_allclose(actual,expected[indices],atol=2e-12,rtol=0)
    out=records[-1]['result']
    for i in range(4):
        D=np.vdot(expected,number[i]@number[i+4]@expected).real
        assert out['doublons'][i]==pytest.approx(D,abs=2e-12)
        for s,field in [(0,'z_up'),(4,'z_down')]:
            assert out[field][i]==pytest.approx(1-2*np.vdot(expected,number[i+s]@expected).real,abs=2e-12)


@pytest.mark.parametrize('shape',[(2,2),(3,3),(4,4),(6,6)])
def test_swap_network_covers_each_physical_bond_once(shape):
    nx,ny=shape;sites=rectangle(nx,ny)
    order=hopping_order(nx,ny)
    assert len(order)==2*nx*ny-nx-ny==len(set(order))
    assert all(sum(abs(a-b) for a,b in zip(sites[i],sites[j]))==1 for i,j in order)


def test_free_circuit_converges_to_independent_continuous_orbitals():
    L=3
    occupations=[dict(site=[x,y],spin=s,n=int((x,y)!=(0,0) and ((x+y)%2==0)==(s=='up')))
                 for x in range(L) for y in range(L) for s in ('up','down')]
    exact=np.array(free_dynamics(L,.73,occupations,[.5])[.5]['local_z'])
    errors=[np.max(abs(np.array(free_circuit(L,.73,occupations,.5,r)['local_z'])-exact)) for r in (16,32)]
    assert 3.9<errors[0]/errors[1]<4.1


def test_reusing_a_result_rejects_changed_input(binary,tmp_path):
    solve(2,.1,'exact',0,tmp_path,threads=1,binary=binary)
    with pytest.raises(AssertionError,match='Stale contract'):
        solve(2,.1,'exact',0,tmp_path,threads=2,binary=binary)


def test_signed_error_decomposition_retains_cancellation():
    target=np.array([.1,.2]);circuit=np.array([.2,.1]);measured=np.array([.11,.19])
    v=decompose(measured,circuit,target)
    assert v['experiment_vs_continuous_rms']==pytest.approx(.01)
    assert v['experiment_vs_circuit_rms']==pytest.approx(.09)
    assert v['circuit_vs_continuous_rms']==pytest.approx(.1)
    assert v['twice_cross_mean']==pytest.approx(-.018)


def test_diagonal_symmetry_acts_on_ideal_observables_not_on_states():
    result=dict(doublons=[0.,.1,.2,.3],z_up=[1.,.8,.4,.2],z_down=[-.2,0.,.2,.6])
    raw=observables(result,L=2,symmetrize=False);sym=observables(result,L=2)
    d=raw['local_d'].reshape(2,2)
    np.testing.assert_allclose(sym['local_d'],((d+d.T)/2).ravel())
    np.testing.assert_allclose(sym['mean_d'],raw['mean_d'])
