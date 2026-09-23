"""Independent controls for the same-task algorithm and component comparisons."""
import gzip
import json
from math import exp,log,sqrt

import numpy as np
import pytest
from numpy.polynomial.chebyshev import chebval
from scipy.linalg import expm
from scipy.special import jv

from geometric_delivery import validate_events
from locality_certificate import log_bound,sufficient_radius
from project_paths import ROOT
from qubitization_study import AMPLIFICATION,OpenHubbard,amplification,bessel_degree,query_budget


def direct_hubbard(L):
    sites=[(x,y) for y in range(L) for x in (range(L) if y%2==0 else range(L-1,-1,-1))]
    n=2*L*L;H=np.zeros((2**n,2**n),complex)
    for state in range(2**n):
        H[state,state]=sum(4*(((state>>(2*i))&1)-.5)*(((state>>(2*i+1))&1)-.5) for i in range(L*L))
        for i,a in enumerate(sites):
            for j,b in enumerate(sites):
                if j<=i or sum(abs(x-y) for x,y in zip(a,b))!=1:continue
                for spin in (0,1):
                    u,v=2*i+spin,2*j+spin
                    if ((state>>u)^(state>>v))&1:
                        between=((1<<v)-1)^((1<<(u+1))-1)
                        H[state^(1<<u)^(1<<v),state]=-(-1)**((state&between).bit_count())
    return H


def test_open_boundary_lcu_matches_independent_fermionic_hamiltonian():
    paulis={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
    model=OpenHubbard(2);H=np.zeros((256,256),complex)
    for word,c in model.terms.items():
        P=np.array([[1.]])
        for letter in reversed(word):P=np.kron(P,paulis[letter])
        H+=c*P
    assert np.max(abs(H-direct_hubbard(2)))<1e-13
    assert model.alpha==12
    assert sum(round(2*abs(v)) for v in model.terms.values())==24
    for L in (4,6,8):assert OpenHubbard(L).alpha==4*L*(L-1)+L*L


@pytest.mark.parametrize('x',[.2,12,64,512])
def test_bessel_tail_dominates_complex_polynomial_error(x):
    d,bound=bessel_degree(x,1e-7)
    coeff=np.array([jv(0,x)]+[2*(-1j)**n*jv(n,x) for n in range(1,d+1)])
    z=np.linspace(-1,1,3001)
    error=max(abs(chebval(z,coeff)-np.exp(-1j*x*z)))
    assert error<=bound+2e-12 and bound<=1e-7
    # The independent full-Hilbert-space spectral calculation is another
    # check of the normalization, rather than a query-count identity.
    if x==12:
        H=direct_hubbard(2);w=np.linalg.eigvalsh(H)/12
        assert max(abs(chebval(w,coeff)-np.exp(-1j*x*w)))<=bound+1e-12


@pytest.mark.parametrize('j',[0,1,2,3])
def test_oblivious_amplification_and_perturbation_bound(j):
    z=.25*np.exp(.713j);b=sqrt(1-abs(z)**2)
    U=np.array([[z,b],[b,-z.conjugate()]]);R=np.diag([1,-1])
    actual=np.linalg.matrix_power(-U@R@U.conj().T@R,j)@U
    expected=amplification(j,0.)
    assert abs(actual[0,0])**2==pytest.approx(expected['acceptance_lower'],abs=1e-14)
    rng=np.random.default_rng(45);ideal=.25*expm(1j*np.array([[.3,.2],[.2,-.4]]))
    E=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2));E*=1e-5/np.linalg.norm(E,2)
    def transform(A):
        return sum(c*np.linalg.matrix_power(A@A.conj().T,(n-1)//2)@A for n,c in enumerate(AMPLIFICATION[j]) if c)
    bound=amplification(j,1e-5)
    assert np.linalg.norm(transform(ideal+E)-transform(ideal),2)<=bound['block_error']


def test_query_and_campaign_error_accounting():
    data=json.loads((ROOT/'results/frontier_v1/qubitization/results.json').read_text())
    for v in data['rows']:
        if not v['variant'].startswith('trotter'):
            budget=v['budget_fractions'];assert sum(budget.values())==pytest.approx(1,abs=2e-13)
            assert min(budget.values())>=0
            assert 0<v['attempt_acceptance_lower']<=v['acceptance_lower']<=1
            assert v['T_count']>=v['direct_T']
            assert v['polynomial_tail']<=2*v['base_block_error']
            assert v['observable_error']<=v['epsilon']*budget['simulation']*(1+1e-12)
            assert v['block_queries']==(v['cos_degree']+v['sin_degree'])*(2*v['amplification_steps']+1)


def test_split_uniform_preparation_is_correct_and_has_no_magnitude_in_select():
    import cirq
    from qualtran.bloqs.state_preparation import PrepareUniformSuperposition
    # Simulate the library PREPARE independently of the inventory adapter.
    for n in (24,128):
        k=(n-1).bit_length();q=cirq.LineQubit.range(k)
        context=cirq.DecompositionContext(cirq.GreedyQubitManager('work',maximize_reuse=True))
        circuit=cirq.Circuit(cirq.decompose(PrepareUniformSuperposition(n).on_registers(target=q),context=context))
        order=q+sorted(circuit.all_qubits()-set(q))
        result=cirq.Simulator(dtype=np.complex128,seed=3).simulate(circuit,qubit_order=order).final_state_vector
        probability=abs(result.reshape(2**k,-1))**2
        assert np.max(abs(probability.sum(axis=1)-np.r_[np.full(n,1/n),np.zeros(2**k-n)]))<1e-12
        assert probability[:,1:].sum()<1e-12
    inv=json.loads((ROOT/'results/frontier_v1/qubitization/results.json').read_text())['inventories']
    row=next(v for v in inv if v['L']==4 and v['variant']=='split_uniform')
    assert row['lcu_terms']==128 and row['encoding']['Rotations']==0


def test_magic_transport_restores_entangled_data_in_both_outcomes():
    # Three arbitrary data qubits, one |T> output port; path 0->1->2, target 3.
    rng=np.random.default_rng(381);data=rng.normal(size=8)+1j*rng.normal(size=8);data/=np.linalg.norm(data)
    psi=np.zeros(16,complex);psi[::2]=data/sqrt(2);psi[1::2]=np.exp(1j*np.pi/4)*data/sqrt(2)
    idx=np.arange(16)
    def swap(v,a,b):return v[idx^((((idx>>a)^(idx>>b))&1)*((1<<a)|(1<<b)))]
    moved=swap(swap(psi,0,1),1,2)
    moved=moved[idx^(((idx>>3)&1)<<2)]
    for outcome in (0,1):
        branch=moved.copy();branch[((idx>>2)&1)!=outcome]=0
        if outcome:
            branch[((idx>>3)&1)==1]*=1j
            branch=branch[idx^(1<<2)]
        branch=swap(swap(branch,1,2),0,1)
        target=np.zeros(16,complex);target[::2]=data*np.exp(1j*np.pi/4*((np.arange(8)>>2)&1))
        assert abs(np.vdot(target,branch))**2==pytest.approx(.5,abs=1e-14)


def test_delivery_reservations_and_counts():
    folder=ROOT/'results/frontier_v1/delivery'
    data=json.loads((folder/'results.json').read_text())
    for row in data['rows']:
        with gzip.open(folder/f"batch{row['batch']}_ports{row['ports']}.json.gz",'rt') as f:events=json.load(f)
        validate_events(events)
        deliveries=[v for v in events if v['kind']=='T_delivery']
        assert len(deliveries)==row['T_count']
        assert sum(v['swaps'] for v in deliveries)==row['delivery_SWAPS']
        assert row['seconds']>=row['abstract_delivery_seconds']
        assert row['logical_clocks']>=row['source_only_clocks']


def test_locality_bound_against_exact_operator_dynamics():
    H=direct_hubbard(2);idx=np.arange(256);D=np.diag(((idx&1)!=0)&((idx&2)!=0)).astype(complex)
    # Radius zero retains only the onsite term, which commutes with D.
    for t in (.001,.01,.03):
        U=expm(-1j*t*H);difference=np.linalg.norm(U.conj().T@D@U-D,2)
        assert difference<=exp(log_bound(0,t))+1e-13
    for t in (.1,1,8):
        r=sufficient_radius(t,.0025)
        assert exp(log_bound(r,t))<=.0025<exp(log_bound(r-1,t))
    assert sufficient_radius(1,.0025)==78


def test_growth_zero_failures_are_finite_simultaneous_bounds():
    data=json.loads((ROOT/'results/frontier_v1/growth_summary.json').read_text())
    assert data['shots']==24000000
    for case in data['cases']:
        if case['failures']==0:
            expected=1-(.05/24)**(1/case['shots'])
            assert case['simultaneous_upper95']==pytest.approx(expected,rel=1e-10)
            assert case['simultaneous_upper95']>case['interval95'][1]>0


def test_exposure_preserves_heldout_failures_and_training_cost():
    data=json.loads((ROOT/'results/frontier_v1/exposure.json').read_text())
    assert len(data['nisq'])==384
    assert any(not v['passes_bias_allowance'] for v in data['nisq'])
    for v in data['nisq']:
        assert v['exposure']==pytest.approx(sum(v['components'].values()))
        assert v['idle_locations']>=0
        if v['method'].startswith('tflo'):
            assert v['training_attempts']==6000000
            assert v['total_time_multiplier_with_training'] is not None
