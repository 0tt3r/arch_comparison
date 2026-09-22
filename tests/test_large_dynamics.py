"""Independent small-sector controls for the large matrix-free implementation."""
import numpy as np
import pytest
import large_dynamics
from large_dynamics import coefficients, compile_solver, solve
from time_calibration import SplitSector


@pytest.fixture(scope='module',autouse=True)
def solver(tmp_path_factory):
    # Tests must not replace an executable that a PBS calculation is using.
    previous=large_dynamics.BINARY
    large_dynamics.BINARY=tmp_path_factory.mktemp('solver')/'hubbard_sector'
    compile_solver()
    yield
    large_dynamics.BINARY=previous


@pytest.mark.parametrize('shape',[(2,2),(3,2),(3,3)])
@pytest.mark.parametrize('algorithm',['exact','plaquette','termwise'])
def test_full_state_against_independent_sector(shape,algorithm,tmp_path):
    nx,ny=shape;model=SplitSector(nx,ny);tau=.71;r=3
    path=tmp_path/'state.bin'
    record=solve(nx,ny,tau,algorithm,0 if algorithm=='exact' else r,tmp_path,2,state_path=path)
    actual=np.fromfile(path,dtype=np.complex128).reshape(model.shape)
    # C++ masks use ascending integer order; scipy uses combination order.
    up=np.argsort(np.argsort(model.up));down=np.argsort(np.argsort(model.down))
    actual=actual[np.ix_(up,down)].ravel()
    expected=model.exact(tau) if algorithm=='exact' else model.trotter(tau,r,algorithm)
    assert np.linalg.norm(actual-expected)<2e-11
    assert abs(record['result']['norm_squared']-1)<2e-12


def test_tail_bound_controls_scalar_exponential():
    for a in (.1,10,64,128,512):
        c,tail=coefficients(a,1)
        x=np.linspace(-1,1,101)
        assert tail<1e-12
        assert np.max(abs(np.polynomial.chebyshev.chebval(x,c)-np.exp(-1j*a*x)))<2e-12
