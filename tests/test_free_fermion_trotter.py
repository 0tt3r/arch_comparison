"""U = 0 Gaussian simulation of the Trotter circuits: agreement with a state vector and the builder."""
import pytest

import free_fermion_trotter as f


def test_gaussian_circuits_match_statevector_and_builder():
    rows = f.validate()
    assert max(v['difference'] for v in rows) < 1e-10


@pytest.mark.parametrize('algorithm', ['termwise', 'plaquette'])
def test_trotter_error_is_second_order(algorithm):
    exact = f.exact_doublons(4, 1.)[2, 2]
    a, b = (abs(f.trotter_doublons(4, 1., r, algorithm)[2, 2] - exact)*r*r for r in (16, 32))
    assert a == pytest.approx(b, rel=.05)
