"""Checks of mathematical identities and resource-accounting boundaries."""
import json
from math import exp, isclose, log, pi
from pathlib import Path

import numpy as np
import pytest

from model import (accuracy, estimate, factory, log_gamma2_dephasing,
                   log_gamma2_depolarizing, minimum_distance, star_exposure, workload)


@pytest.fixture
def config():
    return json.loads(Path(__file__).with_name("assumptions.json").read_text())


@pytest.mark.parametrize("qubits", [1, 2])
def test_pec_inverse_by_pauli_convolution(qubits):
    # Pauli labels modulo phases form (Z_2)^(2k); convolution must give delta_I.
    k, p = 4**qubits, .001
    noise = np.full(k, p / (k - 1))
    noise[0] = 1 - p
    convolution = np.array([[noise[i ^ j] for j in range(k)] for i in range(k)])
    delta = np.eye(k)[0]
    inverse = np.linalg.solve(convolution, delta)
    assert np.allclose(convolution @ inverse, delta, atol=1e-14)
    assert isclose(log_gamma2_depolarizing(p, qubits), 2 * log(abs(inverse).sum()), rel_tol=1e-10)


@pytest.mark.parametrize("L", [2, 4, 6, 12])
def test_termwise_counts_from_explicit_strang_sequence(L):
    N, E, r, tau = L**2, 2 * L * (L - 1), 3, .7
    coefficients = [1.] + [.5] * (4 * E) + [1.] * (N - 1)
    stage = [(i, a * tau / (2 * r)) for i, a in enumerate(coefficients)]
    merged = []
    for key, angle in (stage + stage[::-1]) * r:
        if merged and merged[-1][0] == key:
            merged[-1] = (key, merged[-1][1] + angle)
        else:
            merged.append((key, angle))
    w = workload(L, tau, r, "termwise")
    assert len(merged) == w.rotations
    assert isclose(sum(a for _, a in merged), w.theta, rel_tol=1e-12)
    assert isclose(w.theta, tau * (5 * N - 4 * L))


@pytest.mark.parametrize("L", [2, 4, 6, 10])
def test_open_boundary_plaquette_cover(L):
    covered = []
    for offset in (0, 1):
        for x in range(-offset, L - 1, 2):
            for y in range(-offset, L - 1, 2):
                verts = {(i, j) for i in (x, x + 1) for j in (y, y + 1)
                         if 0 <= i < L and 0 <= j < L}
                covered.extend(tuple(sorted((v, w))) for v in verts for w in verts
                               if v < w and abs(v[0]-w[0]) + abs(v[1]-w[1]) == 1)
        # The shifted tiling also has truncated tiles on the last row/column.
        if offset:
            for z in range(1, L - 1, 2):
                covered.extend([((L-1, z), (L-1, z+1)), ((z, L-1), (z+1, L-1))])
    assert len(covered) == len(set(covered)) == 2 * L * (L - 1)


def test_plaquette_spectrum_and_compiled_angle(config):
    # A four-cycle has two active eigenvalues; each spin needs two angles dt.
    adjacency = np.array([[0,1,0,1], [1,0,1,0], [0,1,0,1], [1,0,1,0]])
    assert np.allclose(np.linalg.eigvalsh(adjacency), [-2,0,0,2])
    a, b = workload(4, 1, 4), workload(4, 1, 8)
    assert a.theta == b.theta == 44
    assert b.fixed_T == 2 * a.fixed_T
    assert b.theta + b.fixed_T*pi/8 > a.theta + a.fixed_T*pi/8
    # Exact product channel has small O(sum q^2) step dependence; leading
    # angle-only exponent is invariant, fixed T contribution is not.
    assert star_exposure(b, config, "star_angle") > star_exposure(a, config, "star_angle")


def test_budget_and_factory_accounting(config):
    eps, n0 = accuracy(config)
    assert eps == .001
    assert n0 == 7377759
    assert factory(2, .001)[2] < factory(1, .001)[2]
    w = workload(4, 1, 4)
    for arm in ("nisq", "star_original", "star_angle", "star_floor", "ft"):
        row = estimate(w, config, arm)
        assert row["log10_Q_min"] >= log(row["qubits_per_copy"], 10) - 1e-12
        assert row["seconds_per_attempt"] > 0
    with pytest.raises(ValueError):
        accuracy({**config, "doublon_reference": 0})


def test_distance_is_minimal_and_pec_tightens_budget(config):
    duration = lambda d: .001 * d
    a = minimum_distance(100, duration, log(1e-4), config)
    b = minimum_distance(100, duration, log(1e-8), config)
    assert b[0] > a[0]
    d, _, _ = a
    previous_volume = 100 * np.ceil(duration(d-2) / config["cycle_seconds"])
    previous_failure = previous_volume * config["surface_A"] * (config["p"] / config["surface_threshold"])**((d-1)/2)
    assert previous_failure > 1e-4


def test_fixed_Q_runtime_has_its_own_configuration(config):
    _, n0 = accuracy(config)
    row = estimate(workload(8, 1, 4), config, "ft")
    replicas = config["device_qubits"] // row["fixed_Q_qubits_per_copy"]
    expected = log(n0 * row["fixed_Q_seconds_per_attempt"] / replicas, 10)
    assert isclose(row["log10_seconds_at_Q"], expected)


def test_dephasing_and_zero_noise():
    q = .01
    inverse = np.linalg.inv(np.array([[1-q, q], [q, 1-q]]))[:, 0]
    assert isclose(exp(log_gamma2_dephasing(q)), abs(inverse).sum()**2)
    assert log_gamma2_depolarizing(0, 2) == log_gamma2_dephasing(0) == 0
