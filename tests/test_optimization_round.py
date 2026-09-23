"""Optimization-round identities: v1 replay, branch costing, ports, and residuals."""
import json
from math import log, sqrt

import numpy as np
import pytest

import optimization_round as o
from discussion_model import sampling_cost
from model import workload
from resource_extensions import synthesis_plan
from star_channel_study import rus_channel


@pytest.fixture(scope='module', autouse=True)
def channels():
    o.load_channels()


def v1_row(L, eps, arm):
    rows = json.loads((o.v1.OUT/'costs.json').read_text())
    return next(v for v in rows if (v['L'], v['tau'], v['absolute_error'], v['p'], v['arm'], v['scenario'])
                == (L, 1, eps, 1e-3, arm, 'component'))


@pytest.mark.parametrize('arm', ['star', 'hybrid', 'ft'])
def test_reference_menu_replays_integrated_v1(arm):
    sw = o.switches('v1_menu')
    row = (o.optimize_ft(4, 1, .01, 1e-3, 'component', sw) if arm == 'ft'
           else o.optimize_analog(4, 1, arm, .01, 1e-3, 'component', sw))
    assert row['log10_Q_min'] == pytest.approx(v1_row(4, .01, arm)['log10_Q_min'], abs=1e-12)


def test_optimal_richardson_allocation_matches_closed_form():
    # Equal variance, runtime proportional to steps: (1/3 + (4/3) sqrt 2)^2 base runtimes.
    t = 1e-3
    log_P = o.combined_log_P([0., 0.], [t, 2*t], [-1/3, 4/3])
    assert np.exp(log_P)/t == pytest.approx((1/3 + 4/3*sqrt(2))**2, rel=1e-12)
    assert np.exp(o.combined_log_P([0.7], [t], [1.]))/t == pytest.approx(np.exp(.7))


def test_single_branch_sampling_matches_discussion_ledger():
    c = o.base_config(1e-3, 'component')
    for q, s, chi, stat in [(6e4, .06, .7, .6), (2e6, 30., 9., .3)]:
        exact = sampling_cost(q, s, chi, {**c, 'stat_fraction': stat}, .01)['log10_Q_min']
        assert float(o.sampling_log_q(np.array([q]), np.array([chi+log(s)]), stat, c, .01)[0]) == pytest.approx(exact, abs=1e-12)


def test_richardson_residual_uses_only_signed_pairs():
    pairs = o.richardson_pairs('plaquette')
    four = [v['K4'] for v in pairs if v['patch'] == '4x4' and v['tau'] == 1.0]
    # Measured 4x4 plaquette errors fall exactly as r^-4 after extrapolation.
    assert min(four) == pytest.approx(.0514, rel=2e-3) and max(four) == pytest.approx(.0514, rel=2e-3)
    assert o.richardson_coefficient(4, 1, 'plaquette') < .1


@pytest.mark.parametrize('batch', [1, 4, 16])
def test_one_port_consumption_equals_serial_synthesis(batch):
    w = workload(8, 1, 10, 'plaquette')
    style = 'individual' if batch == 1 else 'hamming_weight'
    plan = synthesis_plan(w, .01, {'rotation_synthesis': style, 'hwp_batch_size': batch, 'synthesis_fraction': .05})
    assert o.serial_T_layers(w, plan, 1) == plan.get('synthesized_rotations', w.rotations)
    # Ports cannot parallelize beyond the independent rotations in each stage.
    assert o.serial_T_layers(w, plan, 10**6) == sum(r for _, r in o.plaquette_rotation_stages(w))


def test_compact_block_keeps_one_port_and_limited_lanes():
    sw = o.switches('all')
    assert o.port_menu(8, sw, 'compact') == [1] and max(o.port_menu(8, sw, 'fast')) > 1
    assert 'compact' not in o.analog_layouts('star', 8, sw)
    assert max(o.lane_menu(8, 'termwise', sw)) <= 8 and max(o.lane_menu(8, 'plaquette', sw)) <= 64


@pytest.mark.parametrize('theta', [.0123, .049, .2])
def test_pcec_never_has_less_dephasing_than_raw_rus(theta):
    raw = rus_channel(theta, 9, 1e-3, mode='raw', cap=48)
    pcec = rus_channel(theta, 9, 1e-3, mode='noisy_pcec', cap=48)
    assert pcec['q_dephasing'] >= raw['q_dephasing']
    assert abs(raw['coherent_component']) < 1e-6


def test_calibrated_round_reads_the_stim_fit_and_leaves_v1_on_the_reference_law():
    fit = json.loads((o.ROOT/'results/surface_law_v1/fit.json').read_text())['fits']
    try:
        o.use_round('v2')
        for p in (1e-3, 1e-4):
            c = o.base_config(p, 'component')
            assert c['surface_A'] == fit[str(p)]['A_patch_union_bound']
            assert c['p']/c['surface_threshold'] == pytest.approx(1/fit[str(p)]['Lambda'])
    finally:
        o.use_round('v1')
    c = o.base_config(1e-3, 'component')
    assert (c['surface_A'], c['surface_threshold']) == (.1, .01)


def test_source_menu_defaults_to_cultivation_and_matches_litinski_at_p():
    c = o.base_config(1e-3, 'component')
    menu = o.source_menu(c, {}, [0, 1])
    assert [s['name'] for s in menu] == ['cultivation']*2 and menu[1]['tiles'] == 13
    assert menu[1]['error'] == pytest.approx(35*menu[0]['error']**3)
    o.use_round('v3')
    try:
        litinski = o.source_menu(o.base_config(1e-4, 'component'), {'source_menu': 'litinski'}, [0, 1])[2:]
        assert litinski and all(s['tiles'] == 1 for s in litinski)
        # Four outputs per round share its 90.3 cycles; each output then needs two transfers.
        s = next(s for s in litinski if '20-to-4' in s['name'])
        assert o.source_service(s, c, 11, 11) == pytest.approx(90.3/4 + 22)
    finally:
        o.use_round('v1')


def test_fixed_T_analog_moves_every_pi_over_8_into_the_star_instrument():
    w = workload(8, 2, 4, 'plaquette', 4.)
    a = o.analog_fixed_T(w)
    assert a.fixed_T == 0 and a.rotations == w.rotations + w.fixed_T
    assert a.theta == pytest.approx(w.theta + w.fixed_T*np.pi/8)
