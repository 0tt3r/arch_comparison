from math import ceil, sqrt
import pytest
from discussion_model import assumptions, data_block_tiles, estimate
from budget_optimization import optimize
from delivery_sensitivity import choices, settings


def test_paper_area_formulas_and_charged_lanes():
    c=assumptions();n=38
    for block,area in [('compact',1.5*n+3),('intermediate',2*n+4),('fast',2*n+sqrt(8*n)+1)]:
        assert data_block_tiles(4,{**c,'data_block':block})==ceil(area)+6*4
        assert data_block_tiles(4,{**c,'data_block':block,'parallel_lanes':8})==ceil(area)+6*8
    with pytest.raises(ValueError):data_block_tiles(4,{**c,'data_block':'unknown'})


def test_whole_block_rate_does_not_count_lanes_twice():
    c={**assumptions(),'data_block':'fast','ft_clocks_per_T':1,'T_consumption_ports':1}
    for lanes in (1,4,16):
        row=optimize(4,1,'ft',{**c,'parallel_lanes':lanes},.001)
        from model import workload
        rotations=workload(4,1,row['r'],'plaquette',c['u_over_j']).rotations
        expected=(row['T_count']+2*rotations/lanes)*row['distance']*c['cycle_seconds']
        assert row['consumption_seconds']==pytest.approx(expected)


def test_supply_is_charged_to_each_selected_factory():
    c=assumptions();row=optimize(4,1,'ft',c,.001)
    d=row['distance'];f=row['factories'];cleanup=row['cleanup_levels']
    assert row['qubits_per_copy']==2*(row['data_tiles']+f*(1+12*cleanup))*d*d+f*c['source_qubits']
    assert row['supply_seconds']==pytest.approx(ceil(row['T_count']/f)*row['source_service_cycles']*c['cycle_seconds'])
    assert row['seconds_per_attempt']>=max(row['supply_seconds'],row['consumption_seconds'])


def test_port_assumption_is_explicit_in_menu():
    for arm in ('ft','hybrid'):
        for name,c in choices(arm,settings()):
            if name.startswith(('compact_','intermediate_','fast_b')):
                assert c['T_consumption_ports']==1
            if name=='fast_four_ports':assert c['T_consumption_ports']==4


def test_measured_step_menu_cannot_silently_refine():
    from budget_optimization import settings as budget_settings
    c=assumptions();opts={**budget_settings(),'explicit_step_candidates':[32]}
    row=optimize(4,1,'star',c,.001,opts,C_override=.2)
    assert row['r']==32
    assert row['assumed_trotter_error']==pytest.approx(.2/32**2)
    with pytest.raises(ValueError):optimize(4,1,'star',c,.001,{**opts,'explicit_step_candidates':[4]})
