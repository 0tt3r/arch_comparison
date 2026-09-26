"""Cross-check of the Pinnacle ledger (pinnacle.py) against CUDA-Q Logical's Pinnacle architecture.

Compares, value by value, the processing-block rows and magic-engine operating points that CUDA-Q Logical
vendors from arXiv:2602.11457v2, then compiles a one-step 2x2 Hubbard plaquette circuit on a Pinnacle
device (one gb30 processing block, its RUS scratch block, and the p = 1e-3 engine) through placement and
the QEC pipeline, and records what the compiler does with each T gate. Two routes:

  A  Clifford+T synthesis, then Pauli-based computation: every T becomes one engine request and one
     resource-assisted Pauli-product measurement ("inject") on the processing block.
  B  arbitrary rotations lowered by Pinnacle's own repeat-until-success synthesis (the paper's route).

Physical (P3) lowering of the GB codes and the analytical/schedule tiers are attempted and their status
recorded; in cudaq-logical 0.1.1 they are not available for these codes (see result.json).

    .venv-cudaq/bin/python python/cudaq_pinnacle_check.py
"""
import json
import sys
import traceback
from collections import Counter
from math import ceil, isclose

import cudaq.logical as cql
from cudaq.logical.architectures import pinnacle as P

import pinnacle as ledger
from cudaq_compile import program_from_operations
from hubbard_circuits import build_workload
from project_paths import ROOT

OUT = ROOT/'results/cudaq_pinnacle_check_v1'
CYCLE_NS = 240.0
P_PHYS = 1e-3
PRECISION = 1e-6


def tables():
    rows = []
    for code in ledger.CODES:
        block = P.processing_block(getattr(P, code.name))
        inst = P._PINNACLE_GB_INSTANCES[code.name]
        rows.append(dict(item=code.name, ledger=[code.n, code.k, code.d, code.block_qubits, code.cycle_rounds],
                         cudaq=[inst.n, block.logical_capacity, inst.distance, block.physical_qubits,
                                inst.logical_cycle_rounds]))
    by_name = {e.name: e for e in P._MAGIC_ENGINES}
    names = {'p1e-4_1e-9': 'pinnacle_me_p1e_4_pout1e_9', 'p1e-4_1e-11': 'pinnacle_me_p1e_4_pout1e_11',
             'p1e-3_1e-9': 'pinnacle_me_p1e_3_pout1e_9', 'p1e-3_1e-11': 'pinnacle_me_p1e_3_pout1e_11'}
    for e in ledger.ENGINES:
        c = by_name[names[e.name]]
        rows.append(dict(item=e.name, ledger=[e.qubits, e.cycles_per_attempt, e.acceptance, e.output_error],
                         cudaq=[c.physical_qubits, c.cycles_per_attempt, c.acceptance_probability,
                                c.output_infidelity]))
    for r in rows:
        r['match'] = all(isclose(a, b, rel_tol=.02) for a, b in zip(r['ledger'], r['cudaq']))
    return rows


def device(preset, n_logical, engine):
    arch = getattr(P, preset)
    block = P.processing_block(arch)
    b = cql.devices.DeviceBuilder('PinnacleCheck')
    compute = b.logical.add_compute(capacity=n_logical, name='compute')
    supply = b.logical.add_factory(produces=cql.logical.T_STATE, via=engine.producer, capacity=1, buffer_size=1,
                                   name='engine')
    cq = b.qec.bind(compute, architecture=arch)
    fq = b.qec.bind(supply, encoding=arch.encoding)
    opts = dict(granularity=cql.architecture.ResourceGranularity.PATCH,
                footprint=cql.architecture.PhysicalFootprint('qubit', block.physical_qubits,
                                                             'Pinnacle processing block (Table I)'),
                native_actions=cql.architecture.physical_actions.clifford_set(),
                native_instruments=(cql.architecture.physical_instruments.MPP,),
                capabilities=(cql.architecture.NATIVE_PAULI_PRODUCT_ROTATION,))
    blocks = ceil(n_logical/block.logical_capacity)
    data = b.physical.add_resources('gb_block', blocks, name='blocks', **opts)
    scratch = b.physical.add_resources('gb_block', 1, name='rus_scratch', **opts)
    lane = b.physical.add_resources('qubit', engine.physical_qubits, name='engine_qubits',
                                    granularity=cql.architecture.ResourceGranularity.CARRIER)
    b.physical.bind(cq, to=data)
    (sq,) = cq.auxiliary_regions
    b.physical.bind(sq, to=scratch)
    b.physical.bind(fq, to=lane, factory_model=cql.devices.FactoryModel(
        startup_cycles=engine.cycles_per_attempt, output_interval_cycles=engine.cycles_per_attempt,
        evidence=cql.analysis.user_assertion('Pinnacle engine attempt cadence')))
    d = P._PINNACLE_GB_INSTANCES[preset].distance
    b.physical.set_operating_point(timing={'cycle_ns': CYCLE_NS, 'surface_cycle_ns': CYCLE_NS,
                                           'rpp_ns': (d+2)*CYCLE_NS, 'mpp_ns': (d+2)*CYCLE_NS,
                                           'condition_ns': 200.0})
    footprint = dict(blocks=blocks, block_qubits=block.physical_qubits, scratch_blocks=1,
                     engine_qubits=engine.physical_qubits,
                     cudaq_total=(blocks+1)*block.physical_qubits+engine.physical_qubits,
                     ledger_total=ledger.unit_qubits(next(c for c in ledger.CODES if c.name == preset),
                                                     n_logical, next(e for e in ledger.ENGINES
                                                                     if e.qubits == engine.physical_qubits)))
    return b.build(), footprint


def first_line(exc):
    text = f'{type(exc).__name__}: {exc}'
    lines = text.splitlines()
    detail = next((l.strip() for l in lines if "' op " in l), None)
    return (lines[0] + (' | ' + detail if detail else ''))[:500]


def compile_route(route, dev, prog):
    out = dict(route=route)
    try:
        if route == 'A':
            full = cql.compiler.synthesize(prog, gate_set=cql.compiler.clifford_t, precision=PRECISION)
            actions = cql.estimate(full, tier=cql.estimate.Tier.LOGICAL).to_dict()['actions']
            out['T_count'] = actions.get('qlx_standard_t', 0) + actions.get('qlx_standard_tdg', 0)
            p0 = cql.compiler.to_pbc(full)
        else:
            p0 = cql.compile(prog)
        out['p0_actions'] = cql.estimate(p0, tier=cql.estimate.Tier.LOGICAL).to_dict()['actions']
        p1 = cql.compiler.place(p0, device=dev)
        p2 = cql.compile(p1, pipeline=cql.compiler.pipelines.qec(), device=dev,
                         policy={'rpp_precision': PRECISION})
        static = cql.estimate(p2, tier=cql.estimate.Tier.STATIC).to_dict()
        ops = static['operation_counts']
        calls = Counter()
        for name, n in static['gadget_calls'].items():
            calls['inject' if '_inject' in name else name.rsplit('_', 1)[-1] if 'pinnacle' in name else name] += n
        out.update(status='P2 compiled', engine_requests=ops.get('resource_request', 0),
                   resource_rotations=ops.get('resource_rotate_product', 0), gadget_calls=dict(calls))
    except Exception as exc:
        out.update(status='failed before P2', error=first_line(exc))
        return out, None
    for tier in ('ANALYTICAL',):
        try:
            cql.estimate(p2, tier=getattr(cql.estimate.Tier, tier), p_phys=P_PHYS, failure_budget=.1,
                         cycle_time=CYCLE_NS*1e-9)
            out[tier.lower()] = 'ok'
        except Exception as exc:
            out[tier.lower()] = first_line(exc)
    try:
        cql.compile(p2, pipeline=cql.compiler.pipelines.physical(), device=dev)
        out['p3'] = 'ok'
    except Exception as exc:
        out['p3'] = first_line(exc)
    return out, p2


def main():
    sys.setrecursionlimit(50000)
    OUT.mkdir(parents=True, exist_ok=True)
    rows = tables()
    engine = P.magic_engine(p_phys=P_PHYS, target_output_infidelity=1e-9)
    ops, _ = build_workload(2, .5, 1, 'plaquette')
    dev, footprint = device('gb30', 8, engine)
    prog = program_from_operations(ops, 8, 'hubbard_2x2')
    routes = [compile_route(r, dev, prog)[0] for r in ('A', 'B')]
    result = dict(cudaq_logical_version=getattr(cql, '__version__', 'unknown'),
                  workload='2x2 Hubbard, plaquette, one step, tau = 0.5', precision=PRECISION,
                  rotations=sum(o.name == 'rotate' for o in ops), tables=rows,
                  tables_match=all(r['match'] for r in rows), footprint=footprint, routes=routes,
                  notes=["P3 failure diagnostic (native pass, printed to stderr): 'fabric.resource_rotate_product' op "
                         "has an invalid resource-assisted physical Pauli product; with carrier-level qubits instead: "
                         "'fabric.transform_begin' op has an invalid physical frame projection",
                         "Analytical tier: the WSC gadgets' fabric.output_syndrome ops have no typed cost semantics",
                         "Route B needs pygridsynth (installed 2.0.0 into .venv-cudaq)"])
    (OUT/'result.json').write_text(json.dumps(result, indent=2, default=str)+'\n')
    print(json.dumps(result, indent=1, default=str))


if __name__ == '__main__':
    main()
