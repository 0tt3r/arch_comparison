"""Standalone reproducer: `phys-schedule` ScheduleConflict for a small fixed-angle program.

A three-qubit pattern of four pi/8 Pauli rotations (ZZ, X, ZZ, X) is repeated n times,
synthesized to Clifford+T, placed on a surface-code device with one T-state factory, and
compiled through the QEC and physical pipelines, then scheduled. With cudaq-logical
0.1.1.post1 some repetition counts schedule and others are rejected with

    'phys.schedule' op nested schedule entry must lie within its parent interval

Usage:  python repro.py [n ...]        (default: 1 2 3 4 5 6 8)

Everything the original project read from its configuration files is inlined below.
"""
import json
import sys
import traceback
from inspect import Parameter, Signature
from dataclasses import replace
from math import pi
from pathlib import Path

import cudaq.logical as cql
from cudaq.logical.algebra import X, Z

OUT = Path(__file__).resolve().parent/'output'

# Inlined constants (from the project's assumptions.json and factory model).
EPS = 1e-3                 # accuracy target
SYNTHESIS_FRACTION = 0.075
P_PHYS = 1e-3
CYCLE_S = 240e-9           # syndrome cycle
FEEDFORWARD_S = 200e-9
LOGICAL_FRACTION = 0.1
DISTANCE = 3
WIDTH = 8                  # logical data qubits allocated (only 0..2 are used)
FACTORY_TILES = 37         # declared two-level 15-to-1 footprint, in patches
FACTORY_CLOCKS = 178.5127840745416   # declared service time per T, in units of d cycles

# The pattern: exp(-i pi/8 P) for P = Z0Z1, X1, Z1Z2, X2 (cql.ops.rotate takes 2*theta).
PATTERN = [((0, 1), 'ZZ'), ((1,), 'X'), ((1, 2), 'ZZ'), ((2,), 'X')]


@cql.protocol(implements=cql.logical.produce(cql.logical.T_STATE))
def t_state_source() -> cql.types.resource[cql.logical.T_STATE]:
    return cql.produce(cql.logical.T_STATE)


def macro_mpp_lowering(encoding):
    """Native encoded multi-Pauli measurement across separate surface patches.

    The supplied pinnacle Pauli-product-measurement lowering accepted only a single packed
    block, so the project replaced it with this macro (see README, question 2).
    """
    def compile_site(site, context):
        placements = tuple(context.placements)
        x, z, sign = (site.parameters[k] for k in ('x_mask', 'z_mask', 'sign'))
        annotation = cql.types.patch[context.encoding]

        def generated(*blocks):
            product = None
            for i, block in enumerate(blocks):
                letter = 'Y' if (x >> i) & 1 and (z >> i) & 1 else ('X' if (x >> i) & 1 else 'Z')
                term = getattr(cql.types, letter)(block[0])
                product = term if product is None else product @ term
            return cql.ops.mpp(product if sign > 0 else -product)
        generated.__name__ = generated.__qualname__ = f'macro_mpp_x{x}_z{z}_s{sign}'
        params = tuple(Parameter(f'block{i}', Parameter.POSITIONAL_OR_KEYWORD, annotation=annotation)
                       for i in range(len(placements)))
        result_type = tuple[tuple([annotation]*len(placements) + [bool])]
        generated.__signature__ = Signature(params, return_annotation=result_type)
        gadget = cql.gadgets.GadgetDefinition(
            generated, implements=cql.logical.mpp, name=generated.__name__,
            type_hints={**{p.name: annotation for p in params}, 'return': result_type},
            metadata={'construction': 'assumed native encoded multi-Pauli measurement'})
        return cql.qec.GeneratedQECArtifact(gadget, dict(site.parameters))
    return cql.qec.QECLowering(compile_site, objective_family='pauli_product_measurement',
                               objective=cql.logical.mpp, codes=(encoding,), name='native_mpp_macro',
                               plugin='bug_report', version='1', metadata={})


def program(repetitions):
    @cql.program
    def prog() -> None:
        q = cql.allocate(WIDTH, state=cql.types.zero)
        for _ in range(repetitions):
            for qubits, letters in PATTERN:
                product = None
                for i, letter in zip(qubits, letters):
                    factor = {'X': X, 'Z': Z}[letter](q[i])
                    product = factor if product is None else product @ factor
                values = cql.ops.rotate(product, angle=2*(pi/8))
                after = {id(b): a for b, a in zip(product.operands, values)}
                for i, before in zip(qubits, tuple(q[i] for i in qubits)):
                    q[i] = after[id(before)]
        for i in range(WIDTH):
            cql.measure_z(q[i])
    return prog


def device():
    surface = cql.architectures.surface.definitions(DISTANCE)
    arch = cql.architectures.pinnacle.for_code(surface.code, cycle_rounds=DISTANCE)
    arch = replace(arch, name='SurfaceMacro', link_roots=tuple(
        r for r in arch.link_roots if getattr(r, 'objective_family', None) != 'pauli_product_measurement'
    ) + (macro_mpp_lowering(arch.encoding),))
    b = cql.devices.DeviceBuilder('SurfaceMacroReference')
    compute = b.logical.add_compute(capacity=WIDTH, name='compute')
    supply = b.logical.add_factory(produces=cql.logical.T_STATE, via=t_state_source,
                                   capacity=1, buffer_size=1, name='t_factory')
    compute_qec = b.qec.bind(compute, architecture=arch)
    factory_qec = b.qec.bind(supply, encoding=surface.code)
    opts = {'granularity': cql.architecture.ResourceGranularity.PATCH,
            'footprint': surface.square_patch_footprint,
            'native_actions': cql.architecture.physical_actions.clifford_set(),
            'native_instruments': (cql.architecture.physical_instruments.MPP,),
            'capabilities': (cql.architecture.NATIVE_PAULI_PRODUCT_ROTATION,)}
    data = b.physical.add_resources('surface_patch', WIDTH, name='data', **opts)
    scratch = b.physical.add_resources('surface_patch', 1, name='scratch', **opts)
    lane = b.physical.add_resources(
        'factory_lane', 1, name='factory', granularity=cql.architecture.ResourceGranularity.PATCH,
        footprint=cql.architecture.PhysicalFootprint('qubit', FACTORY_TILES*surface.square_patch_footprint.units,
                                                     'declared factory footprint'))
    b.physical.bind(compute_qec, to=data)
    (scratch_qec,) = compute_qec.auxiliary_regions
    b.physical.bind(scratch_qec, to=scratch)
    b.physical.bind(factory_qec, to=lane, factory_model=cql.devices.FactoryModel(
        startup_cycles=FACTORY_CLOCKS*DISTANCE, output_interval_cycles=FACTORY_CLOCKS*DISTANCE,
        evidence=cql.analysis.user_assertion('declared serial factory service time')))
    cycle = CYCLE_S*1e9
    b.physical.set_operating_point(timing={'cycle_ns': cycle, 'surface_cycle_ns': cycle,
                                           'rpp_ns': 9*DISTANCE*cycle, 'mpp_ns': DISTANCE*cycle,
                                           'condition_ns': FEEDFORWARD_S*1e9})
    return b.build()


def run(repetitions):
    out = OUT/f'repetitions_{repetitions}'
    out.mkdir(parents=True, exist_ok=True)
    stage, n_t = 'synthesize', None
    try:
        dev = device()
        full = cql.compiler.synthesize(program(repetitions), gate_set=cql.compiler.clifford_t,
                                       precision=EPS*SYNTHESIS_FRACTION/2)
        counts = cql.estimate(full, tier=cql.estimate.Tier.LOGICAL).to_dict()['actions']
        n_t = counts.get('qlx_standard_t', 0) + counts.get('qlx_standard_tdg', 0)
        stage = 'place'
        p1 = cql.compiler.place(cql.compiler.to_pbc(full), device=dev)
        stage = 'qec pipeline'
        p2 = cql.compile(p1, pipeline=cql.compiler.pipelines.qec(), device=dev)
        stage = 'physical pipeline'
        p3 = cql.compile(p2, pipeline=cql.compiler.pipelines.physical(), device=dev)
        (out/'physical_p3.mlir').write_text(p3.to_mlir())
        stage = 'schedule'
        schedule = cql.compiler.schedule(p3)
        result = dict(repetitions=repetitions, T_count=n_t, status='scheduled',
                      makespan_ns=schedule.makespan_ns)
    except Exception as exc:
        (out/'traceback.txt').write_text(traceback.format_exc() + f'\nCAUSE:\n{exc.__cause__}\n')
        result = dict(repetitions=repetitions, T_count=n_t, status=f'FAILED at {stage}', exception=type(exc).__name__,
                      message=str(exc).splitlines()[0][:300],
                      diagnostic=next((l.strip() for l in (str(exc)+'\n'+str(exc.__cause__)).splitlines()
                                       if 'must lie within' in l), None))
    (out/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    # cudaq-logical 0.1.1's ownership verifier recurses along SSA chains.
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 50000))
    reps = [int(a) for a in sys.argv[1:]] or [1, 2, 3, 4, 5, 6, 8]
    results = [run(n) for n in reps]
    (OUT/'summary.json').write_text(json.dumps(results, indent=2)+'\n')
