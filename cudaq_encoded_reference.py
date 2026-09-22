"""P0-to-P3 FT integration using explicit surface-patch and factory macros.

This checks the compiler integration, not physical NN lattice-surgery timing.
Native multi-Pauli operations and factory service are declared model inputs.
"""
import argparse
from dataclasses import asdict, replace
from inspect import Parameter, Signature
import json
from pathlib import Path
import sys
from math import pi

import cudaq.logical as cql

from cudaq_compile import program_from_operations
from hubbard_circuits import Operation, angle_summary, build_workload
from model import accuracy, factory


@cql.protocol(implements=cql.logical.produce(cql.logical.T_STATE))
def t_state_source() -> cql.types.resource[cql.logical.T_STATE]:
    return cql.produce(cql.logical.T_STATE)


def macro_mpp_lowering(encoding):
    """Explicit native encoded-MPP assumption for cross-patch final readout.

    Pinnacle's supplied WSC measurement only accepts a single packed block.
    This adapter preserves the Pauli product on separate surface patches;
    it does not construct a surgery circuit or certify its distance.
    """
    def compile_site(site,context):
        placements=tuple(context.placements)
        assert len({p.block or p.placement for p in placements})==len(placements)
        assert all(p.logical_index==0 for p in placements)
        x,z,sign=(site.parameters[k] for k in ('x_mask','z_mask','sign'))
        assert x|z==(1<<len(placements))-1
        annotation=cql.types.patch[context.encoding]
        def generated(*blocks):
            product=None
            for i,block in enumerate(blocks):
                letter='Y' if (x>>i)&1 and (z>>i)&1 else ('X' if (x>>i)&1 else 'Z')
                term=getattr(cql.types,letter)(block[0])
                product=term if product is None else product@term
            return cql.ops.mpp(product if sign>0 else -product)
        generated.__name__=f'hubbard_macro_mpp_x{x}_z{z}_s{sign}'
        generated.__qualname__=generated.__name__
        parameters=tuple(Parameter(f'block{i}',Parameter.POSITIONAL_OR_KEYWORD,
                                   annotation=annotation) for i in range(len(placements)))
        result_type=tuple[tuple([annotation]*len(placements)+[bool])]
        generated.__signature__=Signature(parameters,return_annotation=result_type)
        gadget=cql.gadgets.GadgetDefinition(generated,implements=cql.logical.mpp,
            name=generated.__name__,
            type_hints={**{p.name:annotation for p in parameters},'return':result_type},
            metadata={'construction':'assumed native encoded multi-Pauli measurement',
                      'distance_status':'unvalidated','geometry':'unspecified'})
        return cql.qec.GeneratedQECArtifact(gadget,dict(site.parameters))
    return cql.qec.QECLowering(compile_site,objective_family='pauli_product_measurement',
        objective=cql.logical.mpp,codes=(encoding,),name='hubbard_native_mpp_macro',
        plugin='architecture_comparison',version='1',
        metadata={'status':'physical macro assumption, not gate-level surgery'})


def run(L=2,steps=4,distance=3,probe_rotations=0):
    # CUDA-Q 0.1.1's ownership verifier recursively traverses SSA chains.
    # Thousands of synthesized T gates exceed Python's default depth of 1000.
    sys.setrecursionlimit(max(sys.getrecursionlimit(),50000))
    c=json.loads(Path('assumptions.json').read_text())
    ops,_=build_workload(L,1.,steps,'plaquette')
    if probe_rotations:
        ops=[Operation('rotate',(0,),pi/8,'X' if i%2 else 'Z','fixed_T')
             for i in range(probe_rotations)]
    eps,_=accuracy(c)
    precision=eps*c['synthesis_fraction']/(2*max(1,angle_summary(ops)['dynamic_rotations']))
    surface=cql.architectures.surface.definitions(distance)
    architecture=cql.architectures.pinnacle.for_code(surface.code,cycle_rounds=distance)
    architecture=replace(architecture,name='HubbardSurfaceMacro',link_roots=tuple(
        root for root in architecture.link_roots
        if getattr(root,'objective_family',None)!='pauli_product_measurement'
    )+(macro_mpp_lowering(architecture.encoding),))
    builder=cql.devices.DeviceBuilder('HubbardSurfaceMacroReference')
    compute=builder.logical.add_compute(capacity=2*L*L,name='compute')
    supply=builder.logical.add_factory(produces=cql.logical.T_STATE,via=t_state_source,
                                       capacity=1,buffer_size=1,name='t_factory')
    compute_qec=builder.qec.bind(compute,architecture=architecture)
    factory_qec=builder.qec.bind(supply,encoding=surface.code)
    options={'granularity':cql.architecture.ResourceGranularity.PATCH,
             'footprint':surface.square_patch_footprint,
             'native_actions':cql.architecture.physical_actions.clifford_set(),
             'native_instruments':(cql.architecture.physical_instruments.MPP,),
             'capabilities':(cql.architecture.NATIVE_PAULI_PRODUCT_ROTATION,)}
    data=builder.physical.add_resources('surface_patch',2*L*L,name='data',**options)
    scratch=builder.physical.add_resources('surface_patch',1,name='scratch',**options)
    tiles,clocks,error=factory(2,c['p'])
    lane=builder.physical.add_resources('factory_lane',1,name='factory',
          granularity=cql.architecture.ResourceGranularity.PATCH,
          footprint=cql.architecture.PhysicalFootprint('qubit',tiles*surface.square_patch_footprint.units,
                                                      'declared serial two-level 15-to-1 reference'))
    builder.physical.bind(compute_qec,to=data)
    (scratch_qec,)=compute_qec.auxiliary_regions
    builder.physical.bind(scratch_qec,to=scratch)
    builder.physical.bind(factory_qec,to=lane,factory_model=cql.devices.FactoryModel(
          startup_cycles=clocks*distance,output_interval_cycles=clocks*distance,
          evidence=cql.analysis.user_assertion('serial two-level 15-to-1 mean service; analytical_forms.tex')))
    cycle=c['cycle_seconds']*1e9
    builder.physical.set_operating_point(timing={
        'cycle_ns':cycle,'surface_cycle_ns':cycle,'rpp_ns':9*distance*cycle,
        'mpp_ns':distance*cycle,'condition_ns':c['feedforward_seconds']*1e9})
    device=builder.build()
    print('Synthesizing FT workload',flush=True)
    full=cql.compiler.synthesize(program_from_operations(ops,2*L*L,'hubbard'),
                                 gate_set=cql.compiler.clifford_t,precision=precision)
    p0=cql.compiler.to_pbc(full)
    p1=cql.compiler.place(p0,device=device)
    print('Selecting encoded gadgets and factory',flush=True)
    p2=cql.compile(p1,pipeline=cql.compiler.pipelines.qec(),device=device)
    print('Lowering and scheduling physical macros',flush=True)
    p3=cql.compile(p2,pipeline=cql.compiler.pipelines.physical(),device=device)
    schedule=cql.compiler.schedule(p3)
    try:
        result=cql.estimate(schedule,tier=cql.estimate.Tier.SCHEDULE,p_phys=c['p'],
                             failure_budget=eps*c['logical_fraction'],cycle_time=c['cycle_seconds'])
        estimate_record=result.to_dict()
    except Exception as exc:
        estimate_record={'status':'REJECTED by CUDA-Q schedule estimator',
                         'exception':type(exc).__name__,'message':str(exc),
                         'cause':str(exc.__cause__), 'makespan_ns':schedule.makespan_ns}
    counts=cql.estimate(full,tier=cql.estimate.Tier.LOGICAL).to_dict()
    nt=counts['actions'].get('qlx_standard_t',0)+counts['actions'].get('qlx_standard_tdg',0)
    # Deliberately loose serial sanity bound for this macro model, including
    # factory startup for EVERY T and ample time for preparation and readout.
    sanity_ns=(nt+2*L*L+1)*(2*clocks+100)*distance*cycle
    sane=schedule.makespan_ns<=sanity_ns
    record={'L':L,'r':steps,'distance':distance,'precision':precision,
            'probe_rotations':probe_rotations,'T_count':nt,
            'stages':[str(b.stage) for b in (p0,p1,p2,p3)],'schedule':estimate_record,
            'status':(estimate_record.get('status') or ('macro schedule passes timing sanity check' if sane else
                      'REJECTED timing: scheduler exceeds a loose serial upper bound')),
            'scope':'native encoded MPP/RPP macro reference, not nearest-neighbor gate-level QEC',
            'timing_sanity_bound_ns':sanity_ns,'timing_sanity_passed':sane,
            'factory_output_injection_error':error,
            'factory_service_clocks':clocks,'factory_patch_tiles':tiles,
            'patch_footprint':surface.square_patch_footprint.units,
            'rpp_ns':9*distance*cycle}
    name='encoded_reference' if (L,steps,distance)==(2,4,3) else f'encoded_reference_L{L}_r{steps}_d{distance}'
    out=Path('results/cudaq_v1')/(name if not probe_rotations else f'factory_probe_{probe_rotations}')
    out.mkdir(parents=True,exist_ok=True)
    (out/'result.json').write_text(json.dumps(record,indent=2)+'\n')
    (out/'physical_p3.mlir').write_text(p3.to_mlir())
    (out/'schedule.json').write_text(json.dumps([asdict(entry) for entry in schedule.entries],indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--L',type=int,default=2)
    parser.add_argument('--steps',type=int,default=4)
    parser.add_argument('--distance',type=int,default=3)
    parser.add_argument('--probe-rotations',type=int,default=0)
    args=parser.parse_args()
    run(args.L,args.steps,args.distance,args.probe_rotations)
