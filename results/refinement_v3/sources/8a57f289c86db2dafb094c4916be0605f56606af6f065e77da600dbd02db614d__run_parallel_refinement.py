"""Compile selected stage-parallel plaquette circuits and retain alternatives."""
from dataclasses import asdict
import gzip
import hashlib
import json
from pathlib import Path

from hubbard_circuits import build_workload
from model import workload
from nisq_compilation import cudaq_operations, fermionic_blocks
from nisq_observables import parity_acceptance
from parallel_nisq import observable_stages, parallel_route, plaquette_stages
from run_nisq_refinement import describe

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/refinement_v3'


def run(L,prune=True):
    c=json.loads((ROOT/'assumptions.json').read_text());w=workload(L,1.,4,'plaquette',c['u_over_j'])
    ops,observed=build_workload(L,1.,4);prep,_=fermionic_blocks(ops)
    stages=plaquette_stages(L,1.,4);modes=list(range(2*L*L))
    if prune:stages,modes=observable_stages(stages,observed)
    tileorder=[m for group in stages[0].groups for m in sorted({m for b in group for m in b.modes})]
    tileorder += [m for m in modes if m not in tileorder]
    choices=[]
    for label,order in [('site',modes),('spin',sorted(modes,key=lambda m:(m%2,m))),('tile',tileorder)]:
        for inside in (False,True):
            gates,mapping,info=parallel_route(stages,order,prep,'mean',inside)
            record,events=describe(gates,mapping,observed,c,w)
            record.update(info);record['initial_order_policy']=label
            score=record['costs_by_idle_ratio']['1.0']['log10_Q_min']
            choices.append((score,record['observable_exposure']['chi'],gates,events,record))
    _,_,gates,events,record=min(choices,key=lambda v:v[:2])
    record['candidates']=[{'initial_order':v[4]['initial_order_policy'],'optimize_inside':v[4]['optimize_inside'],
        'G2':v[4]['schedule']['two_qubit_gates'],'depth':v[4]['schedule']['depth_layers'],
        'chi':v[4]['observable_exposure']['chi'],'log10_Q':v[0]} for v in choices]
    record.update({'L':L,'tau':1.,'steps':4,'algorithm':'plaquette','observable_pruned':prune,
                   'configuration':c,'observable_modes':observed})
    record['postselection']=parity_acceptance(events,record['schedule'],record['mode_to_wire'],c)
    import cudaq.logical as cql
    from cudaq_compile import logical_counts,program_from_operations
    p0=cql.compile(program_from_operations(cudaq_operations(gates),len(modes),'parallel_hubbard'))
    record['cudaq_logical']=logical_counts(p0)
    assert record['cudaq_logical']['actions']['qlx_standard_pauli_rotation']==sum(g.name in ('rz','rx','rzz') for g in gates)
    directory=OUT/f"parallel_L{L}_{'observable' if prune else 'full'}"
    directory.mkdir(parents=True,exist_ok=True)
    with gzip.open(directory/'native_p0.mlir.gz','wt') as f:f.write(p0.to_mlir())
    with gzip.open(directory/'native_gates.json.gz','wt') as f:json.dump([asdict(g) for g in gates],f)
    record['sources']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
                       ('parallel_nisq.py','run_parallel_refinement.py','nisq_compilation.py','nisq_observables.py',
                        'hubbard_circuits.py','run_nisq_refinement.py','model.py')}
    (directory/'result.json').write_text(json.dumps(record,indent=2)+'\n')
    print(directory.name,'G2',record['schedule']['two_qubit_gates'],'depth',record['schedule']['depth_layers'],
          'chi',record['observable_exposure']['chi'],'log10 Q',record['costs_by_idle_ratio']['1.0']['log10_Q_min'],flush=True)


if __name__=='__main__':
    for size in (2,4,6,8):run(size)
    run(4,False)
