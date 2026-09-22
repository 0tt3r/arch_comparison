"""Compile native RZZ circuits; quantify local PEC and parity postselection."""

from project_paths import ROOT, project_file
import argparse
from dataclasses import asdict
import gzip
import hashlib
import importlib.metadata
import json
from math import log
from pathlib import Path

from hubbard_circuits import build_workload
from model import finish, log_gamma2_dephasing, log_gamma2_depolarizing, workload
from nisq_compilation import (causal_cone, cudaq_operations, fermionic_blocks,
                              fused_route, grid_parity_route, schedule)
from nisq_observables import idle_probability, observable_exposure, parity_acceptance

OUT=ROOT/'results/nisq_v2'
LOOKAHEADS=(0,8,16,32)


def full_exposure(sched,c):
    return (sched['two_qubit_gates']*log_gamma2_depolarizing(c['p'],2)
            +sched['one_qubit_gates']*log_gamma2_depolarizing(c['p']/10,1)
            +sched['idle_qubit_layers']*log_gamma2_depolarizing(idle_probability(c),1)
            +(sched['width']+2)*log_gamma2_dephasing(c['p']))


def describe(gates,mapping,observed,c,w):
    sched,events=schedule(gates,len(mapping),gate_ns=c['gate_seconds']*1e9,
                          reset_ns=c['reset_seconds']*1e9,measurement_ns=c['measurement_seconds']*1e9)
    outputs=[mapping[q] for q in observed]
    exposure=observable_exposure(events,sched,outputs,c)
    costs={}
    for ratio in (0.,.1,1.):
        scenario={**c,'idle_p_per_100ns':ratio*c['p']}
        exp=observable_exposure(events,sched,outputs,scenario)
        costs[str(ratio)]=finish(w,'physical_local_PEC',exp['chi'],sched['physical_copy_qubits'],
                                 sched['seconds_per_attempt'],scenario)
    return {'schedule':sched,'observable_exposure':exposure,'full_circuit_chi':full_exposure(sched,c),
            'mode_to_wire':mapping,'observable_wires':outputs,'costs_by_idle_ratio':costs},events


def best_route(blocks,modes,prep,observed,c,w):
    candidates=[]
    for lookahead in LOOKAHEADS:
        gates,mapping,info=fused_route(blocks,modes,prep,lookahead)
        record,events=describe(gates,mapping,observed,c,w)
        record.update(info);record['lookahead']=lookahead
        row=record['costs_by_idle_ratio']['1.0']
        score=(row['log10_Q_min'],record['observable_exposure']['chi']+log(row['seconds_per_attempt']))
        candidates.append((score,gates,events,record))
    _,gates,events,record=min(candidates,key=lambda item:item[0])
    record['candidate_metrics']=[{'lookahead':row['lookahead'],
        'two_qubit_gates':row['schedule']['two_qubit_gates'],
        'depth_layers':row['schedule']['depth_layers'],'chi':row['observable_exposure']['chi'],
        'log10_Q_min':row['costs_by_idle_ratio']['1.0']['log10_Q_min']}
        for _,_,_,row in candidates]
    return gates,events,record


def compile_case(L,steps,algorithm,c,compile_cudaq=True,variants=('fused_full','fused_cone','grid_tree')):
    ops,observed=build_workload(L,1.,steps,algorithm,c['u_over_j'])
    prep,blocks=fermionic_blocks(ops)
    kept,modes=causal_cone(blocks,observed)
    w=workload(L,1.,steps,algorithm,c['u_over_j'])
    records=[]
    for name in variants:
        print(f'{algorithm} L={L} r={steps}: {name}',flush=True)
        if name=='grid_tree':
            gates,info=grid_parity_route(ops,2*L*L)
            record,events=describe(gates,dict(enumerate(range(2*L*L))),observed,c,w)
            record.update(info)
        else:
            selected,selected_modes=(blocks,list(range(2*L*L))) if name=='fused_full' else (kept,modes)
            gates,events,record=best_route(selected,selected_modes,prep,observed,c,w)
        record.update({'L':L,'tau':1.,'steps':steps,'algorithm':algorithm,'variant':name,
                       'configuration':c,'blocks_full':len(blocks),'blocks_in_fermionic_cone':len(kept),
                       'modes_in_fermionic_cone':len(modes),
                       'native_assumption':'arbitrary RZZ and 1q pulses: 10ns each; p2=p, p1=p/10'})
        record['postselection']=parity_acceptance(events,record['schedule'],record['mode_to_wire'],c)
        directory=OUT/f'{algorithm}_L{L}_r{steps}_{name}'
        directory.mkdir(parents=True,exist_ok=True)
        if compile_cudaq:
            import cudaq.logical as cql
            from cudaq_compile import program_from_operations,logical_counts
            built=cql.compile(program_from_operations(cudaq_operations(gates),record['schedule']['width'],'native_rzz'))
            record['cudaq_logical']=logical_counts(built)
            # Rz, Rx and RZZ use the same QLX rotation action, with explicit support.
            rotations=sum(g.name in ('rz','rx','rzz') for g in gates)
            assert record['cudaq_logical']['actions'].get('qlx_standard_pauli_rotation',0)==rotations
            with gzip.open(directory/'native_p0.mlir.gz','wt') as handle:handle.write(built.to_mlir())
        with gzip.open(directory/'native_gates.json.gz','wt') as handle:
            json.dump([asdict(g) for g in gates],handle)
        record['sources']={name:hashlib.sha256((project_file(name)).read_bytes()).hexdigest()
                           for name in ('nisq_compilation.py','nisq_observables.py','run_nisq_refinement.py',
                                        'hubbard_circuits.py','cudaq_compile.py')}
        (directory/'result.json').write_text(json.dumps(record,indent=2)+'\n')
        print('  native entanglers',record['schedule']['two_qubit_gates'],
              'chi',round(record['observable_exposure']['chi'],3),flush=True)
        records.append(record)
    return records


def cone_sweep():
    rows=[]
    for algorithm in ('plaquette','termwise'):
        for steps in (1,2,4):
            for L in (2,4,6,8,12,16,20,24):
                ops,observed=build_workload(L,1.,steps,algorithm)
                _,blocks=fermionic_blocks(ops);kept,modes=causal_cone(blocks,observed)
                rows.append({'algorithm':algorithm,'L':L,'r':steps,'full_modes':2*L*L,
                             'cone_modes':len(modes),'full_blocks':len(blocks),'cone_blocks':len(kept)})
    (OUT/'cone_sweep.json').write_text(json.dumps(rows,indent=2)+'\n')


def manifest():
    files=[project_file(name) for name in ('nisq_compilation.py','nisq_observables.py','run_nisq_refinement.py',
           'hubbard_circuits.py','cudaq_compile.py','model.py','compiled_model.py',
           'boundary_timing.py','factory_sensitivity.py','plot_nisq_refinement.py',
           'test_model.py','test_hubbard_circuits.py','test_compiled_model.py','test_nisq_refinement.py',
           'nisq_refinement.tex','nisq_refinement.pdf','references.bib','refinement_references.bib',
           'NISQ_REFINEMENT.md','README.md','WORKING_DECISIONS.md','assumptions.json',
           'requirements-cudaq.lock')]
    snapshots=OUT/'sources';snapshots.mkdir(exist_ok=True)
    for path in files:
        if path.exists():
            data=path.read_bytes();sha=hashlib.sha256(data).hexdigest()
            (snapshots/(sha+'__'+path.name)).write_bytes(data)
    for path in OUT.glob('*_L*_r*_*/result.json'):
        for name,sha in json.loads(path.read_text())['sources'].items():
            saved=snapshots/(sha+'__'+name)
            if not saved.exists() or hashlib.sha256(saved.read_bytes()).hexdigest()!=sha:
                raise ValueError(f'Missing original source snapshot for {path.parent.name}: {name} {sha}')
    packages={}
    for name in ('cudaq-logical','cuda-quantum-cu12','numpy','scipy','matplotlib'):
        try:packages[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:pass
    record={'packages_in_compiler_environment':packages,
            'sources':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.exists()},
            'artifacts':{str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='manifest.json'}}
    (OUT/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes',type=int,nargs='+',default=[2,4,6,8])
    parser.add_argument('--steps',type=int,default=4)
    parser.add_argument('--algorithms',nargs='+',default=['plaquette','termwise'])
    parser.add_argument('--variants',nargs='+',choices=['fused_full','fused_cone','grid_tree'],
                        default=['fused_full','fused_cone','grid_tree'])
    parser.add_argument('--counts-only',action='store_true')
    parser.add_argument('--cones-only',action='store_true')
    parser.add_argument('--manifest-only',action='store_true')
    args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if not args.manifest_only:
        c=json.loads((project_file('assumptions.json')).read_text())
        if args.cones_only:cone_sweep()
        else:
            for L in args.sizes:
                for algorithm in args.algorithms:
                    compile_case(L,args.steps,algorithm,c,not args.counts_only,args.variants)
    manifest()
