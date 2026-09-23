"""Compile an onsite stage, including clean-AND arithmetic and feedforward."""
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path

import cudaq.logical as cql

from cudaq_compile import logical_counts,t_count,program_from_operations
from hubbard_circuits import Operation
from hwp_stage import stage_operations,nearest_neighbor_route,schedule_patch_gates
from project_paths import ROOT

OUT=ROOT/'results/implementation_v1/ft_stage'


def program(ops,width):
    @cql.program
    def stage() -> None:
        qubits=list(cql.allocate(width,state=cql.types.zero))
        for op in ops:
            if op.name=='unand':
                a,b,t=op.qubits
                result=cql.measure_x(qubits[t]);qubits[t]=None
                qubits[a],qubits[b]=cql.ops.cond(result,then=lambda x,y:cql.cz(x,y),
                    else_=lambda x,y:(x,y),carries=(qubits[a],qubits[b]))
            else:
                for q in op.qubits:
                    if qubits[q] is None:qubits[q]=cql.allocate(1,state=cql.types.zero)[0]
                if op.name=='cx':
                    a,b=op.qubits;qubits[a],qubits[b]=cql.cx(qubits[a],qubits[b])
                elif op.name=='rz':
                    q=op.qubits[0];qubits[q]=cql.rz(qubits[q],op.angle)
                else:
                    q=op.qubits[0];qubits[q]=getattr(cql,op.name)(qubits[q])
        for q in qubits:
            if q is not None:cql.measure_z(q)
    return stage


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    c=json.loads((ROOT/'configs/implementation_checks.json').read_text())['ft'];rows=[]
    for batch in c['batch_sizes']:
        ops,width=stage_operations(c['L']**2,batch,c['theta'])
        analog=cql.compile(program(ops,width));full=cql.compiler.synthesize(analog,gate_set=cql.compiler.clifford_t,precision=c['synthesis_precision'])
        counts=logical_counts(full)
        for label,build in [('logical',analog),('synthesized',full)]:
            with gzip.open(OUT/f'batch{batch}_{label}.mlir.gz','wt') as f:f.write(build.to_mlir())
        angle_cost={}
        for op in ops:
            if op.name=='rz' and op.angle.hex() not in angle_cost:
                build=cql.compiler.synthesize(program_from_operations([Operation('rz',(0,),op.angle)],1,'rotation'),gate_set=cql.compiler.clifford_t,precision=c['synthesis_precision'])
                angle_counts=logical_counts(build)
                angle_cost[op.angle.hex()]=dict(T=t_count(angle_counts),clifford=sum(
                    angle_counts['actions'].get('qlx_standard_'+name,0) for name in ('h','s','sdg','x','y','z')))
        nt=sum(op.name in ('t','tdg') for op in ops)+sum(angle_cost[op.angle.hex()]['T'] for op in ops if op.name=='rz')
        assert nt==t_count(counts)
        routed,layout=nearest_neighbor_route(ops,width)
        for ports in c['source_ports']:
          for buffer_policy in ('bounded','continuous'):
            timing=schedule_patch_gates(routed,layout['width'],angle_cost,ports,c['source_interval_cycles'],c['code_distance'],buffer_policy)
            with gzip.open(OUT/f'batch{batch}_ports{ports}_{buffer_policy}_schedule.json.gz','wt') as f:json.dump(timing['events'],f)
            rows.append(dict(batch=batch,ports=ports,buffer_policy=buffer_policy,live_logical_qubits=width,packed_patches=layout['width'],
                clean_ancillas=width-2*c['L']**2,T_count=nt,rotation_classes=angle_cost,
                compiled_counts=counts,uncomputations=sum(op.name=='unand' for op in ops),
                authored_CX=sum(op.name=='cx' for op in ops),routed_CX=sum(g['name']=='cx' for g in routed),
                logical_clocks=timing['logical_clocks'],seconds=timing['seconds'],
                required_buffer_upper_estimate=timing['required_buffer_upper_estimate'],
                source_only_clocks=timing['source_only_clocks'],
                scope='Explicit P0 compilation; NN data-patch schedule and shared supply; magic-state transport abstracted, no physical surgery or decoder validation'))
        print('HWP batch',batch,'T',nt,'width',width,flush=True)
    record=dict(rows=rows,configuration=c,cudaq_logical_version=importlib.metadata.version('cudaq-logical'),
        source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('python/cudaq_hwp_stage.py','python/hwp_stage.py','configs/implementation_checks.json')})
    (OUT/'results.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':run()
