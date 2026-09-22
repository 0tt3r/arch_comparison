"""Remove avoidable idle noise before product-state preparation and after readout.

The ASAP body and attempt duration are unchanged. Initialization may finish just
before a wire's first use; each wire may be measured after its last use. Initial
Neel X gates can move to immediately before their next gate. Independent idle
depolarization commutes with those X gates. This is an exact noise-accounting
change for the declared model, not a new dynamics approximation.
"""

from project_paths import ROOT, project_file
import gzip
import json
from math import exp, expm1, log, log1p
from pathlib import Path

from hubbard_circuits import Operation
from model import finish, log_gamma2_depolarizing, workload
from nisq_compilation import schedule
from nisq_observables import idle_probability


def boundary_gaps(gates,sched,events):
    n=sched['width'];depth=sched['depth_layers']
    per_wire=[[] for _ in range(n)]
    for start,gate in events:
        for q in gate.qubits:per_wire[q].append((start,gate))
    head=[];tail=[]
    for entries in per_wire:
        entries.sort(key=lambda pair:pair[0])
        if not entries:
            head.append(depth);tail.append(0);continue
        first,gate=entries[0]
        if gate.name=='x' and len(entries)>1:
            head.append(max(0,entries[1][0]-1))
        else:head.append(first)
        tail.append(depth-entries[-1][0]-1)
    assert sum(head)+sum(tail)<=sched['idle_qubit_layers']
    return head,tail


def refine(record,gates):
    c=record['configuration']
    sched,events=schedule(gates,record['schedule']['width'],gate_ns=c['gate_seconds']*1e9,
                          reset_ns=c['reset_seconds']*1e9,measurement_ns=c['measurement_seconds']*1e9)
    assert json.loads(json.dumps(sched))==json.loads(json.dumps(record['schedule']))
    head,tail=boundary_gaps(gates,sched,events)
    exposure=record['observable_exposure']
    removed_cone=sum(head[q] for q in exposure['input_support'])+sum(tail[q] for q in record['observable_wires'])
    removed_all=sum(head)+sum(tail)
    assert removed_cone<=exposure['locations']['idle']
    idle=idle_probability(c)
    chi=exposure['chi']-removed_cone*log_gamma2_depolarizing(idle,1)
    mapping={int(k):v for k,v in record['mode_to_wire'].items()}
    initial_modes=sorted(mapping)
    post=record['postselection'];logs=[]
    for index,mask in enumerate((1,2,3)):
        removed=(sum(head[q] for q,mode in enumerate(initial_modes) if mask&(1<<(mode%2)))
                 +sum(tail[wire] for mode,wire in mapping.items() if mask&(1<<(mode%2))))
        logs.append(post['log_character_visibilities'][index]-removed*log1p(-4*idle/3))
    acceptance=(1+sum(exp(value) for value in logs))/4
    no_fault_log=post['log_no_fault_probability']-removed_all*log1p(-idle)
    w=workload(record['L'],record['tau'],record['steps'],record['algorithm'],c['u_over_j'])
    return {'policy':'prepare at first use; measure at last use; no readout crosstalk assumed',
            'removed_idle_qubit_layers':removed_all,'removed_cone_idle_layers':removed_cone,
            'head_idle_layers':head,'tail_idle_layers':tail,
            'chi':chi,'postselection_acceptance':acceptance,
            'postselection_bias_upper_bound':-expm1(min(0,no_fault_log-log(acceptance))),
            'cost':finish(w,'physical_local_PEC',chi,sched['physical_copy_qubits'],sched['seconds_per_attempt'],c)}


def main():
    root=ROOT/'results/nisq_v2'
    rows=[]
    for path in sorted(root.glob('*_L*_r*_*/result.json')):
        record=json.loads(path.read_text())
        with gzip.open(path.parent/'native_gates.json.gz','rt') as f:
            gates=[Operation(**{**g,'qubits':tuple(g['qubits'])}) for g in json.load(f)]
        rows.append({'case':path.parent.name,**refine(record,gates)})
    (root/'boundary_timing.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':main()
