"""STAR logical RUS reference: Pauli covariance and parallel retry timing.

This is an explicit logical gadget model. It does not certify the residual noise
or nearest-neighbor surgery layout of a surface-code implementation.
"""

from project_paths import ROOT, project_file
import json
from math import ceil, exp, log, log2, log1p
from pathlib import Path

import numpy as np

from hubbard_circuits import build_workload
from model import accuracy, finish, minimum_distance, star_exposure, workload

I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],dtype=complex)
Y=np.array([[0,-1j],[1j,0]],dtype=complex);Z=np.diag([1.,-1.])
PAULIS=(I,X,Y,Z)
CX=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],dtype=complex)


def injection_kraus(theta,outcome,before=(0,0),after=(0,0)):
    """Ancilla exp(-i theta Z)|+>, CX(data,ancilla), then ancilla Z readout."""
    anc=np.array([np.exp(-1j*theta),np.exp(1j*theta)])/np.sqrt(2)
    embed=np.kron(I,anc[:,None])
    matrix=np.kron(PAULIS[after[0]],PAULIS[after[1]])@CX@np.kron(PAULIS[before[0]],PAULIS[before[1]])@embed
    return matrix[[outcome,outcome+2],:]


def expected_max_attempts(batch,cap=64):
    """Mean max of batch independent geometric(1/2) attempts, each capped."""
    if batch<1:return 0.
    return sum(1-(1-2.**(-k))**batch for k in range(cap))


def rotation_batches(operations,width,lanes):
    """Disjoint full-JW-support batches preserving every per-patch dependency.

    All operations in a layer finish before the next layer. This explicit
    barrier policy is conservative relative to asynchronous dispatch.
    """
    ready=[0]*width;layers={}
    for op in operations:
        layer=max(ready[q] for q in op.qubits)
        for q in op.qubits:ready[q]=layer+1
        layers.setdefault(layer,[]).append(op)
    batches=[];clifford_layers=0
    for layer in sorted(layers):
        rotations=[op for op in layers[layer] if op.name=='rotate']
        if any(op.name!='rotate' for op in layers[layer]):clifford_layers+=1
        for i in range(0,len(rotations),lanes):batches.append(len(rotations[i:i+lanes]))
    return batches,clifford_layers


def parallel_star(L,c,lanes,logical_bias='worst_case'):
    ops,_=build_workload(L,1.,4)
    w=workload(L,1.,4,'plaquette',c['u_over_j']);eps,_=accuracy(c)
    chi=star_exposure(w,c,'star_floor');log_amp=chi/2 if logical_bias=='worst_case' else 0.
    if logical_bias not in ('worst_case','covariant_reference'):raise ValueError(logical_bias)
    count=w.rotations+w.fixed_T
    cap=ceil(log2(count/(eps*c['logical_fraction']/10)))
    log_acceptance=count*log1p(-2.**(-cap))
    batches,clifford=rotation_batches(ops,2*L*L,lanes)
    expected_rounds=sum(expected_max_attempts(n,cap) for n in batches)
    # Nine clocks per attempt gives the accepted 18-clock serial mean.
    clocks=2+clifford+9*expected_rounds
    tiles=ceil(1.5*(2*w.sites+6))+lanes*c['star_injection_tiles']
    def duration(d):
        return clocks*d*c['cycle_seconds']+expected_rounds*c['feedforward_seconds']+c['reset_seconds']+c['measurement_seconds']
    d,seconds,_=minimum_distance(tiles,duration,log(eps*c['logical_fraction'])-log_amp,c)
    row=finish(w,'star_parallel',chi-log_acceptance,2*tiles*d*d,seconds,c,distance=d,lanes=lanes,
        logical_bias=logical_bias,rotation_batches=len(batches),max_batch=max(batches),
        mean_retry_rounds=expected_rounds,retry_cap=cap,abort_probability=-np.expm1(log_acceptance),
        source_chi=chi,clifford_layers=clifford,
        implementation='full-JW-support logical batches; nine-clock attempt and six tiles per lane assumed')
    return row


def run():
    root=ROOT;c=json.loads((project_file('assumptions.json')).read_text());rows=[]
    for eps in (.001,.01):
        scenario={**c,'relative_error':eps/c['doublon_reference']}
        for L in (2,4,6,8):
            for lanes in (1,2,4,8):
                for bias in ('worst_case','covariant_reference'):
                    rows.append({**parallel_star(L,scenario,lanes,bias),'absolute_error':eps})
    out=root/'results/refinement_v3';out.mkdir(exist_ok=True,parents=True)
    (out/'star_parallel.json').write_text(json.dumps(rows,indent=2)+'\n')
    for eps in (.001,.01):
        for L in (2,4,6,8):
            selected=min((r for r in rows if r['L']==L and r['absolute_error']==eps and r['logical_bias']=='worst_case'),key=lambda r:r['log10_Q_min'])
            print('STAR',eps,L,'lanes',selected['lanes'],'Q',round(10**selected['log10_Q_min']),flush=True)


if __name__=='__main__':run()
