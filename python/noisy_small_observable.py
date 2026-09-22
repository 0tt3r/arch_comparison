"""Exact noisy eight-qubit observables, including particle-number postselection.

This bounded density-matrix calculation uses the saved native gate list and
independent primitive Pauli channels. It does not extrapolate to L >= 4.
"""

from project_paths import ROOT, project_file
import gzip
import json
from pathlib import Path
import time

import numpy as np

from hubbard_circuits import Operation
from nisq_compilation import schedule
from nisq_observables import idle_probability
from small_dynamics import HubbardSector,rectangle


class DensityNoise:
    def __init__(self,n):
        self.n=n;self.indices=np.arange(2**n)
        self.flip=[self.indices^(1<<(n-1-q)) for q in range(n)]
        self.bits=[(self.indices>>(n-1-q))&1 for q in range(n)]
        self.sign=[1-2*b for b in self.bits]
        self.same=[b[:,None]==b[None,:] for b in self.bits]

    def depolarize(self,rho,qubits,p):
        strength=4**len(qubits)*p/(4**len(qubits)-1)
        twirled=rho
        for q in qubits:
            ix=self.flip[q]
            twirled=(twirled+twirled[ix[:,None],ix[None,:]])*self.same[q]/2
        return (1-strength)*rho+strength*twirled

    def gate(self,rho,gate):
        q=gate.qubits[0];ix=self.flip[q]
        if gate.name=='x':return rho[ix[:,None],ix[None,:]]
        if gate.name=='h':
            left=(self.sign[q][:,None]*rho+rho[ix,:])/np.sqrt(2)
            return (left*self.sign[q][None,:]+left[:,ix])/np.sqrt(2)
        if gate.name=='rx':
            c,s=np.cos(gate.angle/2),np.sin(gate.angle/2)
            left=c*rho-1j*s*rho[ix,:]
            return c*left+1j*s*left[:,ix]
        if gate.name in ('rz','rzz'):
            values=self.sign[q]
            if gate.name=='rzz':values=values*self.sign[gate.qubits[1]]
            phase=np.exp(-.5j*gate.angle*values)
            return phase[:,None]*rho*phase.conj()[None,:]
        raise ValueError(gate.name)


def evaluate(record,gates):
    n=record['schedule']['width']
    if n>8:raise ValueError('Density calibration is deliberately limited to eight qubits')
    noise=DensityNoise(n);c=record['configuration'];p=c['p']
    rho=np.zeros((2**n,2**n),dtype=complex);rho[0,0]=1
    for q in range(n):
        ix=noise.flip[q];rho=(1-p)*rho+p*rho[ix[:,None],ix[None,:]]
    sched,events=schedule(gates,n,gate_ns=c['gate_seconds']*1e9,
        reset_ns=c['reset_seconds']*1e9,measurement_ns=c['measurement_seconds']*1e9)
    layers={}
    for start,gate in events:layers.setdefault(start,[]).append(gate)
    for layer in range(sched['depth_layers']):
        busy=set()
        for gate in layers.get(layer,[]):
            rho=noise.gate(rho,gate)
            rho=noise.depolarize(rho,gate.qubits,p if len(gate.qubits)==2 else p/10)
            busy.update(gate.qubits)
        for q in set(range(n))-busy:rho=noise.depolarize(rho,(q,),idle_probability(c))
    probabilities=rho.diagonal().real.copy()
    for q in range(n):probabilities=(1-p)*probabilities+p*probabilities[noise.flip[q]]
    mapping={int(k):v for k,v in record['mode_to_wire'].items()}
    up=sum(noise.bits[wire] for mode,wire in mapping.items() if mode%2==0)
    down=sum(noise.bits[wire] for mode,wire in mapping.items() if mode%2==1)
    D=np.prod([noise.bits[q] for q in record['observable_wires']],axis=0)
    selected={'unselected':np.ones(2**n,dtype=bool),'spin_parity':(up%2==0)&(down%2==0),
              'spin_particle_number':(up==2)&(down==2)}
    result={}
    for name,mask in selected.items():
        acceptance=float(probabilities[mask].sum())
        result[name]={'acceptance':acceptance,'doublon':float(np.dot(probabilities*mask,D)/acceptance)}
    assert np.isclose(result['spin_parity']['acceptance'],record['postselection']['acceptance'],atol=2e-11)
    result['trace_error']=float(abs(np.trace(rho)-1))
    result['hermiticity_error']=float(np.max(abs(rho-rho.conj().T)))
    result['min_eigenvalue']=float(np.linalg.eigvalsh(rho).min())
    return result


def main():
    root=(ROOT / 'results/refinement_v3');rows=[];model=HubbardSector(rectangle(2,2))
    exact=model.doublon(model.exact(1),(1,1))
    for r in (4,16,64):
        path=root/('parallel_L2_observable'+('' if r==4 else f'_r{r}'))
        record=json.loads((path/'result.json').read_text())
        with gzip.open(path/'native_gates.json.gz','rt') as f:
            gates=[Operation(**{**g,'qubits':tuple(g['qubits'])}) for g in json.load(f)]
        t=time.monotonic();row=evaluate(record,gates)
        ideal=model.doublon(model.trotter(1,r),(1,1))
        for name in ('unselected','spin_parity','spin_particle_number'):
            row[name]['bias_from_noiseless_circuit']=row[name]['doublon']-ideal
            row[name]['bias_from_exact_dynamics']=row[name]['doublon']-exact
        row.update({'L':2,'r':r,'tau':1.,'ideal_circuit_D':ideal,'exact_dynamics_D':exact,
                    'elapsed_seconds':time.monotonic()-t,'scope':'exact noisy eight-qubit calculation only'})
        rows.append(row);print('Noisy r',r,json.dumps(row),flush=True)
        tmp=root/'noisy_small.json.tmp';tmp.write_text(json.dumps(rows,indent=2)+'\n');tmp.replace(root/'noisy_small.json')


if __name__=='__main__':main()
