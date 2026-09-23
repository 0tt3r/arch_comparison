"""Circuit-level Clifford probe of a specified small rotated-patch growth.

Ideal small-patch preparation isolates growth faults. Fresh qubits to the right
are Z prepared and those below are X prepared; the target patch is stabilized
for d_f rounds. X/Z memory tests constrain Pauli components, not an arbitrary
non-Clifford injection by themselves. Stim supplies NN extraction schedules.
"""
import argparse
import hashlib
import json

import numpy as np
import pymatching
import stim
from scipy.stats import beta

from project_paths import ROOT

OUT=ROOT/'results/implementation_v1/star_growth'


def template(d,basis):
    circuit=stim.Circuit.generated('surface_code:rotated_memory_'+basis,distance=d,rounds=2)
    coords=circuit.get_final_qubit_coordinates()
    data=[q for q,(x,y) in coords.items() if x%2==1 and y%2==1]
    body=[];begun=False
    for instruction in circuit.flattened():
        if instruction.name=='H':begun=True
        if begun and instruction.name not in ('DETECTOR','SHIFT_COORDS','OBSERVABLE_INCLUDE'):
            body.append(instruction)
            if instruction.name=='MR':break
    xanc=set(q.value for ins in body if ins.name=='H' for q in ins.targets_copy())
    support={a:[] for a in xanc}
    for ins in body:
        if ins.name=='CX':
            t=ins.targets_copy()
            for a,b in zip(t[::2],t[1::2]):
                if a.value in xanc:support[a.value].append(b.value)
    return coords,data,body,support


def growth_circuit(di,df,p,basis='x'):
    coordinates,data,body,xsupport=template(df,basis)
    inverse={tuple(xy):q for q,xy in coordinates.items()}
    old_coords,old_data,old_body,_=template(di,basis)
    mapping={q:inverse[tuple(xy)] for q,xy in old_coords.items()}
    old=set(mapping[q] for q in old_data)
    anc=sorted(set(coordinates)-set(data));circuit=stim.Circuit()
    def add(name,targets,noise=False):
        circuit.append(name,targets)
        if noise and p:
            if name=='CX':circuit.append('DEPOLARIZE2',targets,p)
            elif name=='H':circuit.append('DEPOLARIZE1',targets,p/10)
            elif name=='RX':circuit.append('Z_ERROR',targets,p)
            elif name in ('R','MR'):circuit.append('X_ERROR',targets,p)
    add('RX' if basis=='x' else 'R',sorted(old))
    add('R',[mapping[q] for q in old_coords if q not in old_data])
    for _ in range(2):
        for ins in old_body:
            add(ins.name,[mapping[q.value] for q in ins.targets_copy()])
    right=[q for q in data if q not in old and coordinates[q][1]<2*di]
    bottom=[q for q in data if q not in old and coordinates[q][1]>=2*di]
    if right:add('R',right,True)
    if bottom:add('RX',bottom,True)
    add('R',anc,True)
    last_record={};measurements=circuit.num_measurements
    previous_record=None
    for round_index in range(df):
        busy=set()
        for ins in body:
            targets=[q.value for q in ins.targets_copy()]
            if ins.name=='TICK':
                if p:
                    idle=sorted(set(coordinates)-busy)
                    if idle:circuit.append('DEPOLARIZE1',idle,p/10)
                busy=set();circuit.append('TICK');continue
            if ins.name=='MR':
                if p:
                    circuit.append('X_ERROR',targets,p)
                    circuit.append('DEPOLARIZE1',data,p) # 100 ns measurement idle
                last_record={q:measurements+j for j,q in enumerate(targets)}
                measurements+=len(targets)
            add(ins.name,targets,True);busy.update(targets)
        # 4 CX + 2 H layers (60 ns), 100 ns measurement, and an
        # 80 ns reset/feedforward allowance make the declared 240 ns cycle.
        if p:circuit.append('DEPOLARIZE1',sorted(coordinates),.8*p)
        if round_index==0:
            # Exact flows only at the transition. Later extraction rounds use
            # local time differences, avoiding a dense first-round reference.
            relations=[]
            for flow in circuit.without_noise().flow_generators():
                if flow.input_copy().weight==0 and flow.output_copy().weight==0:
                    indices=flow.measurements_copy()
                    if indices:relations.append(sum(1<<i for i in indices))
            # Row-reduce on the new checks first. The generic flow basis can
            # mix several new checks, creating artificial hyperedges. Source-
            # only constraints have zero noise here and need not be detectors.
            pivot=0
            for column in sorted(last_record.values()):
                found=next((i for i in range(pivot,len(relations)) if (relations[i]>>column)&1),None)
                if found is None:continue
                relations[pivot],relations[found]=relations[found],relations[pivot]
                for i in range(len(relations)):
                    if i!=pivot and (relations[i]>>column)&1:relations[i]^=relations[pivot]
                pivot+=1
            for bits in relations[:pivot]:
                circuit.append('DETECTOR',[stim.target_rec(i-measurements) for i in range(measurements) if (bits>>i)&1])
        else:
            for a in anc:
                circuit.append('DETECTOR',[stim.target_rec(last_record[a]-measurements),stim.target_rec(previous_record[a]-measurements)])
        previous_record=last_record
    if p:circuit.append('Z_ERROR' if basis=='x' else 'X_ERROR',data,p)
    add('MX' if basis=='x' else 'M',data)
    final={q:measurements+j for j,q in enumerate(data)};total=circuit.num_measurements
    # Readout closes checks of the measured basis, retaining logical parity.
    if basis=='x':supports=xsupport
    else:
        zanc=set(anc)-set(xsupport);supports={a:[] for a in zanc}
        for ins in body:
            if ins.name=='CX':
                t=ins.targets_copy()
                for a,b in zip(t[::2],t[1::2]):
                    if b.value in zanc:supports[b.value].append(a.value)
    for a,support in supports.items():
        circuit.append('DETECTOR',[stim.target_rec(last_record[a]-total)]+[stim.target_rec(final[q]-total) for q in support])
    logical=[q for q in data if coordinates[q][0 if basis=='x' else 1]==1]
    circuit.append('OBSERVABLE_INCLUDE',[stim.target_rec(final[q]-total) for q in logical],0)
    return circuit


def run():
    OUT.mkdir(parents=True,exist_ok=True);rows=[]
    for di,df in ((3,3),(3,5),(3,7),(5,5),(5,7)):
        for p in (.001,.0001):
            for basis in ('x','z'):
                circuit=growth_circuit(di,df,p,basis)
                ideal=growth_circuit(di,df,0,basis)
                detections,observables=ideal.compile_detector_sampler(seed=8127).sample(128,separate_observables=True)
                assert not detections.any() and not observables.any()
                dem=circuit.detector_error_model(decompose_errors=True)
                decoder=pymatching.Matching.from_detector_error_model(dem)
                seed=20260923+di*10000+df*100+int(p*1e6)+(basis=='z')
                sampler=circuit.compile_detector_sampler(seed=seed)
                shots=failures=0
                # Fixed sample size keeps the binomial confidence interval valid.
                while shots<500000:
                    defects,obs=sampler.sample(25000,separate_observables=True)
                    failures+=int(np.any(decoder.decode_batch(defects)!=obs,axis=1).sum());shots+=len(obs)
                lo=0 if failures==0 else float(beta.ppf(.025,failures,shots-failures+1))
                hi=1 if failures==shots else float(beta.ppf(.975,failures+1,shots-failures))
                name=f'd{di}_to_d{df}_{basis}_p{p:g}.stim';(OUT/name).write_text(str(circuit)+'\n')
                rows.append(dict(di=di,df=df,p=p,basis=basis,shots=shots,failures=failures,rate=failures/shots,
                    confidence_interval=[lo,hi],seed=seed,circuit=name,
                    detectors=circuit.num_detectors,physical_qubits=2*df*df-1,
                    scope='Ideal input encoding; noisy direct growth/stabilization and final readout; component memory failure with MWPM'))
                (OUT/'results.json').write_text(json.dumps(dict(rows=rows,stim=stim.__version__,pymatching=pymatching.__version__,
                    source_sha256=hashlib.sha256((ROOT/'python/star_growth_probe.py').read_bytes()).hexdigest()),indent=2)+'\n')
                print(di,df,p,basis,failures,shots,flush=True)


def check():
    saved=json.loads((OUT/'results.json').read_text())
    assert saved['source_sha256']==hashlib.sha256((ROOT/'python/star_growth_probe.py').read_bytes()).hexdigest()
    assert len(saved['rows'])==20
    for row in saved['rows']:
        di,df,p,basis=[row[k] for k in ('di','df','p','basis')]
        circuit=growth_circuit(di,df,p,basis)
        assert (OUT/row['circuit']).read_text()==str(circuit)+'\n'
        pymatching.Matching.from_detector_error_model(circuit.detector_error_model(decompose_errors=True))
        ideal=growth_circuit(di,df,0,basis)
        det,obs=ideal.compile_detector_sampler(seed=1).sample(256,separate_observables=True)
        assert not det.any() and not obs.any()
        # An inserted logical operator must change the output and evade all
        # stabilizer checks. This guards against vacuous observable recording.
        coords,data,_,_=template(df,basis)
        logical=[q for q in data if coords[q][1 if basis=='x' else 0]==1]
        control=stim.Circuit();inserted=False
        for ins in ideal:
            if ins.name==('MX' if basis=='x' else 'M'):
                control.append('Z_ERROR' if basis=='x' else 'X_ERROR',logical,1);inserted=True
            control.append(ins)
        det,obs=control.compile_detector_sampler(seed=2).sample(256,separate_observables=True)
        assert inserted and not det.any() and obs.all()
        shots,failures=row['shots'],row['failures']
        assert shots==500000 and row['rate']==failures/shots
        interval=[0 if failures==0 else float(beta.ppf(.025,failures,shots-failures+1)),
                  1 if failures==shots else float(beta.ppf(.975,failures+1,shots-failures))]
        assert interval==row['confidence_interval']
    print('PASS: 20 growth circuits, graph decomposition, noiseless and logical-flip controls, fixed-sample intervals')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    check() if args.check else run()
