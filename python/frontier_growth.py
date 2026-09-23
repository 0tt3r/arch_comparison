"""Fixed-shot, restartable growth components toward the data-code distance."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import resource
import socket
import time

import numpy as np
import pymatching
from scipy.stats import beta
import stim

from project_paths import ROOT
from star_growth_probe import growth_circuit

OUT=ROOT/'results/frontier_v1/growth'


def cases():
    c=json.loads((ROOT/'configs/frontier_study.json').read_text())['growth']
    rows=[]
    for di,df in c['pairs']:
        for p in c['p_values']:
            for basis in c['bases']:
                key=f'd{di}_d{df}_{basis}_p{p:g}'
                seed=c['seed']+int(hashlib.sha256(key.encode()).hexdigest()[:8],16)
                rows.append(dict(key=key,di=di,df=df,p=p,basis=basis,shots=c['shots_per_case'],seed=seed))
    return rows


def sources():
    names=['python/frontier_growth.py','python/star_growth_probe.py','configs/frontier_study.json']
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}


def atomic(path,value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def execute(job):
    row,smoke=job;start=time.monotonic();folder=OUT/('smoke' if smoke else 'cases');folder.mkdir(parents=True,exist_ok=True)
    contract={**row,'shots':1000 if smoke else row['shots'],'sources':sources()};path=folder/(row['key']+'.json')
    if path.exists():
        previous=json.loads(path.read_text());assert previous['contract']==contract and previous['status']=='complete'
        return previous
    lock=path.with_suffix('.lock')
    with lock.open('x') as f:f.write(json.dumps(dict(pid=os.getpid(),host=socket.gethostname(),job=os.environ.get('PBS_JOBID'))))
    try:
        circuit=growth_circuit(row['di'],row['df'],row['p'],row['basis'])
        ideal=growth_circuit(row['di'],row['df'],0,row['basis'])
        det,obs=ideal.compile_detector_sampler(seed=7).sample(64,separate_observables=True)
        assert not det.any() and not obs.any()
        model=circuit.detector_error_model(decompose_errors=True)
        decoder=pymatching.Matching.from_detector_error_model(model)
        sampler=circuit.compile_detector_sampler(seed=row['seed']);shots=failures=0
        batch=max(128,min(5000,20000000//circuit.num_detectors))
        while shots<contract['shots']:
            n=min(batch,contract['shots']-shots);det,obs=sampler.sample(n,separate_observables=True)
            failures+=int(np.any(decoder.decode_batch(det)!=obs,axis=1).sum());shots+=n
        delta=.05;cells=len(cases())
        interval=[0 if failures==0 else float(beta.ppf(delta/2,failures,shots-failures+1)),
                  float(beta.ppf(1-delta/2,failures+1,shots-failures))]
        upper=float(beta.ppf(1-delta/cells,failures+1,shots-failures))
        text=str(circuit)+'\n';cpath=folder/(row['key']+'.stim');cpath.write_text(text)
        result=dict(status='complete',contract=contract,shots=shots,failures=failures,rate=failures/shots,
            interval95=interval,simultaneous_upper95=upper,circuit_sha256=hashlib.sha256(text.encode()).hexdigest(),
            detectors=circuit.num_detectors,physical_qubits=2*row['df']**2-1,
            elapsed_seconds=time.monotonic()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            host=socket.gethostname(),pbs_job=os.environ.get('PBS_JOBID'),stim=stim.__version__,pymatching=pymatching.__version__,
            scope='Ideal encoded input; decoded growth + target memory + readout. Not a full injected-state channel.')
        atomic(path,result);print(row['key'],failures,shots,result['elapsed_seconds'],flush=True);return result
    finally:lock.unlink()


def check(smoke=False):
    rows=cases()
    if smoke:rows=[v for v in rows if v['di']==9 and v['df']==25 and v['p']==.001]
    folder=OUT/('smoke' if smoke else 'cases');records=[]
    for row in rows:
        path=folder/(row['key']+'.json');value=json.loads(path.read_text());contract=value['contract']
        assert value['status']=='complete' and contract=={**row,'shots':1000 if smoke else row['shots'],'sources':sources()}
        assert value['shots']==contract['shots'] and not path.with_suffix('.lock').exists()
        assert hashlib.sha256(path.with_suffix('.stim').read_bytes()).hexdigest()==value['circuit_sha256']
        records.append(value)
    result=dict(status='PASS: complete fixed-shot census and source/circuit contracts',cases=len(rows),shots=sum(v['shots'] for v in records))
    atomic(folder/'validation.json',result);print(result);return result


def run(smoke=False,workers=8):
    if not smoke and not os.environ.get('PBS_JOBID'):raise RuntimeError('Run the larger fixed-shot campaign in PBS')
    rows=cases()
    if smoke:rows=[v for v in rows if v['di']==9 and v['df']==25 and v['p']==.001]
    with ProcessPoolExecutor(max_workers=min(workers,len(rows))) as pool:list(pool.map(execute,[(v,smoke) for v in rows]))
    check(smoke)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--smoke',action='store_true');parser.add_argument('--check',action='store_true');parser.add_argument('--workers',type=int,default=8)
    args=parser.parse_args();check(args.smoke) if args.check else run(args.smoke,args.workers)
