"""Restartable U=4 center-hole sector calculation, separately pinned from Neel runs."""
import argparse
import json
import os
from pathlib import Path
import socket
import subprocess
import time

import numpy as np

from large_dynamics import coefficients
from project_paths import ROOT
from small_dynamics import rectangle
from willow_archive import sha, save
from willow_benchmark import load_reduced
from willow_reference import hopping_order, OUT

SOURCE=ROOT/'cpp/willow_sector.cpp'
BINARY=ROOT/'build/willow_sector'


def sources():
    names=['cpp/willow_sector.cpp','cpp/hubbard_sector.cpp','python/willow_sector.py',
           'python/willow_reference.py','python/large_dynamics.py','python/small_dynamics.py',
           'results/willow_v1/reduced.json.gz','configs/willow_reference.json']
    return {name:sha(ROOT/name) for name in names}


def compile_solver(binary=BINARY):
    binary=Path(binary);binary.parent.mkdir(parents=True,exist_ok=True)
    command=['g++','-O3','-std=c++17','-fopenmp',str(SOURCE),'-o',str(binary)]
    subprocess.run(command,check=True)
    return dict(command=command,compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
                binary_sha256=sha(binary),sources=sources())


def masks(L, occupations=None):
    sites=rectangle(L,L)
    if occupations is None:
        hole=(L//2-1,L//2-1)
        occupations=[dict(site=list(xy),spin=s,n=int(xy!=hole and ((sum(xy)%2==0)==(s=='up'))))
                     for xy in sites for s in ('up','down')]
    return tuple(sum(1<<sites.index(tuple(v['site'])) for v in occupations if v['spin']==s and v['n'])
                 for s in ('up','down'))


def solve(L,tau,mode,r,out,threads=16,tolerance=1e-13,occupations=None,binary=BINARY,state_path=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    iu,idown=masks(L,occupations)
    prefix=f'L{L}_{mode}_tau{tau:g}_r{r}'
    contract=dict(L=L,U=4,flux=0,initial_up=iu,initial_down=idown,tau=tau,mode=mode,r=r,
                  threads=threads,tolerance=tolerance,binary_sha256=sha(binary),sources=sources())
    targets=[out/f'{prefix}_k{k}.json' for k in range(1,(r if mode=='trajectory' else 1)+1)]
    existing={}
    for path in targets:
        if path.exists():
            existing[path]=json.loads(path.read_text())
            assert existing[path]['contract']==contract, f'Stale contract: {path}'
            assert existing[path]['status']=='complete'
    if len(existing)==len(targets) and state_path is None:return list(existing.values())
    coeff,tail=coefficients(4*L*(L-1)+L*L,tau,tolerance)
    cp=out/f'{prefix}.coefficients.txt';np.savetxt(cp,np.column_stack((coeff.real,coeff.imag)),fmt='%.17g')
    op=out/f'{prefix}.order.txt';np.savetxt(op,hopping_order(L,L),fmt='%d')
    command=[str(binary),str(L),str(L),str(iu),str(idown),str(tau),mode,str(r),str(cp),str(op)]
    if state_path is not None:command.append(str(state_path))
    env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_PROC_BIND':'spread','OMP_PLACES':'cores'}
    start=time.monotonic();records=[]
    # One producer owns this series. Completed prefix values remain immutable;
    # after interruption the state evolution restarts, avoiding large checkpoints.
    with (out/f'{prefix}.stderr.log').open('a') as err:
        process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=err,text=True,env=env)
        try:
            for k,line in enumerate(process.stdout,1):
                result=json.loads(line)
                assert k<=len(targets)
                assert abs(result['norm_squared']-1)<1e-9
                assert abs(result['tau']-(k*tau if mode=='trajectory' else tau))<1e-12
                record=dict(contract=contract,sample=k,result=result,status='complete',
                    state_tail_bound=(k if mode=='trajectory' else 1)*tail if mode!='circuit' else None,
                    chebyshev_degree=len(coeff)-1,coefficient_sha256=sha(cp),order_sha256=sha(op),
                    host=socket.gethostname(),pbs_jobid=os.getenv('PBS_JOBID'),command=command,
                    elapsed_seconds=time.monotonic()-start,
                    scope='Finite OBC U=4 zero-flux product-state sector. Centered interaction differs from physical H only by a scalar. Numerical tail excludes rounding.')
                path=targets[k-1]
                if path in existing:
                    for name in ('doublons','z_up','z_down'):
                        np.testing.assert_allclose(existing[path]['result'][name],result[name],atol=1e-10,rtol=0)
                    record=existing[path]
                else:save(path,record)
                records.append(record)
                print(prefix,k,'seconds',round(time.monotonic()-start,2),'norm',result['norm_squared'],flush=True)
            assert process.wait()==0,'Solver failed; inspect stderr'
            assert len(records)==len(targets),'Incomplete solver stream'
        finally:
            if process.poll() is None:process.terminate();process.wait()
    return records


def run(task):
    build=json.loads((OUT/'build.json').read_text())
    assert build['sources']==sources() and build['binary_sha256']==sha(BINARY)
    c=next(c for c in load_reduced()['cases'] if (c['L'],c['U'],c['flux'],c['steps'])==(4,4,'0',3))
    assert masks(4,c['initial_occupations'])==masks(4)
    cfg=json.loads((ROOT/'configs/willow_reference.json').read_text())
    threads=cfg['threads'];occupations=c['initial_occupations']
    if not os.getenv('PBS_JOBID'):raise RuntimeError('4x4 calculations require a PBS allocation')
    if task=='smoke':
        solve(4,.05,'exact',0,OUT/'smoke',threads,occupations=occupations)
        solve(4,.05,'circuit',3,OUT/'smoke',threads,occupations=occupations)
    elif task=='production':
        solve(4,cfg['time_spacing'],'trajectory',cfg['time_samples'],OUT/'continuous',threads,occupations=occupations)
        for tau in cfg['times']:
            for r in cfg['steps']:solve(4,tau,'circuit',r,OUT/'circuit',threads,occupations=occupations)
        solve(4,1.,'exact',0,OUT/'direct_check',threads,tolerance=1e-14,occupations=occupations)
    else:raise ValueError(task)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepare',action='store_true');p.add_argument('--task',choices=['smoke','production'])
    args=p.parse_args()
    if args.prepare:save(OUT/'build.json',compile_solver())
    if args.task:run(args.task)
