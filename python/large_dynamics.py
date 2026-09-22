"""Bounded, restartable matrix-free calibration; large cases belong on PBS."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import socket
import resource
import subprocess
import time

import numpy as np
from scipy.special import gammaln, jv

from project_paths import ROOT

SOURCE = ROOT/'cpp/hubbard_sector.cpp'
BINARY = ROOT/'build/hubbard_sector'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def coefficients(scale, tau, tolerance=1e-12):
    """Chebyshev coefficients and a conservative analytic Bessel-tail bound.

    |J_k(a)| <= (a/2)^k/k! exp(a²/[4(k+1)]). Subsequent bounds have
    ratio at most a/[2(k+1)], so a geometric tail bounds the state error.
    Floating-point error is checked separately against independent small cases.
    """
    a=abs(scale*tau)
    if a==0:return np.array([1.,0.],complex),0.
    degree=max(1,math.ceil(a))
    while True:
        k=degree+1;ratio=a/(2*(k+1))
        tail=2*math.exp(k*math.log(a/2)-gammaln(k+1)+a*a/(4*(k+1)))/(1-ratio)
        if tail<tolerance:break
        degree+=1
    k=np.arange(degree+1)
    c=2*(-1j)**k*jv(k,scale*tau);c[0]/=2
    return c,float(tail)


def compile_solver():
    BINARY.parent.mkdir(exist_ok=True)
    command=['g++','-O3','-std=c++17','-fopenmp',str(SOURCE),'-o',str(BINARY)]
    subprocess.run(command,check=True)
    return {'command':command,'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
            'source_sha256':sha(SOURCE),'binary_sha256':sha(BINARY)}


def atomic_json(path,data):
    temporary=path.with_suffix('.json.tmp');temporary.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def solve(nx,ny,tau,algorithm,r,out,threads=1,tolerance=1e-12,state_path=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    name=f'{nx}x{ny}_tau{tau:g}_{algorithm}_r{r}_t{threads}'
    source_hash=sha(SOURCE);binary_hash=sha(BINARY)
    contract={'nx':nx,'ny':ny,'tau':tau,'algorithm':algorithm,'r':r,'threads':threads,
              'tolerance':tolerance,'source_sha256':source_hash,'binary_sha256':binary_hash,
              'driver_sha256':sha(__file__)}
    path=out/f'{name}.json'
    if path.exists() and state_path is None:
        saved=json.loads(path.read_text())
        if saved['contract']!=contract:raise ValueError(f'Stale result contract: {path}')
        return saved
    edges=(nx-1)*ny+(ny-1)*nx;scale=2*edges+nx*ny
    coeff,tail=coefficients(scale,tau,tolerance)
    coeff_path=out/f'{name}.coefficients.txt'
    np.savetxt(coeff_path,np.column_stack([coeff.real,coeff.imag]),fmt='%.17g')
    command=[str(BINARY),str(nx),str(ny),str(tau),algorithm,str(r),str(coeff_path)]
    if state_path is not None:command.append(str(state_path))
    start=time.time();usage_before=resource.getrusage(resource.RUSAGE_CHILDREN)
    proc=subprocess.run(command,check=True,text=True,capture_output=True,
                        env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_PROC_BIND':'spread','OMP_PLACES':'cores'})
    result=json.loads(proc.stdout)
    usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    record={'contract':contract,'result':result,'chebyshev_degree':len(coeff)-1,
            'chebyshev_state_tail_bound':tail if algorithm=='exact' else None,
            'coefficient_sha256':sha(coeff_path),'host':socket.gethostname(),'pbs_jobid':os.getenv('PBS_JOBID'),
            'elapsed_seconds':time.time()-start,'command':command,'status':'complete',
            'child_cpu_seconds':usage.ru_utime+usage.ru_stime-usage_before.ru_utime-usage_before.ru_stime,
            'children_peak_rss_kib':usage.ru_maxrss,
            'scope':'Finite OBC patch, Neel state, U/J=4; site order is row snake; no bulk error bound'}
    if algorithm!='benchmark' and abs(result['norm_squared']-1)>1e-8:raise ValueError('Norm check failed')
    atomic_json(path,record)
    print(name, 'seconds',round(record['elapsed_seconds'],2),'norm',result['norm_squared'],flush=True)
    return record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--compile',action='store_true')
    p.add_argument('--patch',default='4x4');p.add_argument('--tau',type=float,default=1.)
    p.add_argument('--algorithms',nargs='+',default=['exact','plaquette','termwise'])
    p.add_argument('--steps',nargs='+',type=int,default=[16,32,64]);p.add_argument('--threads',type=int,default=16)
    p.add_argument('--out',type=Path,default=ROOT/'results/dynamics_4x4')
    args=p.parse_args()
    if args.compile:print(json.dumps(compile_solver(),indent=2),flush=True)
    nx,ny=map(int,args.patch.split('x'))
    if nx*ny>=16 and not os.getenv('PBS_JOBID'):
        raise RuntimeError('Run the full 4x4 allocation on a PBS compute node')
    for algorithm in args.algorithms:
        for r in ([0] if algorithm=='exact' else args.steps):
            solve(nx,ny,args.tau,algorithm,r,args.out,args.threads)


if __name__=='__main__':main()
