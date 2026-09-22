"""Bounded exact-sector Hubbard dynamics and local product-formula calibration.

Finite-reference truncation differences are measurements, not infinite-lattice
error certificates. All spatial patches use the same Neel pattern and U/J.
"""
import argparse
from itertools import combinations
import json
from pathlib import Path
import time

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


def rectangle(nx,ny):
    return [(x,y) for y in range(ny) for x in (range(nx) if y%2==0 else range(nx-1,-1,-1))]


def sector_masks(n,k):
    return np.array([sum(1<<i for i in occupied) for occupied in combinations(range(n),k)],dtype=np.int64)


def spin_hopping(masks,edges):
    lookup={int(mask):i for i,mask in enumerate(masks)}
    rows=[];cols=[];values=[]
    for col,mask in enumerate(masks):
        mask=int(mask)
        for i,j in edges:
            if ((mask>>i)^(mask>>j))&1:
                lo,hi=sorted((i,j));between=((1<<hi)-1)^((1<<(lo+1))-1)
                sign=(-1)**((mask&between).bit_count())
                rows.append(lookup[mask^(1<<i)^(1<<j)]);cols.append(col);values.append(-sign)
    return sparse.csr_matrix((values,(rows,cols)),shape=(len(masks),len(masks)),dtype=float)


class HubbardSector:
    def __init__(self,sites,u=4.,max_dimension=20000):
        self.sites=list(sites);self.n=len(sites);lookup={xy:i for i,xy in enumerate(sites)}
        initial_up=sum(1<<i for i,(x,y) in enumerate(sites) if (x+y)%2==0)
        initial_down=((1<<self.n)-1)^initial_up
        self.up=sector_masks(self.n,initial_up.bit_count())
        self.down=sector_masks(self.n,initial_down.bit_count())
        self.shape=(len(self.up),len(self.down));self.dimension=np.prod(self.shape).item()
        if self.dimension>max_dimension:raise ValueError(f'Dimension {self.dimension} exceeds bounded calibration cap')
        edges=[[],[]]
        for i,(x,y) in enumerate(sites):
            for xy,color in [((x+1,y),x%2),((x,y+1),y%2)]:
                if xy in lookup:edges[color].append((i,lookup[xy]))
        self.hopping=[]
        for group in edges:
            hu=spin_hopping(self.up,group);hd=spin_hopping(self.down,group)
            self.hopping.append((sparse.kron(hu,sparse.eye(len(self.down)),format='csr')+
                                 sparse.kron(sparse.eye(len(self.up)),hd,format='csr')))
        self.potential=np.zeros(self.shape)
        for i in range(self.n):
            self.potential+=u*np.outer(((self.up>>i)&1)-.5,((self.down>>i)&1)-.5)
        self.potential=self.potential.ravel()
        self.hamiltonian=self.hopping[0]+self.hopping[1]+sparse.diags(self.potential,format='csr')
        self.initial=np.zeros(self.dimension,dtype=complex)
        a=np.flatnonzero(self.up==initial_up)[0];b=np.flatnonzero(self.down==initial_down)[0]
        self.initial[a*self.shape[1]+b]=1

    def doublon(self,psi,site):
        i=self.sites.index(tuple(site))
        mask=np.outer((self.up>>i)&1,(self.down>>i)&1).ravel()
        return float(np.dot(np.abs(psi)**2,mask))

    def exact(self,tau):
        return expm_multiply((-1j*tau)*self.hamiltonian,self.initial,traceA=-1j*tau*self.potential.sum())

    def trotter(self,tau,steps):
        dt=tau/steps;psi=self.initial.copy()
        onsite=np.exp(-.5j*dt*self.potential)
        pink=-.5j*dt*self.hopping[0];gold=-1j*dt*self.hopping[1]
        for _ in range(steps):
            psi=onsite*psi
            psi=expm_multiply(pink,psi,traceA=0.)
            psi=expm_multiply(gold,psi,traceA=0.)
            psi=expm_multiply(pink,psi,traceA=0.)
            psi=onsite*psi
        return psi


def run(out):
    out.mkdir(parents=True,exist_ok=True);rows=[];spatial=[]
    # Use a shared observed coordinate and nested patches for spatial comparisons.
    patches={'cross5':[(1,0),(0,1),(1,1),(2,1),(1,2)],
             'square3':rectangle(3,3),'square2':rectangle(2,2),'rectangle2x3':rectangle(2,3)}
    models={name:HubbardSector(sites) for name,sites in patches.items()}
    for tau in (.1,.25,.5,1.):
        exact={}
        for name,model in models.items():
            started=time.monotonic();psi=model.exact(tau);value=model.doublon(psi,(1,1));exact[name]=value
            for steps in (1,2,4,8,16,32,64):
                approx=model.trotter(tau,steps);D=model.doublon(approx,(1,1))
                rows.append({'patch':name,'sites':model.n,'dimension':model.dimension,'tau':tau,'steps':steps,
                    'doublon_exact':value,'doublon_trotter':D,'absolute_trotter_error':abs(D-value),
                    'relative_trotter_error':abs(D-value)/value if value else None,
                    'state_infidelity':max(0.,1-abs(np.vdot(psi,approx))**2),
                    'norm_error':abs(np.linalg.norm(approx)-1),
                    'passes_reference_split_trotter_budget':abs(D-value)<=.001*.25/2,
                    'passes_actual_signal_split_trotter_budget':abs(D-value)<=.01*value*.25/2})
            print(name,'tau',tau,'D',round(value,8),'elapsed',round(time.monotonic()-started,2),flush=True)
        spatial.append({'tau':tau,'observed_site':[1,1],'reference':'square3','truncated':'cross5',
                        'reference_D':exact['square3'],'truncated_D':exact['cross5'],
                        'finite_reference_difference':abs(exact['square3']-exact['cross5']),
                        'reference_split_truncation_budget':.001*.25/2,
                        'actual_signal_split_truncation_budget':.01*exact['square3']*.25/2,
                        'status':'finite-reference comparison; outer-boundary error unbounded'})
    curves=[]
    times=np.linspace(0.,1.,41)
    for name,model in models.items():
        states=expm_multiply(-1j*model.hamiltonian,model.initial,start=0.,stop=1.,num=len(times),
                             traceA=-1j*model.potential.sum())
        curves.extend({'patch':name,'tau':float(t),'doublon':model.doublon(psi,(1,1)),
                       'norm_error':float(abs(np.linalg.norm(psi)-1))} for t,psi in zip(times,states))
    data={'trotter':rows,'truncation':spatial,'curves':curves,'systematic_budget_split':{
        'total_fraction':.25,'trotter_share':.5,'truncation_share':.5,
        'status':'explicit refinement scenario; replaces rather than adds to old Trotter reserve'},
        'scope':'exact state evolution within fixed spin sectors of patches up to nine sites; no L=4 calibration'}
    temporary=out/'small_dynamics.json.tmp';temporary.write_text(json.dumps(data,indent=2)+'\n')
    temporary.replace(out/'small_dynamics.json')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('results/refinement_v3'))
    args=parser.parse_args();run(args.out)
