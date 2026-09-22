"""Finite-patch, ordering-specific Strang error through tau=8.

Spin-sector factors speed up the repeated Gaussian hopping stages. Exact
references use the independent sparse full-sector Hamiltonian. No bulk error
bound is inferred from the finite patches.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh

from project_paths import ROOT
from small_dynamics import HubbardSector, rectangle, spin_hopping


class SplitSector(HubbardSector):
    def __init__(self, nx, ny):
        super().__init__(rectangle(nx,ny),max_dimension=900000)
        lookup={xy:i for i,xy in enumerate(self.sites)}
        self.edges=[];groups=[[],[]]
        for i,(x,y) in enumerate(self.sites):
            for xy,color in [((x+1,y),x%2),((x,y+1),y%2)]:
                if xy in lookup:
                    edge=(i,lookup[xy]);self.edges.append(edge);groups[color].append(edge)
        self.groups=groups
        self.local_potential=[self.u_site(i) for i in range(self.n)]
        self._spectra={}
        self._pairs={}

    def u_site(self,i):
        return 4*np.outer(((self.up>>i)&1)-.5,((self.down>>i)&1)-.5)

    def spectral_stage(self,group,dt,spin):
        key=(group,spin)
        if key not in self._spectra:
            masks=self.up if spin==0 else self.down
            self._spectra[key]=eigh(spin_hopping(masks,self.groups[group]).toarray())
        values,vectors=self._spectra[key]
        return (vectors*np.exp(-1j*dt*values))@vectors.T

    def bond_pairs(self,edge,spin):
        key=(edge,spin)
        if key not in self._pairs:
            masks=self.up if spin==0 else self.down
            h=spin_hopping(masks,[edge]).tocoo()
            take=h.row<h.col
            self._pairs[key]=(h.row[take],h.col[take],h.data[take])
        return self._pairs[key]

    def bond_stage(self,psi,edge,dt,spin):
        # Each fixed-particle spin sector is a collection of disjoint 2x2
        # hopping blocks; apply their exponentials without a many-body matrix.
        a,b,h=self.bond_pairs(edge,spin)
        view=psi if spin==0 else psi.T
        x=view[a].copy();y=view[b].copy()
        phase=-1j*np.sin(dt)*h[:,None]
        view[a]=np.cos(dt)*x+phase*y
        view[b]=np.cos(dt)*y+phase*x

    def trotter(self,tau,steps,algorithm='plaquette'):
        dt=tau/steps;psi=self.initial.reshape(self.shape).copy()
        if algorithm=='plaquette':
            onsite=np.exp(-.5j*dt*self.potential.reshape(self.shape))
            # Dense spin-sector factors are economical on the <=9-site patches.
            if self.n<=9:
                pink=[self.spectral_stage(0,dt/2,s) for s in (0,1)]
                gold=[self.spectral_stage(1,dt,s) for s in (0,1)]
                for _ in range(steps):
                    psi*=onsite;psi=pink[0]@psi@pink[1].T
                    psi=gold[0]@psi@gold[1].T
                    psi=pink[0]@psi@pink[1].T;psi*=onsite
                return psi.ravel()
            # The sparse implementation avoids dense 924x924 matrix products.
            return super().trotter(tau,steps)
        if algorithm!='termwise':raise ValueError(algorithm)
        # Matches build_workload: onsite[0], ordered bonds (XX, YY per spin),
        # onsite[1:], then the reverse list. XX and YY on one bond commute,
        # and their paired hopping gate conserves the fixed spin sector.
        first=np.exp(-.5j*dt*self.local_potential[0])
        rest=np.exp(-1j*dt*(self.potential.reshape(self.shape)-self.local_potential[0]))
        if self.n<=9:
            sweeps=[]
            for spin,masks in enumerate((self.up,self.down)):
                forward=np.eye(len(masks),dtype=complex)
                for edge in self.edges:self._left_bond(forward,edge,dt/2,spin)
                sweeps.append(forward)
            # Real symmetric bond Hamiltonians give reverse = forward.T.
            for _ in range(steps):
                psi*=first;psi=sweeps[0]@psi@sweeps[1].T
                psi*=rest;psi=sweeps[0].T@psi@sweeps[1];psi*=first
        else:
            for _ in range(steps):
                psi*=first
                for edge in self.edges:
                    for spin in (0,1):self.bond_stage(psi,edge,dt/2,spin)
                psi*=rest
                for edge in reversed(self.edges):
                    for spin in (1,0):self.bond_stage(psi,edge,dt/2,spin)
                psi*=first
        return psi.ravel()

    def _left_bond(self,matrix,edge,dt,spin):
        a,b,h=self.bond_pairs(edge,spin)
        x=matrix[a].copy();y=matrix[b].copy();phase=-1j*np.sin(dt)*h[:,None]
        matrix[a]=np.cos(dt)*x+phase*y;matrix[b]=np.cos(dt)*y+phase*x


def atomic_json(path,data):
    tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');tmp.replace(path)


def run(out,patches,times):
    out.mkdir(parents=True,exist_ok=True)
    path=out/'calibration.json'
    data=json.loads(path.read_text()) if path.exists() else {'rows':[],'scope':'Finite patches and local doublon; no bulk accuracy certificate',
        'hamiltonian':'OBC, half-filled Neel, U/J=4, central coordinate (1,1)',
        'step_rule':'r0=ceil(16 max(1,tau)); r0,2r0,4r0; additional refinement if coefficient spread exceeds 5 percent'}
    done={(v['patch'],v['tau'],v['algorithm']) for v in data['rows']}
    for nx,ny in patches:
        patch=f'{nx}x{ny}';model=None
        for tau in times:
            missing=[a for a in ('plaquette','termwise') if (patch,tau,a) not in done]
            if not missing:continue
            if model is None:model=SplitSector(nx,ny)
            started=time.monotonic();exact=model.exact(tau);D=model.doublon(exact,(1,1))
            print(patch,tau,'exact',D,'seconds',time.monotonic()-started,flush=True)
            for algorithm in missing:
                measurements=[];r0=int(np.ceil(16*max(1,tau)))
                for r in (r0,2*r0,4*r0):
                    started=time.monotonic();psi=model.trotter(tau,r,algorithm)
                    signed=model.doublon(psi,(1,1))-D
                    measurements.append(dict(r=r,signed_error=signed,coefficient=signed*r*r,
                                             norm_error=float(abs(np.linalg.norm(psi)-1)),seconds=time.monotonic()-started))
                    print(patch,tau,algorithm,r,'C',signed*r*r,'seconds',measurements[-1]['seconds'],flush=True)
                # The coefficient can cross zero; report absolute as well as relative agreement.
                last=[x['coefficient'] for x in measurements[-2:]]
                spread=abs(last[1]-last[0])/max(abs(last[0]),abs(last[1]),1e-10)
                if spread>.05:
                    r=8*r0;started=time.monotonic();psi=model.trotter(tau,r,algorithm)
                    signed=model.doublon(psi,(1,1))-D
                    measurements.append(dict(r=r,signed_error=signed,coefficient=signed*r*r,
                                             norm_error=float(abs(np.linalg.norm(psi)-1)),seconds=time.monotonic()-started))
                    last=[x['coefficient'] for x in measurements[-2:]]
                    spread=abs(last[1]-last[0])/max(abs(last[0]),abs(last[1]),1e-10)
                data['rows'].append(dict(patch=patch,sites=nx*ny,dimension=model.dimension,tau=tau,algorithm=algorithm,
                    exact_doublon=D,exact_norm_error=float(abs(np.linalg.norm(exact)-1)),measurements=measurements,
                    coefficient=max(abs(v) for v in last),relative_coefficient_spread=spread,
                    converged=spread<=.05))
                data['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
                atomic_json(path,data)
    return data


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--patches',nargs='+',default=['2x2','3x3'])
    p.add_argument('--times',nargs='+',type=float,default=[.1,.25]+list(np.arange(.5,8.01,.5)))
    p.add_argument('--out',type=Path,default=ROOT/'results/discussion_v3')
    args=p.parse_args();run(args.out,[tuple(map(int,x.split('x'))) for x in args.patches],args.times)
