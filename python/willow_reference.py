"""Matched Willow hopping order and finite-sector reference calculations.

The hopping order follows Appendix D2 of 2510.26845v1. Routing permutations
are tracked in a fermionic frame; the reverse half restores the original frame.
The interacting production task is U=4, zero flux, center hole on 4x4.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from project_paths import ROOT
from small_dynamics import rectangle
from willow_archive import sha, save
from willow_benchmark import load_reduced, published_arrays, combine_paulis

OUT = ROOT/'results/willow_reference_v1'


def hopping_order(nx, ny):
    """Each physical nearest-neighbor hopping appears once, in paper order."""
    sites=rectangle(nx,ny); index={xy:i for i,xy in enumerate(sites)}
    live={xy:i for i,xy in enumerate(sites)}
    edges={tuple(sorted((i,index[xy]))) for i,(x,y) in enumerate(sites)
           for xy in ((x+1,y),(x,y+1)) if xy in index}
    done=set(); order=[]
    def hop(a,b):
        e=tuple(sorted((a,b)))
        if e in edges and e not in done:
            order.append(e);done.add(e)
    def vertical():
        for y in range(ny-1):
            x=nx-1 if y%2==0 else 0
            hop(live[x,y],live[x,y+1])
    def swaps(offset):
        for y in range(ny):
            for x in range(offset,nx-1,2):
                a,b=(x,y),(x+1,y)
                hop(live[a],live[b]);live[a],live[b]=live[b],live[a]
    vertical()
    for _ in range(2*nx):
        if done==edges:break
        swaps(0)
        if nx%2==0:vertical()
        if done==edges:break
        swaps(1);vertical()
    assert done==edges, 'Swap network omitted a hopping term'
    return order


def free_circuit(L, flux, occupations, tau, steps):
    sites=rectangle(L,L);order=hopping_order(L,L)
    init=np.zeros((L*L,2))
    for v in occupations:init[sites.index(tuple(v['site'])),0 if v['spin']=='up' else 1]=v['n']
    R=np.eye(L*L,dtype=complex)
    co=np.cos(tau/(2*steps));si=np.sin(tau/(2*steps))
    for _ in range(steps):
        for i,j in order+order[::-1]:
            x,y=sites[i];xx,yy=sites[j]
            phase=np.exp(1j*flux*x*(yy-y)) if x==xx else 1.
            a,b=R[i].copy(),R[j].copy()
            R[i]=co*a+1j*si*phase*b
            R[j]=co*b+1j*si*np.conjugate(phase)*a
    density=np.abs(R)**2@init
    z=np.empty((L,L,2))
    for i,(x,y) in enumerate(sites):z[x,y]=1-2*density[i]
    return combine_paulis(z,z[:,:,0]*z[:,:,1])


def circuit_audit():
    data=load_reduced();rows=[]
    for c in data['cases']:
        if c['U']!=0:continue
        arrays=published_arrays(c)
        for s in c['shots']:
            tau=s['tau'];v=free_circuit(c['L'],0 if c['flux']=='0' else np.pi,c['initial_occupations'],tau,c['steps'])
            p=arrays['exact',tau];L=c['L']
            raw=max(abs(np.array(v['local_z'])-p['local_z']))
            z=np.array(v['local_z']).reshape(L,L,2)
            sym=max(abs(((z+z.transpose(1,0,2))/2).ravel()-p['local_z']))
            d=np.array(v['local_d']).reshape(L,L)
            dsym=max(abs(((d+d.T)/2).ravel()-p['local_d']))
            rows.append(dict(L=L,flux=c['flux'],steps=c['steps'],tau=tau,
                max_Z_difference=float(raw),max_symmetrized_Z_difference=float(sym),
                max_symmetrized_D_difference=float(dsym)))
    result=dict(rows=rows,maximum_unsymmetrized_Z_difference=max(v['max_Z_difference'] for v in rows),
        maximum_symmetrized_Z_difference=max(v['max_symmetrized_Z_difference'] for v in rows),
        maximum_symmetrized_D_difference=max(v['max_symmetrized_D_difference'] for v in rows))
    assert len(rows)==288
    assert result['maximum_symmetrized_Z_difference']<1e-12
    # Archived ZZ-derived doublons have approximately 1e-8 roundoff, while Z
    # itself reproduces to machine precision. Keep both residuals visible.
    assert result['maximum_symmetrized_D_difference']<1e-7
    result['status']='PASS: paper hopping order reproduces all archived U=0 circuit references after their diagonal averaging'
    save(OUT/'circuit_audit.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--circuit-audit',action='store_true')
    args=parser.parse_args()
    if args.circuit_audit:circuit_audit()
