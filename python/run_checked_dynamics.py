"""Finish the finite-lattice step menu, refining only if the tail has not settled."""
import json
from large_dynamics import ROOT, solve
from small_dynamics import rectangle


def run():
    out=ROOT/'results/dynamics_4x4/tau1'
    exact=solve(4,4,1.,'exact',0,out,16)
    i=rectangle(4,4).index((2,2));D=exact['result']['doublons'][i]
    convergence=[]
    for algorithm in ('plaquette','termwise'):
        records=[solve(4,4,1.,algorithm,r,out,16) for r in (8,16,32)]
        def coefficients():
            return [abs(v['result']['doublons'][i]-D)*v['contract']['r']**2 for v in records[-2:]]
        a,b=coefficients();spread=abs(a-b)/max(a,b,1e-10)
        if spread>.05:
            records.append(solve(4,4,1.,algorithm,64,out,16))
            a,b=coefficients();spread=abs(a-b)/max(a,b,1e-10)
        if spread>.05:raise RuntimeError(f'{algorithm}: still unconverged; do not publish a fitted coefficient')
        convergence.append({'algorithm':algorithm,'steps':[v['contract']['r'] for v in records],
                            'relative_tail_coefficient_spread':spread})
    from large_dynamics import atomic_json
    atomic_json(ROOT/'results/dynamics_4x4/convergence.json',{'rows':convergence,
        'criterion':'Two finest |D_r-D_exact| r^2 coefficients agree within 5%; all r=8,16,32 measured',
        'scope':'4x4 tau=1 central coordinate (2,2)'})


if __name__=='__main__':run()
