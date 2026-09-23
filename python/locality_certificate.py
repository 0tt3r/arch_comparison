"""A conservative centered truncation certificate in the even-fermion algebra.

Remove onsite terms by a product interaction picture. A two-spin hopping bond
has norm 2J. The commutator recursion starts on at most four bonds; each bond
meets at most seven bonds (including itself). For centered O=D-I/2, ||O||=1/2,
Duhamel and the LR recursion bound the omitted boundary by
  (2/7) sum_{n>=R+1} (28*tau)^n/n!.
R is the graph distance retained around the observable. This upper bound is
state independent and is not a prediction of physical propagation speed.
"""
import hashlib
import json
from math import ceil, exp, lgamma, log, log1p

from project_paths import ROOT

OUT=ROOT/'results/frontier_v1'


def log_bound(radius,tau):
    if radius<0 or tau<0:raise ValueError('Nonnegative radius and time required')
    if tau==0:return float('-inf')
    x=28*tau;n=radius+1
    if n+1<=x:return 0.
    # The first omitted term times a geometric majorant; no tail underflow.
    value=log(2/7)+n*log(x)-lgamma(n+1)-log1p(-x/(n+1))
    return min(0.,value+1e-12)


def sufficient_radius(tau,epsilon):
    if not 0<epsilon<1:raise ValueError('An absolute allowance between zero and one is required')
    lo=-1;hi=max(1,ceil(80*tau+20))
    while log_bound(hi,tau)>log(epsilon):hi*=2
    while hi-lo>1:
        mid=(lo+hi)//2
        if log_bound(mid,tau)<=log(epsilon):hi=mid
        else:lo=mid
    return hi


def run():
    config=json.loads((ROOT/'configs/frontier_study.json').read_text())['locality'];rows=[]
    for tau in config['times']:
        for epsilon in config['allowances']:
            radius=sufficient_radius(tau,epsilon)
            rows.append(dict(tau=tau,spatial_allowance=epsilon,radius=radius,
                square_side=2*radius+1,square_sites=(2*radius+1)**2,
                diamond_sites=1+2*radius*(radius+1),certified_error_bound=exp(log_bound(radius,tau)),
                preceding_radius_bound=exp(log_bound(radius-1,tau)) if radius else None,
                scope='Sufficient state-independent patch under a deliberately loose commutator bound; not a necessary radius or measured bulk convergence'))
    curves=[dict(tau=tau,radius=r,bound=exp(log_bound(r,tau))) for tau in config['times'] for r in range(config['max_radius']+1)]
    sources=['python/locality_certificate.py','configs/frontier_study.json']
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'locality_certificate.json').write_text(json.dumps(dict(rows=rows,curves=curves,
        derivation='Centered doublon norm 1/2; hopping bond norm 2; four first edges and at most seven successors; Duhamel/LR tail (2/7) sum_{n>=R+1}(28 tau)^n/n!. Geometric majorant used numerically.',
        source_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources}),indent=2)+'\n')
    for row in rows:print(row['tau'],row['spatial_allowance'],row['radius'],row['certified_error_bound'])


if __name__=='__main__':run()
