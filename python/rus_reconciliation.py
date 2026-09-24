"""Reconcile raw-RUS error accumulation with Toshio et al., PRX 15, 021057 (2025), Sec. IV and App. E.

Their model: a resource prepared for angle a > 0 is wrong with probability
P_L(a) = a^(2/k) P_ud, P_ud = k p / 15, and then over-rotates by
Delta(a) = -a^(1-2/k). Completion at trial K composes the channels
E_{2^(K-1) theta} o E_{-2^(K-2) theta} o ... o E_{-theta} (their Eq. 19), with
angles above pi/8 reduced by a logical S (wrap to |a| <= pi/8), and the averaged
channel is sum_K 2^-K E^K. The only undetermined input is Delta at a negative
angle, i.e. the coherent error of a failed trial:

- conjugate: a failed teleportation implements minus the resource's actual
  angle, so its error is -Delta(|a|) (physical);
- same-sign: Delta evaluated at |a|, the positive-angle formula.

We compute the worst-case error sqrt(x^2 + y^2) of the averaged channel under
both, as alpha' = eps/(|theta| p), and the coherent-cancellation value
alpha_RUS of their Eq. (24) as a check on the model (their Fig. 11: about 0.4 k
at small angles).
"""
import argparse
import json
from math import pi, sin

import numpy as np

from project_paths import ROOT

OUT = ROOT/'results/rus_reconciliation_v1'


def wrap(a):
    """Signed analog angle after removing multiples of pi/4 with logical S gates."""
    return (a + pi/8) % (pi/4) - pi/8


def trial(a, k, p):
    """(P_L, Delta) for a resource prepared at |a|; Delta refers to the positive-angle formula."""
    m = abs(a)
    if m == 0:
        return 0., 0.
    return m**(2/k)*k*p/15, -m**(1-2/k)


def averaged(theta, k, p, convention, cap=60):
    """x, y of sum_K 2^-K E^K and the PCEC effective rate of their Eq. (24)."""
    x = y = pcec = 0.
    for K in range(1, cap+1):
        xk = yk = ek = 0.
        for n in range(K):
            a = wrap(2**n*theta)
            P, D = trial(a, k, p)
            s = np.sign(a) or 1.
            success = n == K-1
            if convention == 'conjugate':
                d = s*D if success else -s*D     # a failure implements -(actual resource angle)
            else:
                d = s*D if success else s*D      # same-sign: failures repeat the success error
            xk += P*sin(d)**2
            yk += .5*P*sin(2*d)
            ek += 2*P*sin(D)**2                  # PCEC removes the coherent part and doubles dephasing
        w = 2.**-K
        x, y, pcec = x + w*xk, y + w*yk, pcec + w*ek
    return x, y, pcec


def table(ks=(2, 3, 4, 6), thetas=np.geomspace(1e-4, pi/8, 41), p=1e-3):
    rows = []
    for k in ks:
        for theta in thetas:
            row = dict(k=k, theta=float(theta))
            for convention in ('conjugate', 'same_sign'):
                x, y, pcec = averaged(theta, k, p, convention)
                row[convention] = dict(alpha=float(np.hypot(x, y)/(theta*p)), x=x, y=y)
            row['alpha_pcec'] = float(pcec/(theta*p))
            rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    rows = table()
    for r in rows[::8]:
        print(f"k={r['k']} theta={r['theta']:.1e}  PCEC/k={r['alpha_pcec']/r['k']:.3f}  "
              f"raw conjugate/k={r['conjugate']['alpha']/r['k']:.3f} (|y|/x={abs(r['conjugate']['y'])/r['conjugate']['x']:.1e})  "
              f"raw same-sign/k={r['same_sign']['alpha']/r['k']:.3f}")
    if args.write:
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT/'alpha.json').write_text(json.dumps(rows, indent=1)+'\n')


if __name__ == '__main__':
    main()
