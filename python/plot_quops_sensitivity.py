"""Angle-sweep figure from results/quops_sensitivity_v1 (STAR) and results/quops_v1 (FT, physical)."""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from plot_capacity import COLORS
from project_paths import ROOT

import quops_study as _q
OUT = ROOT/_q.settings().get('sensitivity_output', 'results/quops_sensitivity_v1')


def main():
    t = json.loads((OUT/'capability.json').read_text())
    base = json.loads((_q.OUT/'capability.json').read_text())
    variants = sorted({x['variant'] for x in t if x['variant'].startswith('angle_')}, key=lambda v: float(v.split('_')[1]))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    markers = {1e-3: 'o', 1e-4: 's'}
    fits = {}
    for p in (1e-3, 1e-4):
        mean, sat = [], []
        for v in variants:
            s = max(x['score'] for x in t if (x['variant'], x['arm'], x['p']) == (v, 'star', p))
            mean.append(float(v.split('_')[1])/2); sat.append(s)
        mean, sat = np.array(mean), np.array(sat)
        a1.loglog(mean, sat, markers[p], color=COLORS['star'], ms=7, label=f'STAR saturation, $p={p:g}$')
        k = np.median((sat*mean*p)[mean < .02])
        fits[p] = k
        x = np.geomspace(1e-3, .3, 50)
        a1.loglog(x, k/(p*x), color=COLORS['star'], lw=1, ls='--')
    a1.axvspan(.008, .035, color='#cccccc', alpha=.4, lw=0)   # mean Trotter angle at r = 8-30 steps per unit time
    a1.text(.0085, 2.2e6, 'Hubbard\nTrotter angles', fontsize=8, color='#555555', va='top')
    a1.axvline(np.pi/16, color='#999999', ls=':', lw=1)
    a1.text(np.pi/16*1.05, 2.2e6, 'QUOPS', fontsize=8, color='#555555', va='top')
    a1.set(xlabel=r'Mean injected angle $\langle|\theta|\rangle$ (rad)', ylabel='Largest passing size, any Q',
           title=r'STAR saturation follows $s\approx0.41/(p\langle|\theta|\rangle)$', ylim=(1e3, 3e6))
    a1.legend(frameon=False, fontsize=8, loc='lower left'); a1.grid(alpha=.2, which='both')
    p = 1e-3
    shades = np.linspace(.35, 1, len(variants))
    for v, alpha in zip(variants, shades):
        r = [x for x in t if (x['variant'], x['arm'], x['p']) == (v, 'star', p)]
        Q = np.array([x['Q'] for x in r]); y = np.array([x['score'] for x in r], float); y[y == 0] = np.nan
        a2.step(Q, y, where='post', color=COLORS['star'], alpha=alpha, lw=1.8,
                label=rf'STAR, $\langle|\theta|\rangle={float(v.split("_")[1])/2:.3g}$')
    b = [x for x in base if x['p'] == p and x['mode'] == 'standard']
    Q = np.array([x['Q'] for x in b])
    for arm, label in (('ft', 'FTQC'), ('nisq', 'Physical qubits')):
        y = np.array([x[arm] for x in b], float); y[y == 0] = np.nan
        a2.step(Q, y, where='post', color=COLORS[arm], lw=2.2, label=label)
    a2.set(xscale='log', yscale='log', xlim=(min(1e3, _q.settings().get('q_min', 1e2)), 1e9), ylim=(1e2, 1e8), xlabel='Available physical qubits Q',
           ylabel='Largest passing size s', title=r'Random circuits, $p=10^{-3}$, standard threshold')
    a2.legend(frameon=False, fontsize=7.5, loc='upper left'); a2.grid(alpha=.2)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT/f"figures/quops_angle_sweep{_q.settings().get('figure_suffix', '')}.{ext}", dpi=180, bbox_inches='tight')
    (OUT/'angle_law.json').write_text(json.dumps({str(k): v for k, v in fits.items()}, indent=1)+'\n')
    print(fits)


if __name__ == '__main__':
    main()
