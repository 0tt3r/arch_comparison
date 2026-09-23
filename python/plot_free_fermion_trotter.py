"""Free-fermion Trotter coefficient versus lattice size (results/free_fermion_trotter_v1)."""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from project_paths import ROOT

COLORS = {1.: '#2a6fb0', 2.: '#e0892b', 4.: '#3a9a5b'}


def main():
    rows = json.loads((ROOT/'results/free_fermion_trotter_v1/coefficients.json').read_text())['rows']
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    for ax, algorithm in zip(axes, ('termwise', 'plaquette')):
        for tau, color in COLORS.items():
            group = [v for v in rows if v['algorithm'] == algorithm and v['family'] == 'fixed'
                     and v['tau'] == tau and v['r'] == 32]
            L = [v['L'] for v in group]
            ax.plot(L, [max(v['central_coefficient'], 1e-4) for v in group], 'o-', ms=3.5, lw=1.5,
                    color=color, label=rf'$\tau={tau:g}$')
            ax.axvline(4*tau+2, color=color, lw=.8, ls=':')
        ax.set(yscale='log', ylim=(1e-4, 10), xlabel='Lattice size $L$', title=f'{algorithm}, $U=0$')
        ax.grid(alpha=.2)
    axes[0].set_ylabel(r'central $|D_r-D|\,r^2$ ($r=32$)')
    axes[0].legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT/f'figures/free_fermion_trotter.{ext}', dpi=180, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
