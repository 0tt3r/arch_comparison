"""Exact Trotter coefficients on strips versus length (results/strip_dynamics_v1)."""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from project_paths import ROOT

WIDTH = {2: ('#2a6fb0', 'o'), 3: ('#e0892b', 's'), 4: ('#3a9a5b', 'D')}
STEP = {1.: 16, 2.: 32}


def main():
    rows = json.loads((ROOT/'results/strip_dynamics_v1/coefficients.json').read_text())
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    for ax, tau in zip(axes, (1., 2.)):
        for algorithm, ls in (('termwise', '-'), ('plaquette', '--')):
            for width, (color, marker) in WIDTH.items():
                group = sorted((v for v in rows if v['tau'] == tau and v['algorithm'] == algorithm
                                and v['ny'] == width and v['r'] == STEP[tau] and v['nx'] >= width),
                               key=lambda v: v['nx'])
                if group:
                    ax.plot([v['nx'] for v in group], [v['central'] for v in group], ls=ls, marker=marker,
                            ms=4, lw=1.4, color=color,
                            label=f'width {width}' if algorithm == 'termwise' else None)
        ax.set(xlabel='Patch length $L_x$', title=rf'$U/J=4$, $\tau={tau:g}$', ylim=(0, None))
        ax.grid(alpha=.2)
    axes[0].set_ylabel(r'central $|D_r-D|\,r^2$')
    axes[0].plot([], [], 'k-', label='termwise')
    axes[0].plot([], [], 'k--', label='plaquette')
    axes[0].legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT/f'figures/strip_dynamics.{ext}', dpi=180, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
