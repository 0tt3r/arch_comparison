"""QUOPS capability figure from saved results (results/quops_v1)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from plot_capacity import COLORS
from project_paths import ROOT
from quops_study import MODES, OUT, read, settings

ARMS = ('nisq', 'star', 'ft')

LABELS = {'nisq': 'Physical qubits', 'star': 'STAR', 'ft': 'FTQC'}
STYLE = {'standard': '-', 'mitigated': '--'}
# Measured physical-qubit scores and challenge targets, Proctor et al. (arXiv:2609.12146).
MEASURED = [('Willow', 216), ('Helios-1', 1504)]
TARGETS = [('RSA-2048', 2.5e8), ('FeMoco', 3.4e8)]


def published_factories(p):
    """FT with Litinski-type factories re-optimized under the calibrated law, alongside cultivation."""
    rows = read(ROOT/settings().get('sensitivity_output', 'results/quops_sensitivity_v1')/'capability.json')
    return {v['Q']: v['score'] for v in rows if v['variant'] == 'calibrated_factories' and v['arm'] == 'ft' and v['p'] == p}


def main():
    table = read(OUT/'capability.json')
    cfg = settings()
    fig, axes = plt.subplots(1, len(cfg['physical_errors']), figsize=(12, 4.8), sharey=True)
    for ax, p in zip(np.atleast_1d(axes), cfg['physical_errors']):
        for arm in ARMS:
            for mode in MODES:
                rows = [v for v in table if v['p'] == p and v['mode'] == mode]
                Q = np.array([v['Q'] for v in rows])
                y = np.array([v[arm] for v in rows], float)
                y[y == 0] = np.nan
                ax.step(Q, y, where='post', color=COLORS[arm], ls=STYLE[mode], lw=2 if mode == 'standard' else 1.4)
                if arm == 'ft' and mode == 'standard':
                    menu = published_factories(p)
                    better = np.array([menu[v['Q']] for v in rows], float)
                    better[better == 0] = np.nan
                    ax.fill_between(Q, y, np.fmax(better, y), step='post', color=COLORS['ft'], alpha=.15, lw=0)
        for name, value in MEASURED:
            ax.axhline(value, color='#999999', lw=.8, ls=':')
            ax.text(1.3*settings().get('q_min', 1e2), value*1.15, f'{name} (measured)', fontsize=7.5, color='#666666')
        for (name, value), offset in zip(TARGETS, (0.55, 1.2)):
            ax.axhline(value, color='#bbbbbb', lw=.8, ls='-.')
            ax.text(1.3*settings().get('q_min', 1e2), value*offset, f'{name} target', fontsize=7.5, color='#666666',
                    va='top' if offset < 1 else 'bottom')
        ax.set(xscale='log', yscale='log', xlim=(settings().get('q_min', 1e2), 1e9), ylim=(1e1, 1e12),
               title=rf'$p={p:g}$', xlabel='Available physical qubits Q')
        ax.grid(alpha=.2)
    np.atleast_1d(axes)[0].set_ylabel('QUOPS score (largest passing size s)')
    handles = [Line2D([], [], color=COLORS[a], lw=2, label=LABELS[a]) for a in ARMS]
    handles += [Line2D([], [], color='#555555', ls='-', label='standard threshold'),
                Line2D([], [], color='#555555', ls='--', label='1/γ² throughput scenario')]
    from matplotlib.patches import Patch
    handles.append(Patch(color=COLORS['ft'], alpha=.2, label='FTQC with distillation factories'))
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(.5, -.02))
    fig.suptitle('QUOPS capability on a nearest-neighbor array, one week')
    fig.text(.5, -.08, 'Standard: mean polarization ≥ 1/√e. Dashed: declared scenario in which a 1/γ² throughput penalty must fit in one week.\n'
             'Routing and lattice-surgery timing are declared models; measured scores are at each device\'s own error rates.',
             ha='center', fontsize=8.5, color='#555555')
    fig.tight_layout(rect=(0, .06, 1, 1))
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT/f"figures/quops_capability{settings().get('figure_suffix', '')}.{ext}", dpi=180, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
