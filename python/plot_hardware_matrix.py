"""Capability of every architecture on the four hardware types (2D grid or all-to-all; uniform or non-uniform noise).

Hubbard (tau = L/4, full T-source menu) and QUOPS (standard threshold, FT with cultivation or factories), at p = 1e-3.
Sources: rounds v5 (grid, uniform), v3 (grid, non-uniform), v6 (all-to-all, uniform), v7 (all-to-all, non-uniform);
QUOPS results quops_sd6_v1, quops_v1, quops_a2a_sd6_v1, quops_a2a_v1 and their sensitivity runs (factory menu).
"""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from discussion_model import capacity
from plot_capacity import COLORS
from project_paths import ROOT

LABELS = {'nisq_effective': 'Physical (assumed mitigation)', 'star': 'STAR', 'hybrid': 'Hybrid',
          'ft': 'Surface-code FT', 'pinnacle': 'Pinnacle (qLDPC)', 'nisq': 'Physical'}
PANELS = [('2D grid, uniform noise', 'v5', 'quops_sd6_v1', 'quops_sensitivity_sd6_v1'),
          ('2D grid, non-uniform noise', 'v3', 'quops_v1', 'quops_sensitivity_v1'),
          ('All-to-all, uniform noise', 'v6', 'quops_a2a_sd6_v1', 'quops_sensitivity_a2a_sd6_v1'),
          ('All-to-all, non-uniform noise', 'v7', 'quops_a2a_v1', 'quops_sensitivity_a2a_v1')]
P, EPS, CONFIG = 1e-3, .01, 'calibrated_factories'
Q = np.geomspace(10, 1e9, 500)


def read(path):
    return json.loads((ROOT/path).read_text())


def hubbard():
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    arms = ('nisq_effective', 'star', 'hybrid', 'ft', 'pinnacle')
    for ax, (title, rnd, _, _) in zip(axes.flat, PANELS):
        rows = [v for v in read(f'results/optimization_round_{rnd}/costs.json')
                if v['configuration'] == CONFIG and v['family'] == 'tau_L_over_4' and v['p'] == P
                and v['absolute_error'] == EPS and v['scenario'] == 'component']
        for arm in arms:
            group = [v for v in rows if v['arm'] == arm]
            if not group:
                continue
            y = np.array([capacity(group, q) for q in Q], float)
            ax.step(Q, y, where='post', color=COLORS['nisq' if arm.startswith('nisq') else arm], lw=2,
                    ls='--' if arm.startswith('nisq') else '-')
        ax.set(xscale='log', xlim=(10, 1e9), ylim=(0, 33), yticks=[4, 8, 16, 24, 32], title=title)
        ax.grid(alpha=.2)
    for ax in axes[1]:
        ax.set_xlabel('Available physical qubits Q')
    for ax in axes[:, 0]:
        ax.set_ylabel('Largest lattice side L (τ = L/4)')
    handles = [Line2D([], [], color=COLORS['nisq' if a.startswith('nisq') else a], lw=2,
                      ls='--' if a.startswith('nisq') else '-', label=LABELS[a]) for a in arms]
    fig.legend(handles=handles, loc='lower center', ncol=5, frameon=False, fontsize=9, bbox_to_anchor=(.5, -.01))
    fig.suptitle(f'Hubbard dynamics on four hardware types, p = {P:g}, ε = {EPS:g}, one week')
    fig.tight_layout(rect=(0, .04, 1, 1))
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT/f'figures/hardware_matrix_hubbard.{ext}', dpi=180, bbox_inches='tight')
    plt.close(fig)


def quops():
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    arms = ('nisq', 'star', 'ft', 'pinnacle')
    for ax, (title, _, main, sens) in zip(axes.flat, PANELS):
        # Capability from the saved per-shape minima, so every panel shares one qubit grid.
        costs = [v for v in read(f'results/{main}/costs.json') if v['p'] == P and v['mode'] == 'standard']
        factory = sorted((v['Q'], v['score']) for v in read(f'results/{sens}/capability.json')
                         if v['variant'] == 'calibrated_factories' and v['arm'] == 'ft' and v['p'] == P)
        for arm in arms:
            rows = [v for v in costs if v['arm'] == arm]
            if not rows:
                continue
            y = np.array([max([v['s'] for v in rows if v['qubits'] <= q], default=0) for q in Q], float)
            if arm == 'ft':
                y = np.maximum(y, [max([sc for qf, sc in factory if qf <= q*(1+1e-9)], default=0) for q in Q])
            y[y == 0] = np.nan
            ax.step(Q, y, where='post', color=COLORS[arm], lw=2)
        ax.set(xscale='log', yscale='log', xlim=(10, 1e9), ylim=(10, 1e10), title=title)
        ax.grid(alpha=.2)
    for ax in axes[1]:
        ax.set_xlabel('Available physical qubits Q')
    for ax in axes[:, 0]:
        ax.set_ylabel('QUOPS score')
    handles = [Line2D([], [], color=COLORS[a], lw=2, label=LABELS[a]) for a in arms]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(.5, -.01))
    fig.suptitle(f'Random circuits (QUOPS) on four hardware types, p = {P:g}, one week')
    fig.tight_layout(rect=(0, .04, 1, 1))
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT/f'figures/hardware_matrix_quops.{ext}', dpi=180, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    hubbard()
    quops()
