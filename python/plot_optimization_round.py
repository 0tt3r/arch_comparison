"""Figures and summary for the optimization round, from saved results only."""
from math import log10

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

from discussion_model import capacity
from optimization_round import ENCODED, OUT, read, settings
from plot_capacity import COLORS
from project_paths import ROOT

LABELS = {'star': 'STAR', 'hybrid': 'Hybrid', 'ft': 'FTQC',
          'nisq_effective': 'NISQ: assumed effective mitigation', 'nisq_pec': 'NISQ: two-qubit PEC reference'}
NISQ_STYLE = {'nisq_effective': (COLORS['nisq'], '--'), 'nisq_pec': ('#8c622c', ':')}
CONFIG_LABELS = {'v1_menu': 'v1 menu', 'raw_rus': '+ raw RUS', 'lanes': '+ lanes', 'ft_delivery': '+ FT delivery',
                 'distances': '+ separate distances', 'richardson': '+ Richardson', 'all': 'all five'}


def save(fig, name):
    for suffix in ('pdf', 'png'):
        fig.savefig(ROOT/f'figures/{name}.{suffix}', dpi=180, bbox_inches='tight',
                    metadata={'CreationDate': None, 'ModDate': None} if suffix == 'pdf' else {})
    plt.close(fig)


def select(rows, **keys):
    return [v for v in rows if all(v[k] == x for k, x in keys.items())]


def staircase(group, Q):
    return np.array([capacity(group, q) for q in Q], float)


def grid(group):
    thresholds = [10**v['log10_Q_min'] for v in group if v['log10_Q_min'] is not None and 2 <= v['log10_Q_min'] <= 9]
    return np.unique(np.r_[np.geomspace(1e2, 1e9, 350), thresholds])


def plot_family(rows, family, name, title):
    cfg = settings()
    sizes = cfg['families'][family]['sizes']
    proportional = cfg['families'][family]['kind'] == 'proportional'
    fig, axes = plt.subplots(2, 2, figsize=(11.7, 8.2), sharex=True, sharey=True)
    for i, eps in enumerate(cfg['precisions']):
        for j, p in enumerate(cfg['physical_errors']):
            ax = axes[i, j]
            group = select(rows, family=family, p=p, absolute_error=eps)
            Q = grid(group)
            for arm in ENCODED:
                new = staircase(select(group, configuration='all', arm=arm, scenario='component'), Q)
                old = staircase(select(group, configuration='v1_menu', arm=arm, scenario='component'), Q)
                if arm != 'ft':
                    zero = staircase(select(group, configuration='all', arm=arm, scenario='zero_growth'), Q)
                    ax.fill_between(Q, np.minimum(new, zero), np.maximum(new, zero), step='post',
                                    color=COLORS[arm], alpha=.14, lw=0)
                ax.step(Q, old, where='post', color=COLORS[arm], lw=1, alpha=.55, ls=(0, (4, 2)))
                ax.step(Q, new, where='post', color=COLORS[arm], lw=2.2)
            for arm, (color, style) in NISQ_STYLE.items():
                ax.step(Q, staircase(select(group, configuration='all', arm=arm), Q), where='post',
                        color=color, ls=style, lw=1.6)
            top = max(sizes)
            ticks = [s for s in sizes if s % 8 == 0 or s == sizes[0]]
            ax.set(xscale='log', xlim=(1e2, 1e9), ylim=(0, top+1), yticks=ticks,
                   title=rf'$p={p:g},\ \ \epsilon={eps:g}$')
            if proportional:
                ax.set_yticklabels([f'{s} (τ={s*cfg["families"][family]["scale"]:g})' for s in ticks])
            ax.grid(alpha=.2)
            if i == 1:
                ax.set_xlabel('Available physical qubits Q')
    fig.supylabel('Largest lattice side L' + ('  (evolution time τ = L/4)' if proportional else ''), x=.035)
    handles = [Line2D([], [], color=COLORS[a], lw=2.2, label=LABELS[a]+' · all five changes') for a in ENCODED]
    handles += [Line2D([], [], color=NISQ_STYLE[a][0], ls=NISQ_STYLE[a][1], lw=1.6, label=LABELS[a]) for a in NISQ_STYLE]
    handles += [Line2D([], [], color='#666666', lw=1, ls=(0, (4, 2)), label='Same arm · integrated_v1 menu'),
                Patch(color='#999999', alpha=.3, label='STAR/hybrid zero-growth sensitivity')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .02), ncol=3, frameon=False, fontsize=9)
    fig.suptitle(title, y=.995)
    fig.text(.5, -.03, 'One week · capacity is the largest tested even L within Q · bands are scenarios, not confidence intervals.\n'
             'Lanes, ports and transfers are declared delivery models; NISQ effective bias is assumed. Interacting bulk convergence is not established.',
             ha='center', fontsize=8.5, color='#555555')
    fig.subplots_adjust(bottom=.2, hspace=.25, wspace=.08)
    save(fig, name)


def plot_leaders(table):
    cfg = settings()
    configs = list(cfg['configurations'])
    budgets = cfg['qubit_budgets']
    panels = [(f, p, e) for f in cfg['families'] for p in cfg['physical_errors'] for e in cfg['precisions']]
    fig, axes = plt.subplots(2, 4, figsize=(15.5, 7.4), sharey=True)
    for ax, (family, p, eps) in zip(axes.flat, panels):
        for y, config in enumerate(configs):
            for x, Q in enumerate(budgets):
                v = next(t for t in table if (t['configuration'], t['scenario'], t['family'], t['p'], t['epsilon'], t['Q'])
                         == (config, 'component', family, p, eps, Q))
                leader = v['encoded_leader']
                if leader in COLORS:
                    ax.add_patch(plt.Rectangle((x, y), 1, 1, color=COLORS[leader], alpha=.85, lw=2, ec='white'))
                    ax.text(x+.5, y+.5, str(v[leader]), ha='center', va='center', fontsize=7.5, color='white')
                elif leader.startswith('tie'):
                    ax.add_patch(plt.Rectangle((x, y), 1, 1, fc='#d9d9d9', ec='white', lw=2, hatch='///'))
                    ax.text(x+.5, y+.5, str(max(v[a] for a in ENCODED)), ha='center', va='center', fontsize=7.5, color='#333333')
                else:
                    ax.add_patch(plt.Rectangle((x, y), 1, 1, fc='#f4f4f4', ec='white', lw=2))
        ax.set(xlim=(0, len(budgets)), ylim=(len(configs), 0), xticks=np.arange(len(budgets))+.5,
               yticks=np.arange(len(configs))+.5)
        ax.set_xticklabels([f'$10^{{{int(log10(q))}}}$' for q in budgets], fontsize=8)
        ax.set_yticklabels([CONFIG_LABELS[c] for c in configs], fontsize=8.5)
        label = 'τ = 1' if family == 'tau_1' else 'τ = L/4'
        ax.set_title(f'{label} · p = {p:g} · ε = {eps:g}', fontsize=9.5)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel('Available physical qubits Q', fontsize=9)
    handles = [Patch(color=COLORS[a], alpha=.85, label=LABELS[a]) for a in ENCODED]
    handles += [Patch(fc='#d9d9d9', hatch='///', label='Tie'), Patch(fc='#f4f4f4', label='No encoded arm fits')]
    fig.legend(handles=handles, loc='lower center', ncol=5, frameon=False, fontsize=9, bbox_to_anchor=(.5, -.01))
    fig.suptitle('Which change moves the encoded capacity leader? Cells show the leader and its lattice side L', y=1.0)
    fig.tight_layout(rect=(0, .04, 1, .97))
    save(fig, 'optimization_round_leaders')


def fmt(v):
    return '—' if v is None else f'{v:.2f}'


def summary(rows, table, changes, validation):
    cfg = settings()
    lines = ['# Optimization round: which choices change the winner', '',
             'This round implements the five changes agreed in [OPUS55_REVIEW.md](../../OPUS55_REVIEW.md) §0,',
             'each alone and all together, against a `v1_menu` configuration that reproduces',
             f'`results/integrated_v1` and the NISQ rows of `results/physics_sweep_v1` ({validation["v1_replays"]} replays checked).',
             'Improvement factors below are measured outputs of the declared model, not assumptions.', '',
             '[Headline figure (τ = L/4)](../../figures/optimization_round_capacity.pdf) · '
             '[Fixed τ = 1 benchmark](../../figures/optimization_round_fixed_tau.pdf) · '
             '[Leader map](../../figures/optimization_round_leaders.pdf) · [Methods](../../docs/OPTIMIZATION_ROUND.md)', '',
             '## Encoded capacity leader by configuration', '',
             'Entries give the leader and its lattice side at Q = 10⁶ and 10⁸. The time-growing family (τ = L/4) is the headline;',
             'fixed τ = 1 is retained as an implementation benchmark.', '']
    for family in cfg['families']:
        label = 'τ = L/4' if family != 'tau_1' else 'τ = 1'
        lines += [f'### {label}', '', '| p | ε | Q | ' + ' | '.join(CONFIG_LABELS[c] for c in cfg['configurations']) + ' |',
                  '|---:|---:|---:|' + '---|'*len(cfg['configurations'])]
        for p in cfg['physical_errors']:
            for eps in cfg['precisions']:
                for Q in (10**6, 10**8):
                    cells = []
                    for config in cfg['configurations']:
                        v = next(t for t in table if (t['configuration'], t['scenario'], t['family'], t['p'], t['epsilon'], t['Q'])
                                 == (config, 'component', family, p, eps, Q))
                        caps = '/'.join(str(v[a]) for a in ENCODED)
                        cells.append(f'{v["encoded_leader"]} ({caps})')
                    lines.append(f'| {p:g} | {eps:g} | 10^{int(log10(Q))} | ' + ' | '.join(cells) + ' |')
        lines += ['', 'Cells: leader (STAR/hybrid/FT sides).', '']
    counts = {}
    for c in changes:
        if c['leader_changed']:
            counts[c['configuration']] = counts.get(c['configuration'], 0) + 1
    total = len(cfg['families'])*len(cfg['physical_errors'])*len(cfg['precisions'])*len(cfg['qubit_budgets'])
    lines += ['## How often each change moves the leader', '',
              f'Out of {total} (family, p, ε, Q) cells, relative to the v1 menu:', '',
              '| Configuration | Cells with a different encoded leader |', '|---|---:|']
    for config in cfg['configurations']:
        if config != 'v1_menu':
            lines.append(f'| {CONFIG_LABELS[config]} | {counts.get(config, 0)} |')
    lines += ['', '## What each change buys at representative points', '',
              'log₁₀ Q_min for one week; — means no feasible design within the declared menu.', '',
              '| Family | p | ε | L | Arm | ' + ' | '.join(CONFIG_LABELS[c] for c in cfg['configurations']) + ' |',
              '|---|---:|---:|---:|---|' + '---:|'*len(cfg['configurations'])]
    for family, L in (('tau_1', 8), ('tau_1', 16), ('tau_L_over_4', 8), ('tau_L_over_4', 16)):
        for p in cfg['physical_errors']:
            for eps in cfg['precisions']:
                for arm in ENCODED + ('nisq_effective',):
                    cells = []
                    for config in cfg['configurations']:
                        match = select(rows, configuration=config, family=family, L=L, p=p, absolute_error=eps, arm=arm, scenario='component')
                        cells.append(fmt(match[0]['log10_Q_min']) if match else 'n/a')
                    lines.append(f'| {"τ=1" if family == "tau_1" else "τ=L/4"} | {p:g} | {eps:g} | {L} | {arm} | ' + ' | '.join(cells) + ' |')
    lines += ['', '## Selected designs with all five changes', '',
              'What the optimizer chose, which is the mechanism behind each shift.', '',
              '| Family | p | ε | L | Arm | log₁₀ Q | r | extrapolation | d | RUS | lanes | ports | factories | HWP batch | qubits/copy | s/attempt | χ |',
              '|---|---:|---:|---:|---|---:|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|']
    for family, L in (('tau_1', 8), ('tau_L_over_4', 8), ('tau_L_over_4', 16)):
        for p in cfg['physical_errors']:
            for eps in cfg['precisions']:
                for arm in ENCODED:
                    m = select(rows, configuration='all', family=family, L=L, p=p, absolute_error=eps, arm=arm, scenario='component')
                    if not m or m[0]['log10_Q_min'] is None:
                        continue
                    v = m[0]
                    lines.append(f'| {"τ=1" if family == "tau_1" else "τ=L/4"} | {p:g} | {eps:g} | {L} | {arm} | {fmt(v["log10_Q_min"])} | {v["r"]} | '
                                 f'{v["extrapolation"]} | {v.get("distance")} | {v.get("rus_mode") or "—"} | {v.get("lanes") or "—"} | '
                                 f'{v.get("T_consumption_ports")} | {v.get("factories")} | {v.get("hwp_batch_size") or "—"} | '
                                 f'{v["qubits_per_copy"]:,.0f} | {v["seconds_per_attempt"]:.3g} | {v["chi"]:.2f} |')
    lines += ['', '## Scope', '', cfg['scope'], '',
              'Raw-RUS gains rest on first-order coherent cancellation in the restricted instrument; '
              'Toshio et al. App. E finds raw RUS slightly worse than PCEC. See [RAW_RUS_NOTE.md](../../docs/RAW_RUS_NOTE.md).', '',
              'Both encoded arms use the reference surface-code law (A = 0.1, Λ = 10 at p = 10⁻³). The Stim calibration in '
              '[SURFACE_LAW.md](../../docs/SURFACE_LAW.md) gives Λ = 7.4, about +4–6 in distance, and is not yet propagated.', '',
              f'Validation: {validation["status"]}. Reproduce with `make optimization-round` (compute node recommended; '
              '`pbs/optimization_round.pbs`), redraw with `make optimization-round-figures`, audit with `make optimization-round-check`.', '']
    (OUT/'summary.md').write_text('\n'.join(lines))


def main():
    rows, table = read(OUT/'costs.json'), read(OUT/'capacity.json')
    plot_family(rows, 'tau_L_over_4', 'optimization_round_capacity',
                'Time-growing workload (τ = L/4): capacity after the optimization round')
    plot_family(rows, 'tau_1', 'optimization_round_fixed_tau',
                'Fixed τ = 1 implementation benchmark: capacity after the optimization round')
    plot_leaders(table)
    summary(rows, table, read(OUT/'changes.json'), read(OUT/'validation.json'))


if __name__ == '__main__':
    main()
