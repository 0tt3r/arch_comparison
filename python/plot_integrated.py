"""Current integrated capacity and NISQ requirements, from saved calculations."""
from math import log10

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from discussion_model import capacity
from integrated_comparison import ARMS, OUT, read, settings
from plot_capacity import COLORS
from project_paths import ROOT

LABELS = {'star':'STAR', 'hybrid':'Hybrid', 'ft':'FTQC'}


def save(fig, name):
    for suffix in ('pdf','png'):
        fig.savefig(ROOT/f'figures/{name}.{suffix}', dpi=180, bbox_inches='tight',
                    metadata={'CreationDate':None,'ModDate':None} if suffix=='pdf' else {})
    plt.close(fig)


def plot_capacity(rows):
    cfg = settings()
    previous = read(ROOT/'results/physics_sweep_v1/costs.json')
    fig,axes = plt.subplots(2,2,figsize=(11.7,8),sharex=True,sharey=True)
    for i,eps in enumerate(cfg['precisions']):
      for j,p in enumerate(cfg['physical_errors']):
        ax = axes[i,j]
        group = [v for v in rows if v['p']==p and v['absolute_error']==eps and v['tau']==1]
        thresholds = [10**v['log10_Q_min'] for v in group
                      if v['log10_Q_min'] is not None and 2<=v['log10_Q_min']<=9]
        Q = np.unique(np.r_[np.geomspace(1e2,1e9,350), thresholds])
        for arm in ARMS:
            sides = np.array([[capacity([v for v in group if v['arm']==arm and v['scenario']==scenario],q)
                               for q in Q] for scenario in cfg['scenarios']])
            ax.fill_between(Q,sides.min(axis=0),sides.max(axis=0),step='post',color=COLORS[arm],alpha=.14)
            central = [v for v in group if v['arm']==arm and v['scenario']=='component']
            ax.step(Q,[capacity(central,q) for q in Q],where='post',color=COLORS[arm],lw=2)
        for scenario,color,style in [('effective',COLORS['nisq'],'--'),('pec_two_qubit','#8c622c',':')]:
            old = [v for v in previous if (v['p'],v['absolute_error'],v['arm'],v['scenario'])==(p,eps,'nisq',scenario)]
            ax.step(Q,[capacity(old,q) for q in Q],where='post',color=color,ls=style,lw=1.7)
        old = [v for v in previous if (v['p'],v['absolute_error'],v['arm'],v['scenario'])==(p,eps,'star','held')]
        ax.step(Q,[capacity(old,q) for q in Q],where='post',color='#777777',ls=(0,(1,3)),lw=1)
        ax.set(xscale='log',xlim=(1e2,1e9),ylim=(0,49),yticks=[0,4,8,16,24,32,40,48],
               title=rf'$p={p:g},\quad \epsilon={eps:g}$')
        ax.grid(alpha=.2)
        if j==0: ax.set_ylabel('Largest modeled lattice side L')
        if i==1: ax.set_xlabel('Available physical qubits Q')
    handles = [Line2D([],[],color=COLORS[a],lw=2,label=LABELS[a]) for a in ARMS]
    handles += [Line2D([],[],color=COLORS['nisq'],ls='--',label='Physical: assumed effective mitigation'),
                Line2D([],[],color='#8c622c',ls=':',label='Physical: two-qubit PEC reference'),
                Line2D([],[],color='#777777',ls=(0,(1,3)),label='Previous STAR assumptions')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.065),ncol=3,frameon=False,fontsize=9)
    fig.suptitle('Updated component scenarios: the encoded ordering remains conditional',y=.99)
    fig.text(.5,.015,'τ = 1 · one week · shared finite-patch accuracy model\n'
             'Bands span declared component/supply scenarios, not confidence intervals. Physical-qubit effective bias remains assumed.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.17,1,.96))
    save(fig,'integrated_capacity')


def plot_requirements(rows):
    fig,axes = plt.subplots(1,2,figsize=(11.7,4.8),sharey=True)
    colors = {4:'#2079a5',6:'#bc702f',8:'#8550a4'}
    for ax,eps in zip(axes,(.01,.001)):
        for L,color in colors.items():
            selected = [v for v in rows if (v['kind'],v['p'],v['tau'],v['epsilon'],v['L'],v['Q'],v['training_attempts'])
                        == ('F_max',.001,1,eps,L,1000000,0) and v['F_max'] is not None]
            for speedup,style in [(1,'-'),(10,'--')]:
                central = sorted([v for v in selected if v['scenario']=='component' and v['speedup']==speedup],key=lambda v:v['bias_fraction'])
                ax.plot([v['F_max'] for v in central],[v['bias_fraction'] for v in central],style,
                        color=color,lw=2,label=f'{L} × {L}' if speedup==1 else None)
                if speedup==1:
                    fractions = sorted(set(v['bias_fraction'] for v in central))
                    groups = [[v['F_max'] for v in selected if v['speedup']==1 and v['bias_fraction']==f] for f in fractions]
                    ax.fill_betweenx(fractions,[min(g) for g in groups],[max(g) for g in groups],color=color,alpha=.12)
        ax.set(xscale='log',ylim=(0,1),xlabel='Largest affordable variance multiplier F',
               title=rf'$\epsilon={eps:g}$')
        ax.grid(alpha=.2)
    axes[0].set_ylabel('Allowed remaining offset / total error target')
    axes[0].legend(frameon=False,loc='lower left')
    fig.suptitle('What must physical-qubit mitigation achieve to beat the selected encoded implementations?',y=.99,fontsize=12)
    fig.text(.5,.04,'Q = 10⁶ · p = 10⁻³ · τ = 1 · Gaussian planning at D* = 0.1 · no training charge in this panel\n'
             'Solid: match runtime. Dashed: 10× faster. Both also meet one week. Affordable region lies left of each curve.\n'
             'Steps are reselected along each curve; F must be achievable at that step count. Training sensitivity is in the tables.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.2,1,.94))
    save(fig,'nisq_requirements')


def report(rows, thresholds):
    cfg = settings()
    lines = ['# Integrated comparison and NISQ accuracy requirements','',
        'This is the current integrated comparison. Historical sweeps remain available for replay. '
        'The encoded bands are scenario envelopes; they do not establish statistical error bars or an unconditional architecture ordering.','',
        '[Capacity figure](../../figures/integrated_capacity.pdf) · [NISQ requirements](../../figures/nisq_requirements.pdf) · [Methods](../../docs/INTEGRATED_COMPARISON.md)','',
        '## Capacity at tau=1 and one million qubits','',
        'Entries are central lattice sides with the full scenario range in brackets. The shared accuracy-transfer and delivery assumptions still apply.','',
        '| p | Error | STAR | Hybrid | FTQC | Unique leader across scenarios |',
        '|---:|---:|---:|---:|---:|---|']
    for v in read(OUT/'capacity.json'):
        if v['tau'] != 1 or v['Q'] != 1000000:
            continue
        values = [f"{v[a+'_component']} [{v[a+'_min']}, {v[a+'_max']}]" for a in ARMS]
        lines.append(f"| {v['p']:g} | {v['epsilon']:g} | " + ' | '.join(values)
                     + f" | {v['robust_capacity_leader']} |")
    lines += ['',
        'The table leaves close rankings unresolved. The full capacity table also identifies regions with a unique leader across this menu; '
        'that statement remains conditional on the shared component-composition and layout assumptions.','',
        '## Current 4×4, tau=1 operating points','',
        '| p | Error | Arm | Scenario | Q for one week | Q per copy | seconds/attempt | steps | injection seed / data d | log variance multiplier |',
        '|---:|---:|---|---|---:|---:|---:|---:|---|---:|']
    for v in sorted(rows, key=lambda v:(cfg['physical_errors'].index(v['p']),
                      cfg['precisions'].index(v['absolute_error']),
                      list(cfg['scenarios']).index(v['scenario']), ARMS.index(v['arm']))):
        if v['L']!=4 or v['tau']!=1 or v['log10_Q_min'] is None:
            continue
        lines.append(f"| {v['p']:g} | {v['absolute_error']:g} | {v['arm']} | {v['scenario']} | "
                     f"{10**v['log10_Q_min']:,.0f} | {v['qubits_per_copy']:,} | {v['seconds_per_attempt']:.4g} | "
                     f"{v['r']} | {v.get('injection_distance','—')} / {v['distance']} | {v['chi']:.4g} |")
    lines += ['', 'STAR and hybrid use the exact existing generator-fault submodel at the workload angles, '
              'including noisy PCEC, repeat-until-success service time, detector rejection, coherent residue and retry tails. '
              'Seed 5 uses its measured growth-campaign upper bound; no claim is made about untested seed 7. '
              'The zero-growth sensitivity is available for seed 9, whose zero-failure result is an upper limit rather than a zero-error measurement. '
              'Seed 25 avoids growth but pays its preparation acceptance. Receivers use distance-13 or distance-25 component anchors. '
              'Greater data distances add protected handoff time and volume; their injection channels are not newly calibrated.', '',
              'FT selects individual synthesis or Hamming-weight batches 4, 8, 16, and compact or fast one-port blocks. '
              'These are checked arithmetic counts with an aggregate delivery model, not full routed schedules. '
              'Hybrid fixed-angle supply uses the same cultivation quality and source menu. '
              'The slow scenario increases detector rejection and multiplies source preparation cycles by four.', '',
              '## NISQ requirements at Q=10^6, p=10^-3, tau=1', '',
              'The remaining offset is set to 10% of epsilon in this table. Training sensitivity charges six million additional '
              'attempts of the same duration; it is not a proposed experimental training budget. '
              'Comparison is to the fastest of the selected encoded implementations that fits Q, under the component scenario. '
              'Those implementations were optimized for one-week area, not reoptimized for runtime at Q.', '',
              '| Error | L | Required speedup | Training shots | Encoded comparator | Largest F | Steps | kappa threshold at that step count |',
              '|---:|---:|---:|---:|---|---:|---:|---:|']
    for v in thresholds:
        if (v['kind'],v['p'],v['tau'],v['Q'],v['scenario'],v.get('bias_fraction')) != ('F_max',.001,1,1000000,'component',.1):
            continue
        if v['F_max'] is None:
            lines.append(f"| {v['epsilon']:g} | {v['L']} | {v['speedup']} | {v['training_attempts']} | {v['reference_arm']} | none | — | — |")
        else:
            lines.append(f"| {v['epsilon']:g} | {v['L']} | {v['speedup']} | {v['training_attempts']} | {v['reference_arm']} | "
                         f"{v['F_max']:.3g} | {v['r']} | {v['kappa_max']:.3g} |")
    lines += ['', 'F includes postselection through variance per attempted circuit. Unknown residual bias is an axis, not a fitted constant. '
              'Smoothing, training and other information-sharing strategies must meet the displayed variance and offset requirements '
              'with their own resource accounting. A larger affordable F at greater size or stricter accuracy can arise from a slower encoded comparator; '
              'it does not establish easier mitigation at the deeper physical circuit. '
              'No large noisy-dynamics calculation or experimental pipeline fit was performed.', '',
              '## Interpretation', '',
              'A capacity leader is marked robust only when its smallest capacity across the declared scenarios exceeds every '
              'other encoded arm’s largest capacity. This is robustness to this menu, not a theorem about architectures. '
              'Ties and overlaps remain unresolved. NISQ has no unconditional size boundary in these results. '
              'The data-distance floor, transferred finite-lattice Trotter coefficients, growth composition, source quality and '
              'uncompiled delivery remain named assumptions.', '',
              'Run `make integrated` to resume/recompute the new version and regenerate its figures. '
              '`make integrated-check` audits source contracts, component composition and accounting.','']
    (OUT/'summary.md').write_text('\n'.join(lines))


def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    rows,thresholds = read(OUT/'costs.json'),read(OUT/'break_even.json')
    plot_capacity(rows)
    plot_requirements(thresholds)
    report(rows,thresholds)


if __name__=='__main__':
    main()
