"""Figures and a compact report from the saved published-benchmark replay."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from published_benchmarks import OUT, read
from project_paths import ROOT


def save(fig, name):
    for extension in ('pdf','png'):
        fig.savefig(ROOT/f'figures/{name}.{extension}',dpi=180,bbox_inches='tight',
                    metadata={'CreationDate':None,'ModDate':None} if extension=='pdf' else {})
    plt.close(fig)


def topqad_plot(data):
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.8))
    ax=axes[0]
    for hw,color in [('target','#bc702f'),('desired','#2079a5')]:
        rows=[v for v in data['topqad_star'] if v['hardware']==hw]
        ax.scatter([v['pec_free_seconds_published']/60 for v in rows],
                   [v['pec_free_seconds_replayed']/60 for v in rows],c=color,label=hw.capitalize(),s=36)
        for v in rows:
            if abs(v['pec_free_ratio']-1)>.05:
                ax.annotate(f"U={v['U']}, L={v['L']}",
                    (v['pec_free_seconds_published']/60,v['pec_free_seconds_replayed']/60),
                    xytext=(5,6),textcoords='offset points',fontsize=8)
    ax.plot([1,25],[1,25],color='k',ls='--',lw=1)
    ax.set(xscale='log',yscale='log',xlim=(2,25),ylim=(2,25),
           xlabel='Published STAR time without PEC [min]',ylabel='Replayed time without PEC [min]',
           title='Table IV: area and time accounting')
    ax.legend(frameon=False)
    ax=axes[1]
    v=next(v for v in data['topqad_star'] if (v['L'],v['U'],v['hardware'])==(8,4,'target'))
    points=[('STAR',v['qubits_published'],v['seconds_published'],'#bc702f')]
    for mode,label in [('space','FT: space optimal'),('time','FT: time optimal')]:
        v=next(v for v in data['topqad_ft'] if (v['L'],v['U'],v['hardware'],v['mode'])==(8,4,'target',mode))
        points.append((label,v['qubits_published'],v['seconds_published'],'#7452b4'))
    for label,Q,t,color in points:
        ax.scatter(Q,t/60,c=color,s=70)
        ax.annotate(label,(Q,t/60),xytext=(0,10),textcoords='offset points',ha='center',fontsize=9)
    ax.set(xscale='log',yscale='log',xlim=(1.5e5,1.3e6),ylim=(.2,20),
           xlabel='Physical qubits',ylabel='Published total runtime [min]',
           title='Same GSEE task: 8 × 8, U = 4, energy error 0.32')
    for ax in axes:ax.grid(alpha=.2)
    fig.text(.5,.01,'Replay uses published distances and RUS clocks; it does not validate injection scheduling.\n'
             'Right: original paper estimates, including its different STAR and FT algorithms; not our dynamics task.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.12,1,1));save(fig,'published_topqad')


def phasecraft_plot(data):
    fig,axes=plt.subplots(2,2,figsize=(11.4,7),sharex=True)
    allrows=data['observables']
    for col,U in enumerate((0,4)):
        rows=[v for v in allrows if v['U']==U and v['observable']=='mean_doublon']
        ax=axes[0,col]
        for method,color,style in [('Raw','#999999','o'),('TFLO','#d08a39','.'),('Mitigated','#177aa0','o-')]:
            selected=sorted([v for v in rows if v['method']==method],key=lambda v:v['tau'])
            ax.errorbar([v['tau'] for v in selected],[v['value'] for v in selected],
                        yerr=[v['reported_se'] or 0 for v in selected],fmt=style,color=color,
                        ms=3,lw=1,label=method if method!='Mitigated' else 'TFLO + smoothing')
        references=[('FLO','k','--')] if U==0 else [('TDVP chi=1024','#ce6994',':'),('TDVP chi=2048','k','--')]
        for method,color,style in references:
            selected=sorted([v for v in rows if v['method']==method],key=lambda v:v['tau'])
            ax.plot([v['tau'] for v in selected],[v['value'] for v in selected],style,color=color,lw=1.5,label=method)
        ax.set_title(f'U = {U}: '+('training consistency' if U==0 else 'comparison with finite-bond MPS'))
        ax.set_ylabel('Spatially averaged doublon density')
        ax.legend(frameon=False,fontsize=8,loc='lower right')
        ref_method='FLO' if U==0 else 'TDVP chi=2048'
        ref={round(v['tau'],8):v['value'] for v in rows if v['method']==ref_method}
        ax=axes[1,col]
        for method,color in [('Raw','#999999'),('TFLO','#d08a39'),('Mitigated','#177aa0')]:
            selected=sorted([v for v in rows if v['method']==method and v['tau']>0],key=lambda v:v['tau'])
            ax.plot([v['tau'] for v in selected],[v['value']-ref[round(v['tau'],8)] for v in selected],
                    '.-',color=color,lw=1,ms=4)
        ax.axhspan(-.01,.01,color='#555555',alpha=.1,label='±0.01 guide')
        ax.axhline(0,color='k',lw=.7)
        ax.set(xlabel='Evolution time Jt',ylabel='Difference from reference')
        ax.legend(frameon=False,fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.15)
    fig.suptitle('Published 56-qubit experiment: 160 shots per time point, with shared information',y=1.01)
    fig.text(.5,.01,'7 × 4 periodic lattice · dimer/defect initial state · four Trotter steps · Zenodo 17799843\n'
             'Error bars are published one-sigma uncertainties. U=4 differences include Trotter and reference errors; they are not residual-bias measurements.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.09,1,.97));save(fig,'published_phasecraft')


def report(result, data):
    star=result['topqad_star'];ft=result['topqad_ft']
    ordinary=[v for v in star if abs(v['pec_free_ratio']-1)<=.05]
    lines=['# Published benchmarks: what is reproduced and what remains assumed','',
      '[TopQAD figure](../../figures/published_topqad.pdf) · [Experimental figure](../../figures/published_phasecraft.pdf) · [Methods](../../docs/PUBLISHED_BENCHMARKS.md)','',
      'The published examples are compatible with the resource laws once task, estimator, parallelism and accuracy are kept explicit. '
      'This pass does not fit an experimental mitigation exponent or replace the integrated capacity curves.','',
      '## TopQAD / STAR: Table IV','',
      f'All {len(star)} STAR footprints reproduce the published two-significant-digit values from (4L²+1)2d². '
      f'{len(ordinary)} of 16 PEC-free runtimes agree within {100*max(abs(v["pec_free_ratio"]-1) for v in ordinary):.1f}%. '
      f'The {sum(v["mode"]=="time" for v in ft)} FT time-optimal entries agree within '
      f'{100*max(abs(v["runtime_ratio"]-1) for v in ft if v["mode"]=="time"):.1f}% with one T consumed per d syndrome cycles.','',
      'The replay supplies the paper’s distances, RUS clocks and T counts. It reconstructs the remaining area/time ledger, '
      'including both Hadamard-test quadratures, the multilevel time grid, integer steps and control clocks. '
      'The growth inferred from the total/PEC-free ratio is a diagnostic, not a forward validation of the growth fit.','',
      '| U | L | Hardware | PEC-free replay / published | Total time with anchor growth / published | Implied growth per state | Flag |',
      '|---:|---:|---|---:|---:|---:|---|']
    for v in star:
        g='no nonnegative solution' if v['implied_growth_per_state'] is None else f"{v['implied_growth_per_state']:.2g}"
        lines.append(f"| {v['U']} | {v['L']} | {v['hardware']} | {v['pec_free_ratio']:.3f} | {v['anchor_transfer_ratio']:.3f} | {g} | {v['flags']} |")
    lines += ['',
      'For target hardware, the inferred growth contribution clusters near 3×10⁻⁶. That numerical agreement is suggestive, '
      'but the inference uses the published runtime and cannot validate it independently. '
      'Fitting only U=4,L=4 and transferring the same growth parameter reproduces all seven other target-hardware total runtimes within 10%. '
      'This is a calibrated workload-transfer check, not independent gadget validation. Two desired-hardware rows remain inconsistent under the declared replay: '
      'U=8,L=8 lists 750 clocks, which predicts 470 ms rather than 297 ms; U=8,L=6 has a PEC factor below the angle-only estimate. '
      'Neither row is silently repaired. The preceding STAR supplementary-code endpoint returned HTTP 403 during this pass; '
      'a complete independent acceptance/growth/schedule reproduction remains open.','',
      'At L=8,U=4 under target hardware, STAR uses 230,000 qubits and 8.6 minutes. FT uses 570,000 qubits and 1.9 minutes '
      'in the space-optimal layout or 890,000 qubits and 21.2 seconds in the time-optimal layout. '
      'This is an area-versus-time tradeoff on the paper’s GSEE task, not an ordering for our local dynamics observable.','',
      '### Change one assumption on the fixed 4×4 GSEE workload','',
      '| Single change | Runtime [minutes] | Ratio to anchor |',
      '|---|---:|---:|']
    for v in result['fixed_workload_sensitivities']:
        lines.append(f"| {v['scenario']} | {v['seconds']/60:.2f} | {v['relative_to_anchor']:.3f} |")
    lines += ['',
      'These isolate model sensitivities. The substituted growth bound has a different noise model, preparation and growth protocol; '
      'it is not a calibrated migration of the published hardware. All other clocks, distances and schedule choices stay fixed.','',
      '## Phasecraft / Quantinuum: direct data reduction','',
      'The pinned archive contains 3,200 shots per interaction value, 160 shots at each of 20 times, on 56 qubits. '
      'The mean doublon reconstructed from raw shots agrees with every saved raw mean to floating-point precision. '
      'The following RMS differences use all 20 hardware times; local entries also use all 28 sites.','',
      '| U | Observable | Method | Reference | RMS difference | Median published standard error |',
      '|---:|---|---|---|---:|---:|']
    for v in result['phasecraft_comparisons']:
        if v['reference'] not in ('FLO','TDVP chi=2048') or v['group'] not in ('Raw','TFLO','Mitigated'):continue
        se='—' if v['median_reported_se'] is None else f"{v['median_reported_se']:.4f}"
        lines.append(f"| {v['U']} | {v['observable']} | {v['group']} | {v['reference']} | {v['rmse']:.4f} | {se} |")
    lines += ['',
      'For U=4 the raw mean-doublon variance per shot is about 0.0035, compared with about 0.149 for an average individual '
      'site and the discussion model’s Bernoulli reference 0.09. Spatial averaging and correlations already change the shot budget substantially; '
      'smoothing also shares information across time. Thus a single-point gatewise-PEC estimate is not a prediction of the experimental shot count. '
      'The U=0 comparison uses training information. U=4 differences combine experimental, Trotter and finite-bond reference errors, '
      'and cannot be interpreted as held-out residual bias or as proof of failure of mitigation.','',
      '## Confidence after this pass','',
      '| Item | Evidence now | Remaining assumption |',
      '|---|---|---|',
      '| NISQ experimental footprint and shot count | Direct archive census and exact reconstruction of raw means | Different state, geometry and observable from our baseline |',
      '| NISQ useful observable accuracy | Published processed curves and reference comparisons | U=4 exact truth and transfer to a chosen local target |',
      '| STAR published area and PEC-free ledger | 16 areas and 14 consistent runtimes | Supplied RUS clocks and distances; full noisy injection/growth |',
      '| STAR mitigation overhead | Inferred growth scale from published PEC ratio | Independent forward channel prediction |',
      '| FT time-optimal throughput | 16 one-port timing checks | Published synthesis counts and private factory/layout optimizer |','',
      'Willow shot data are also public: [Zenodo 17881822](https://zenodo.org/records/17881822), an 8.8 GB archive with raw shots '
      'and processed expectations. Its README and metadata were checked; that archive was not reduced here. '
      'This corrects the earlier uncertainty about availability and makes a hardware-matched follow-up possible.','']
    (OUT/'summary.md').write_text('\n'.join(lines))


def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    result=read(OUT/'benchmarks.json');data=read(OUT/'phasecraft_reduced.json')
    topqad_plot(result);phasecraft_plot(data);report(result,data)


if __name__=='__main__':main()
