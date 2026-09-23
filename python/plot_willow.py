"""Experimental operating points and randomized-circuit sampling sensitivities."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from project_paths import ROOT
from willow_archive import OUT,read
from willow_benchmark import load_reduced,design_forecast


def save(fig,name):
    for suffix in ('pdf','png'):
        fig.savefig(ROOT/f'figures/{name}.{suffix}',dpi=180,bbox_inches='tight',
                    metadata={'CreationDate':None,'ModDate':None} if suffix=='pdf' else {})
    plt.close(fig)


def operating_points(result,data):
    fig,axes=plt.subplots(2,2,figsize=(10.7,7.4))
    for col,L in enumerate((4,6)):
        ax=axes[0,col]
        curves=[v for v in result['curves'] if v['L']==L and v['U']==4 and v['flux']=='0'
                and v['steps']==3 and v['observable']=='mean_d']
        case=next(v for v in data['cases'] if (v['L'],v['U'],v['flux'],v['steps'])==(L,4,'0',3))
        ax.plot([s['tau'] for s in case['shots']],
                [s['statistics']['unfiltered']['mean'][L*L] for s in case['shots']],
                ':',color='#aaa',label='Unfiltered shots')
        for method,color,style,label in [('raw','#888','--','Archive raw (PS + symmetry)'),
             ('mitigated','#167da0','o-','Published mitigation'),
             ('TDVP chi=1024','#d183a0',':','MPS 1024'),('TDVP chi=2048','black','--','MPS 2048')]:
            rows=sorted([v for v in curves if v['method']==method],key=lambda v:v['tau'])
            ax.plot([v['tau'] for v in rows],[v['value'] for v in rows],style,color=color,lw=1.4,ms=3,label=label)
        ax.set(xlabel='Evolution time Jt',ylabel='Spatially averaged doublon',ylim=(-.005,.27),
               title=f'{L} × {L}, U = 4, zero flux, 3 steps')
        if col==0:ax.legend(frameon=False,fontsize=8,loc='lower right')
        ax=axes[1,col]
        for U,color in [(0,'#188895'),(4,'#cb7d2b'),(8,'#874ba1')]:
            reference='continuous FLO' if U==0 else 'TDVP chi=2048'
            for method,style in [('raw',':'),('mitigated','o-')]:
                rows=sorted([v for v in result['comparisons'] if v['L']==L and v['U']==U and v['flux']=='0'
                    and v['observable']=='local_d' and v['reference']==reference and v['method']==method],key=lambda v:v['steps'])
                ax.plot([v['steps'] for v in rows],[v['rmse'] for v in rows],style,color=color,
                        ms=4,label=f'U = {U}' if method=='mitigated' else None)
        ax.axhline(.01,color='k',lw=.8,alpha=.5)
        ax.set(xlabel='Trotter steps',ylabel='Local-doublon RMS reference difference',yscale='log',xticks=[1,2,3],ylim=(.003,.25))
        if col==0:ax.legend(frameon=False,fontsize=9)
        ax.set_title('Solid: processed estimates · dotted: archive raw',fontsize=10)
    for ax in axes.flat:ax.grid(alpha=.16)
    fig.suptitle('Willow operating points: 200,000 attempted shots per time, 20 randomized circuits',fontsize=12)
    fig.text(.5,.005,'Lower panels use all 24 times and every site. U=0: exact continuous FLO, but mitigation used training data.\n'
             'U>0: finite-bond MPS reference agreement, not certified accuracy. Derived-observable error covariances are not supplied.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.07,1,.96));save(fig,'willow_operating_points')


def sampling_design(result,data):
    fig,axes=plt.subplots(1,2,figsize=(10.7,4.6))
    for L,color in [(4,'#177da4'),(6,'#b77832')]:
        ax=axes[0];med=[];low=[];high=[]
        for r in (1,2,3):
            values=[s['statistics']['strict']['acceptance'] for c in data['cases']
                    if c['L']==L and c['steps']==r for s in c['shots']]
            med.append(np.median(values));low.append(min(values));high.append(max(values))
        ax.plot([1,2,3],med,'o-',color=color,label=f'{L} × {L}')
        ax.fill_between([1,2,3],low,high,color=color,alpha=.12)
        ax=axes[1]
        selection='strict' if L==4 else 'relaxed2'
        v=next(v for v in result['sampling'] if v['L']==L and v['filter']==selection and v['observable']=='central_d')
        s=next(s for c in data['cases'] if (c['L'],c['U'],c['flux'],c['steps'])==(L,4,'0',3)
               for s in c['shots'] if s['tau']==1.)['statistics'][selection]
        i=(L//2-1)*L+(L//2-1);b=s['between_mean_variance'][i];w=s['within_mean_variance'][i]
        h=np.logspace(0,2,100)
        ax.loglog(h,design_forecast(b,w,20,shot_multiplier=h),color=color,lw=1.6)
        ax.loglog(h,design_forecast(b,w,20,twirl_multiplier=h),'--',color=color,lw=1.6)
        ax.scatter([1],[v['cluster_se']],color=color,zorder=4,s=22)
    axes[0].set(xlabel='Trotter steps',ylabel='Strict number-sector acceptance',xticks=[1,2,3],yscale='log',
                title='Postselection pays for growing circuit depth')
    axes[0].legend(frameon=False)
    axes[1].set(xlabel='Attempted-shot budget / 200,000',ylabel='Central-doublon standard error',
                title='Planning: repeat shots or redraw circuits?')
    axes[1].legend(handles=[Line2D([0],[0],color='k',label='Same 20 randomized circuits'),
                           Line2D([0],[0],color='k',ls='--',label='More independent randomized circuits')],
                   frameon=False,fontsize=8,loc='lower left')
    for ax in axes:ax.grid(alpha=.18)
    fig.text(.5,.015,'Left: median and range over interactions, fluxes and times. Right: U=4, zero flux, 3 steps, Jt=1, published number filters.\n'
             'Right-hand curves extrapolate 20 observed circuit groups; they omit training, mitigation and remaining systematic error.',
             ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.11,1,1));save(fig,'willow_sampling_design')


def report(result,validation):
    lines=['# Willow: measured operating points and what the sampling budget means','',
      '[Operating-point figure](../../figures/willow_operating_points.pdf) · [Sampling figure](../../figures/willow_sampling_design.pdf) · [Methods](../../docs/WILLOW_BENCHMARK.md)','',
      'This pass compares the public experiment with the questions in our resource model. It does not fit a universal mitigation exponent or change the integrated capacity curves.','',
      f'The reduction covers {validation["cases"]} declared center-hole workloads: 4×4 and 6×6, both fluxes, U=0,4,8, and one to three steps. '
      f'It includes {validation["shot_batches"]} time points and {validation["attempted_shots"]:,} attempted shots. '
      'Each time point has 20 randomized compilations with 10,000 shots each. A U=0 training time series costs additional shots, and the published smoothing uses shared time-series information.','',
      '## Resource scale','',
      '| L | Steps | Physical qubits | Two-qubit gates | Strict acceptance, median | Hypothetical two-qubit PEC variance factor, p=0.001 |',
      '|---:|---:|---:|---:|---:|---:|']
    for v in result['resources']:
        lines.append(f'| {v["L"]} | {v["steps"]} | {v["physical_qubits"]} | {v["gates"]} | {v["strict_acceptance_median"]:.3f} | {v["gatewise_2q_PEC_variance_at_p001"]:.3g} |')
    lines += ['',
      'The paper counts 34 and 74 physical qubits, including ancillas; the models have 32 and 72 fermionic modes. '
      'Each 24-time series uses 4.8 million shots. The PEC column is a specified inverse-channel calculation at our reference p, not a prediction of the experimental pipeline. '
      'The experiment uses different hardware noise and a trained estimator, and its useful accuracy must be assessed separately.','',
      '## Accuracy scale, without promoting approximate references to truth','',
      'The table below uses three steps and zero flux, all 24 nonzero times and all sites. Derived doublon uncertainties are covariance-agnostic upper bounds from the supplied Z/ZZ errors; they are not measured total errors.','',
      '| L | U | Reference | Local doublon RMS difference | Mean doublon RMS difference | Local Z median reported SE |',
      '|---:|---:|---|---:|---:|---:|']
    for L in (4,6):
        for U in (0,4,8):
            ref='continuous FLO' if U==0 else 'TDVP chi=2048'
            rows={v['observable']:v for v in result['comparisons'] if v['L']==L and v['U']==U
                  and v['flux']=='0' and v['steps']==3 and v['method']=='mitigated' and v['reference']==ref}
            lines.append(f'| {L} | {U} | {ref} | {rows["local_d"]["rmse"]:.4f} | {rows["mean_d"]["rmse"]:.4f} | {rows["local_z"]["median_reported_se_or_bound"]:.4f} |')
    lines += ['',
      'U=0 is a training consistency check. The archived `exact` values change with Trotter depth and describe the free circuit. '
      'We independently compute continuous U=0 dynamics by orbital evolution, checked against the existing many-body sector solver. '
      'At U=4 and 8, the exact groups are empty: the exact interacting comparisons shown in the paper are not supplied in this pinned archive. '
      'The MPS comparisons therefore include reference and Trotter errors, not just experimental error.','',
      '## Sampling design matters before an exponential cost law','',
      'Shots from one randomized circuit share its noise realization. With K circuit groups, n retained shots per group and conditional observable variance v, '
      'a useful planning decomposition is Var(mean) ≈ (a + v/n)/K, where a is variation between randomized-circuit expectations. '
      'Increasing n at fixed K does not remove the a/K contribution; redrawing independent circuits can. '
      'The fitted nonnegative a is a diagnostic from only 20 groups, not a certified irreducible floor.','',
      '| L | Central doublon SE now | Same circuits, 100× shots | 100× independent circuits |',
      '|---:|---:|---:|---:|']
    for v in result['sampling']:
        if v['filter']==('strict' if v['L']==4 else 'relaxed2') and v['observable']=='central_d':
            lines.append(f'| {v["L"]} | {v["cluster_se"]:.4f} | {v["same_twirls_shots_x100_se"]:.4f} | {v["new_twirls_x100_se"]:.5f} |')
    lines += ['',
      'These examples use U=4, zero flux, three steps and Jt=1, before TFLO and smoothing. They account for number-selection rejection '
      '(exact sector at 4×4; at most two total occupation-number deviations at 6×6). '
      'They forecast sampling variance only, assuming new circuit groups have the same distribution and are independent; they do not forecast mitigation accuracy.','',
      '## Archive and pipeline checks','',
      f'Raw published Z and onsite ZZ-derived doublons reproduce from equal-weight circuit groups, number selection and diagonal symmetry in '
      f'{validation["raw_pipeline_replayed"]} of 36 cases using the declared strict/relaxed/no-selection menu. '
      'The recovered number-deviation thresholds are zero at 4×4, two at 6×6 with zero flux, and one at 6×6 with pi flux. '
      'Thus the archive label `raw` already includes processing; it must not be compared to unfiltered bits as though they were the same estimator.','',
      'The initial-state JSON repeats a singlet configuration under the center-hole labels, including out-of-range coordinates for 4×4. '
      'It is retained and flagged, not used. The consistent product-state occupations are read from the zero-time TDVP observables and checked for one missing particle. '
      'All selected ZIP members have checked sizes/CRCs and SHA-256 hashes. About 1.06 GB was fetched by byte ranges; the full 8.8 GB archive MD5 was not verified.','',
      '## Consequence for the architecture comparison','',
      'The physical-qubit experiment is a real low-footprint operating point with a shared, trained estimator and observable-dependent accuracy. '
      'Our gate-count scale can describe its circuit size, but a gatewise PEC budget cannot describe its estimator. '
      'Neither the experimental curves nor their reference differences certify every local observable at 0.01 or 0.001. '
      'The useful next accuracy check is the missing matched interacting circuit/continuous reference for this pipeline, not a larger arbitrary held-out training campaign. '
      'For STAR/FT, full injection/growth and delivery remain separate unresolved physical inputs.','']
    (OUT/'summary.md').write_text('\n'.join(lines))


def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    result=read(OUT/'benchmarks.json');data=load_reduced()
    operating_points(result,data);sampling_design(result,data);report(result,read(OUT/'validation.json'))


if __name__=='__main__':main()
