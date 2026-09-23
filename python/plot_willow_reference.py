"""Observed experimental discrepancy beside independently calculated circuit error."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from project_paths import ROOT
from willow_archive import read
from willow_reference import OUT


def run():
    data=read(OUT/'results.json');validation=read(OUT/'validation.json')
    colors={'experiment_vs_continuous_rms':'#263a46','experiment_vs_circuit_rms':'#167da0',
            'circuit_vs_continuous_rms':'#c17a29'}
    labels={'experiment_vs_continuous_rms':'Processed − continuous',
            'experiment_vs_circuit_rms':'Processed − ideal circuit',
            'circuit_vs_continuous_rms':'Ideal circuit − continuous'}
    fig,grid=plt.subplots(2,2,figsize=(10.5,7.2));axes=grid.ravel()
    for name,color,style in [('continuous','#263a46','-'),('ideal circuit','#c17a29','--'),
                             ('mitigated','#167da0','o'),('TDVP chi=2048','#a5a5a5',':')]:
        rows=sorted([v for v in data['curves'] if v['steps']==3 and v['method']==name],key=lambda v:v['tau'])
        for ax,observable in zip(axes[:2],('mean_d','central_d')):
            ax.plot([v['tau'] for v in rows],[v[observable] for v in rows],style,color=color,lw=1.5,ms=3,
                    label='Published processed' if name=='mitigated' else name.capitalize())
    axes[0].set(title='Spatial mean doublon, three steps',xlabel='Evolution time Jt',ylabel='Doublon density')
    axes[0].legend(frameon=False,fontsize=8)
    axes[1].set(title='Central-hole site doublon, three steps',xlabel='Evolution time Jt',ylabel='Local doublon')
    times=sorted([v for v in data['decomposition'] if v['steps']==3 and v['method']=='mitigated'],key=lambda v:v['tau'])
    steps=sorted([v for v in data['metrics'] if v['observable']=='local_d' and v['method']=='mitigated'],key=lambda v:v['steps'])
    for key in colors:
        axes[2].plot([v['tau'] for v in times],[v[key] for v in times],color=colors[key],label=labels[key])
        axes[3].plot([v['steps'] for v in steps],[v[key] for v in steps],'o-',color=colors[key],label=labels[key])
    axes[2].set(title='Local doublons, three steps',xlabel='Evolution time Jt',ylabel='RMS across sites',yscale='log')
    axes[3].set(title='All 24 times and 16 sites',xlabel='Trotter steps',ylabel='Local-doublon RMS difference',yscale='log',xticks=[1,2,3])
    axes[3].legend(frameon=False,fontsize=8)
    for ax in axes[2:]:ax.axhline(.01,color='#888',ls=':',lw=.8)
    for ax in axes:ax.grid(alpha=.15)
    fig.suptitle('Matched Willow reference: 4 × 4, U = 4, zero flux, central-hole initial state',fontsize=12)
    fig.text(.5,.012,'Circuit references include the experiment’s diagonal symmetry averaging. Experimental differences retain sampling and processing error.\n'
             'The three RMS curves do not add: the signed discrepancies can reinforce or cancel. No 6 × 6 extrapolation.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.07,1,.96))
    for ext in ('pdf','png'):
        fig.savefig(ROOT/f'figures/willow_matched_reference.{ext}',dpi=180,bbox_inches='tight',
                    metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else {})
    plt.close(fig)
    text=['# Willow with a matched interacting reference','',
      '[Figure](../../figures/willow_matched_reference.pdf) · [Methods](../../docs/WILLOW_REFERENCE.md)','',
      'We independently compute continuous dynamics and ideal experimental circuits for the 4×4, U=4, zero-flux center-hole input. '
      'All 24 nonzero saved times and all three step counts are included. This adds the missing interacting reference to the preceding archive analysis; '
      'it does not refit the published mitigation pipeline or extrapolate its performance.','',
      '## What the circuit reconstruction establishes','',
      f'The Appendix D2 swap-network order reproduces all {validation["free_circuit_reference_points"]} archived U=0 circuit-reference points '
      f'at both 4×4 and 6×6 and both fluxes, with maximum local-Z discrepancy {validation["maximum_free_circuit_Z_difference"]:.2g}. '
      'This requires diagonal symmetry averaging: the unsymmetrized Z discrepancy reaches 0.345. '
      'The onsite opposite-spin ZZ-derived doublons reproduce within 7.5e-9, the precision of the archived quantities. '
      'The U>0 onsite evolution and complete small-system state were separately checked against a full Fock-space matrix construction.','',
      '## Observed differences for the processed experimental estimates','',
      'All values below are local-doublon RMS differences across 24 times and 16 sites. The experimental discrepancy contains finite sampling and processing error; '
      'it is not an estimate of ensemble residual bias.','',
      '| Steps | Processed vs ideal circuit | Ideal circuit vs continuous | Processed vs continuous | Fraction of local points within 0.01 of continuous |',
      '|---:|---:|---:|---:|---:|']
    for v in steps:
        text.append(f'| {v["steps"]} | {v["experiment_vs_circuit_rms"]:.6f} | {v["circuit_vs_continuous_rms"]:.6f} | '
                    f'{v["experiment_vs_continuous_rms"]:.6f} | {v["fraction_total_within_01"]:.3f} |')
    text+=['',
      'The two contributions respond oppositely to refinement: the ideal dynamics improves from one to three steps, '
      'while the processed-data discrepancy from the ideal circuit grows. At three steps, the total RMS is smaller '
      'than the deterministic circuit error alone because the signed contributions partially cancel. This is a measured '
      'depth-versus-accuracy tradeoff for this pipeline, not a projection beyond three steps.','',
      'The signed identity is','',
      '    observed − continuous = (observed − ideal circuit) + (ideal circuit − continuous).','',
      'Consequently, the squared RMS total includes twice the mean product of the two signed contributions. '
      'The saved decomposition retains this cross term, so accidental cancellation is not credited as a generally reliable error budget.','',
      '## Approximate classical references','',
      '| MPS bond dimension | Local doublon RMS vs continuous | Maximum local difference |',
      '|---:|---:|---:|']
    for v in data['mps_comparison']:
        if v['steps']==3 and v['observable']=='local_d':
            text.append(f'| {v["method"].split("=")[-1]} | {v["continuous_reference_rms"]:.6g} | {v["continuous_reference_max"]:.6g} |')
    text+=['','## Numerical scope','',
      f'The (N_up,N_down)=(7,8) sector has 147,232,800 amplitudes. The {validation["reference_records"]} records comprise '
      '24 continuous time points, 72 ideal circuit points and one independent direct-time check. '
      f'Maximum norm error is {validation["maximum_norm_error"]:.2g}; the accumulated analytical Chebyshev state-tail bound is '
      f'{validation["maximum_accumulated_state_tail"]:.2g}. The independent direct evolution to Jt=1 agrees with incremental evolution '
      f'within {validation["direct_vs_incremental_max_observable_difference"]:.2g} in the reported observables. '
      'The tail bound excludes floating-point rounding; independent small-state comparisons and this direct-time check control the numerical implementation.','',
      'This finite-system calculation addresses the actual experimental input and ordering, which differ from the half-filled Néel/plaquette resource baseline. '
      'It does not certify every point at 0.01 or 0.001, determine a universal mitigation exponent, or validate larger interacting lattices.','']
    (OUT/'summary.md').write_text('\n'.join(text))


if __name__=='__main__':run()
