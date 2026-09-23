"""Replay accounting, plot new evidence, and write the frontier study report."""
import argparse
import gzip
import hashlib
import json
from math import exp,log10

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from geometric_delivery import validate_events
from project_paths import ROOT

OUT=ROOT/'results/frontier_v1'


def read(name):return json.loads((OUT/name).read_text())

def save(fig,name):
    for suffix in ('pdf','png'):fig.savefig(ROOT/f'figures/{name}.{suffix}',bbox_inches='tight',dpi=180)
    plt.close(fig)


def validate():
    for name in ('qubitization/results.json','delivery/results.json','exposure.json','locality_certificate.json'):
        data=read(name)
        for path,digest in data['source_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    q=read('qubitization/results.json');assert len(q['rows'])==128 and len(q['inventories'])==8
    from discussion_model import assumptions
    from qubitization_study import cost_candidate
    base={**assumptions(),'data_block':'fast','ft_clocks_per_T':1,'T_consumption_ports':1}
    for v in q['rows']:
        assert abs(sum(v['budget_fractions'].values())-1)<1e-12
        if v['variant'].startswith('trotter'):continue
        inv=next(i for i in q['inventories'] if (i['L'],i['variant'])==(v['L'],v['variant']))
        b=v['budget_fractions'];replay=cost_candidate(inv,v,v['epsilon'],{**base,'qsp_clifford_serial_fraction':v['clifford_serial_fraction']},b['synthesis'],b['logical'],v['cleanup_levels'])
        for k in ('T_count','distance','factories','log10_Q_min','seconds_per_attempt'):assert abs(replay[k]-v[k])<1e-10,(k,v)
    delivery=read('delivery/results.json')
    for v in delivery['rows']:
        with gzip.open(OUT/f"delivery/batch{v['batch']}_ports{v['ports']}.json.gz",'rt') as f:events=json.load(f)
        validate_events(events)
        assert sum(x['kind']=='T_delivery' for x in events)==v['T_count']
    from exposure_diagnostics import growth_summary
    growth_summary()
    growth=read('growth_summary.json');assert all(v['simultaneous_upper95']>0 for v in growth['cases'])
    return dict(status='PASS: same-task budgets and source-menu replay, geometric reservations, saved channel data, growth census and PBS F/0',
        algorithm_costs=128,encoding_inventories=8,delivery_schedules=8,held_out_estimators=384,growth_cases=24,growth_shots=24000000,
        scope='Calculation and component validation; not a full physical injection or routed qubitization implementation')


def algorithms():
    rows=read('qubitization/results.json')['rows'];fig,axes=plt.subplots(1,2,figsize=(10.5,4),layout='constrained')
    for ax,epsilon in zip(axes,(.01,.001)):
        for variant,serial,label,color,style in [('trotter_individual',None,'Trotter: individual rotations','#555555','-'),
            ('trotter_hwp',None,'Trotter: Hamming-weight scenario','#178c70','-'),
            ('split_uniform',0.,'QSVT: T-service timing','#784bb1','-'),
            ('split_uniform',1.,'QSVT: serial Clifford sensitivity','#784bb1',':')]:
            selected=[v for v in rows if v['epsilon']==epsilon and v['tau']==1 and v['clifford_serial_fraction']==serial and
                      (v['variant']==variant if variant.startswith('trotter') else not v['variant'].startswith('trotter'))]
            budgets=np.unique(np.r_[np.geomspace(1e4,1e10,400),[10**v['log10_Q_min'] for v in selected]])
            capacity=[max([v['L'] for v in selected if v['log10_Q_min']<=log10(Q)+1e-12],default=0) for Q in budgets]
            ax.step(budgets,capacity,where='post',label=label,color=color,ls=style,lw=2)
        ax.set(xscale='log',xlim=(1e4,1e10),ylim=(0,8.6),yticks=[0,2,4,6,8],xlabel='Available physical qubits Q',ylabel='Largest tested lattice side L',title=rf'$\tau=1,\ \epsilon={epsilon:g}$')
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8,loc='upper left');save(fig,'frontier_algorithms')


def exposures():
    data=read('exposure.json');fig,axes=plt.subplots(1,3,figsize=(13,3.8),layout='constrained')
    names=['raw','number','tflo_raw','tflo_number'];colors=['#777777','#d18b00','#2463a8','#008c79']
    for method,color,label in zip(names,colors,['Raw','Number selection','Trained affine','Selection + training']):
        rows=[v for v in data['nisq'] if v['method']==method and v['noise']=='matched']
        for passed,marker in ((True,'o'),(False,'x')):
            points=[v for v in rows if v['passes_bias_allowance']==passed]
            axes[0].scatter([v['exposure'] for v in points],[v['log_measured_variance_multiplier'] for v in points],s=15,c=color,marker=marker,label=label if passed else None,alpha=.8)
        axes[1].scatter([v['exposure'] for v in rows],[max(1e-7,abs(v['bias_from_ideal_circuit'])) for v in rows],s=14,c=color,alpha=.8)
    axes[0].plot([0,11],[0,11],'k--',lw=1,label='PEC inverse-norm reference')
    axes[0].set(xlabel=r'Pauli exposure $\chi_{\rm PEC}$',ylabel=r'Measured $\log[F_{\rm variance}]$',title='Variance alone does not certify mitigation')
    axes[0].legend(fontsize=7)
    axes[1].axhline(.005,c='k',ls='--',lw=1);axes[1].set(yscale='log',xlabel=r'Pauli exposure $\chi_{\rm PEC}$',ylabel='Measured held-out bias',title='Same exposure, different residual bias')
    for d,color in ((9,'#157e85'),(25,'#ab4a45')):
        for mode,style in [('raw',':'),('ideal_pcec','--'),('noisy_pcec','-')]:
            rows=[v for v in data['star'] if v['p']==.001 and v['d']==d and v['mode']==mode]
            axes[2].loglog([v['theta'] for v in rows],[v['normalized_error'] for v in rows],style,color=color,label=f'd={d}, '+('raw' if mode=='raw' else 'ideal PCEC' if mode=='ideal_pcec' else 'noisy PCEC'))
    axes[2].set(xlabel=r'Angle $|\theta|$',ylabel=r'$q/(p|\theta|)$',title='STAR generator-fault submodel')
    axes[2].legend(fontsize=7)
    for ax in axes:ax.grid(alpha=.2)
    save(fig,'frontier_exposure')


def components():
    growth=read('growth_summary.json');delivery=read('delivery/results.json');fig,axes=plt.subplots(1,2,figsize=(10.8,4),layout='constrained')
    pairs=[(5,13),(9,13),(13,13),(5,25),(9,25),(13,25)]
    for basis,color,offset in [('x','#287ba0',-.12),('z','#b45e30',.12)]:
        rows=[next(v for v in growth['cases'] if (v['contract']['di'],v['contract']['df'],v['contract']['p'],v['contract']['basis'])==(di,df,.001,basis)) for di,df in pairs]
        for i,v in enumerate(rows):
            y=v['rate'] or v['simultaneous_upper95']
            axes[0].scatter(i+offset,y,c=color,marker='o' if v['failures'] else 'v',s=35,label=f'{basis.upper()}-basis probe' if i==0 else None)
            axes[0].plot([i+offset]*2,[v['rate'] if v['failures'] else y/1.7,v['simultaneous_upper95']],color=color,lw=1)
    axes[0].axhline(1e-5,color='gray',ls=':',label=r'$10^{-5}$ per component')
    axes[0].set(yscale='log',xticks=range(6),xticklabels=[f'{di}→{df}' for di,df in pairs],ylabel='Campaign error probability',xlabel='Initial → final distance',title=r'Growth + memory + readout, $p=10^{-3}$')
    axes[0].legend(fontsize=8)
    for ports,color in [(1,'#8253a5'),(4,'#178c70')]:
        rows=[v for v in delivery['rows'] if v['ports']==ports]
        axes[1].plot([v['batch'] for v in rows],[1000*v['seconds'] for v in rows],'o-',c=color,label=f'{ports} port(s), explicit SWAP delivery')
        axes[1].plot([v['batch'] for v in rows],[1000*v['abstract_delivery_seconds'] for v in rows],'o--',c=color,label=f'{ports} port(s), abstract delivery')
    axes[1].set(xscale='log',yscale='log',xticks=[1,4,8,16],xticklabels=['1','4','8','16'],xlabel='Hamming-weight batch size',ylabel='Onsite-stage duration (ms)',title='16-site compiled stage; bounded buffers')
    axes[1].legend(fontsize=7)
    for ax in axes:ax.grid(alpha=.2)
    save(fig,'frontier_components')


def locality():
    rows=read('locality_certificate.json')['rows'];local=json.loads((ROOT/'results/implementation_v1/locality.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(10.5,3.8),layout='constrained')
    for allowance,style in [(.0025,'-'),(.00125,'--')]:
        pts=[v for v in rows if v['spatial_allowance']==allowance]
        axes[0].loglog([v['tau'] for v in pts],[v['radius'] for v in pts],'o'+style,label=f'Allowance {allowance:g}')
    axes[0].set(xlabel=r'Time $\tau$',ylabel='Sufficient graph radius R',title='Loose state-independent certificate')
    for patch in ('cross5','3x3','4x3'):
        pts=sorted([v for v in local['comparisons'] if v['patch']==patch and v['finite_reference']=='4x4'],key=lambda v:v['tau'])
        axes[1].loglog([v['tau'] for v in pts],[v['difference'] for v in pts],'o-',label=patch+' versus 4x4')
    axes[1].axhline(.0025,c='gray',ls=':');axes[1].set(xlabel=r'Time $\tau$',ylabel='Measured finite-reference difference',title='Finite-patch checks are a different quantity')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
    save(fig,'frontier_locality')


def report():
    q=read('qubitization/results.json');growth=read('growth_summary.json');delivery=read('delivery/results.json')
    lines=['# Same-task algorithms and measured components','',
        'This pass keeps the previous capacity baseline frozen. New costs are separately identified algorithm and layout scenarios; new channel measurements constrain individual mechanisms.','',
        '## Qubitized dynamics on the same open-boundary Hamiltonian','',
        'The centered Hamiltonian has alpha = 4L(L−1) + L² at U/J=4: alpha=64 on 4×4. Repeating each onsite Pauli twice makes all coefficients magnitude 1/2. At 4×4 the 128-term PREPARE is seven Hadamards. A controlled block encoding uses 508 T and no arbitrary rotations, versus 508 T plus 256 rotations with generic amplitude preparation. The local adapter repairs SELECT signs/zero padding without modifying pyLIQTR.','',
        'The cosine/sine branches and their LCU give accepted amplitude 1/4 before amplification. Counts optimize zero through three amplification steps, simulation/synthesis/logical budgets, source count, and cleanup. Conditional fault bias and rejected attempts are charged. Large-degree phase generation and physical SELECT routing remain unvalidated. The two timing choices bracket extra Clifford serialization only; they are not rigorous hardware bounds.','',
        '| 4×4, tau=1 | Absolute error | T per attempted circuit | Q for one week |','|---|---:|---:|---:|']
    for eps in (.01,.001):
        for v in q['rows']:
            if v['L']==4 and v['tau']==1 and v['epsilon']==eps:
                name=v['variant']+(f", Clifford serialization={v['clifford_serial_fraction']:g}" if v['clifford_serial_fraction'] is not None else '')
                lines.append(f"| {name} | {eps:g} | {v['T_count']:,} | {10**v['log10_Q_min']:,.0f} |")
    comparisons=[]
    for v in q['rows']:
        if v['variant'].startswith('trotter') or v['clifford_serial_fraction']!=0:continue
        alternatives=[a for a in q['rows'] if a['variant'].startswith('trotter') and (a['L'],a['tau'],a['epsilon'])==(v['L'],v['tau'],v['epsilon'])]
        comparisons.append(v['log10_Q_min']<min(a['log10_Q_min'] for a in alternatives))
    lines+=['',f'QSVT beats the best of these two Trotter scenarios at {sum(comparisons)}/{len(comparisons)} tested size/time/precision points. This does not rule out better structured block encodings or coherent observable estimation. At tau=8 the Trotter coefficient still comes only from 3×3; its finite-size transfer is explicit.','',
        '## Exposure diagnostics','',
        'The measured NISQ variance multipliers are much smaller than whole-circuit PEC inverse-norm costs, but many held-out bias tests fail. Passing points and failures remain on the same axes. Training adds six million shots per fit and is retained in the table. These data do not collapse to a universal variance-plus-bias law in one exposure coordinate. STAR q/(p theta) is plotted from the enumerated generator-fault channels, with ideal and noisy PCEC separated; it is a submodel test, not full injection calibration.','',
        '## Larger decoded growth components','',
        'PBS smoke 7790942 and production 7790943 both completed F/0. The production census has 24 cases of one million shots (24 million total), all with pinned seeds, circuit hashes, ideal controls, and simultaneous one-sided 95% upper bounds over all 24 components.','',
        '| Initial→final distance | p | X failures | Z failures | Pauli-union upper bound |','|---|---:|---:|---:|---:|']
    for v in growth['pairs']:
        lines.append(f"| {v['di']}→{v['df']} | {v['p']:g} | {v['failures']['x']} | {v['failures']['z']} | {v['pauli_total_upper']:.3g} |")
    lines+=['','The X/Z union assumes both experiments describe one Pauli logical channel. For zero-failure pairs its bound is 1.235e-5, so this shot count does not certify a 1e-5 total floor. The distance-5 seed is visibly worse at p=1e-3. The starting state is ideally encoded; growth, target memory, and readout are included. No subtraction of the memory control or claim about a full injected-state channel is made.','',
        '## Explicit magic delivery','',
        'Magic states enter through perimeter ports, walk by nearest-neighbor SWAPs, are consumed, and return the data permutation by inverse SWAPs. Traversed data patches and surgery mediators remain reserved until return. One bounded output buffer per port and the entire square footprint (225 or 289 tiles) are charged, with the reference source area added. This is a valid conservative logical construction, not an optimized factory layout or physical surgery simulation.','',
        '| Batch | Ports | Abstract delivery (ms) | Explicit delivery (ms) | Array tiles |','|---:|---:|---:|---:|---:|']
    for v in delivery['rows']:lines.append(f"| {v['batch']} | {v['ports']} | {v['abstract_delivery_seconds']*1000:.3f} | {v['seconds']*1000:.3f} | {v['array_tiles']} |")
    lines+=['','The batch ranking changes with placement, transport, and port count. An aggregate HWP T-count saving therefore remains conditional until the complete chosen layout is scheduled. The previous capacity figure is not silently replaced by this isolated onsite stage.','',
        '## Centered locality certificate','',
        'After removing onsite U by an interaction picture, each two-spin hopping bond has norm 2J. The centered doublon has norm 1/2. Counting four first bonds and at most seven successors gives a sufficient omitted-boundary bound min(1, (2/7) sum[n≥R+1] (28 tau)^n/n!). At tau=1, the geometric majorant first meets 0.0025 at R=78 (a 157×157 enclosing square). This very loose certificate cannot justify replacing the full lattice by a 6×6 or 8×8 patch. It also does not say the physical light cone has radius 78. Finite-reference differences remain the useful necessary check, with outer-boundary error unbounded.','',
        '## What is settled and what remains scientific work','',
        'All six calculation deliverables in docs/FRONTIER_WORK.md are implemented and checked. The remaining evidence gaps are end-to-end noisy arbitrary-angle injection/growth, transfer of trained mitigation beyond eight qubits, tight state-dependent bulk truncation, and full routed SELECT/HWP schedules. These require additional protocol or physics input; passing this software audit does not settle them.','']
    (OUT/'summary.md').write_text('\n'.join(lines))


def main(check=False):
    validation=validate();(OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    if not check:algorithms();exposures();components();locality();report()
    print(validation)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');main(parser.parse_args().check)
