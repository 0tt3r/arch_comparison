"""Replay, summarize, and plot the four bounded implementation checks."""
import argparse
import csv
import hashlib
import importlib.metadata
import json
from math import log10

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

from nisq_estimator_study import evaluate_estimator
from project_paths import ROOT

OUT=ROOT/'results/implementation_v1'
METHODS={'raw':'Raw','number':'Number selection','tflo_raw':'Trained affine','tflo_number':'Selection + trained affine'}
COLORS=['#777777','#d18b00','#2463a8','#008c79']


def read(name):return json.loads((OUT/name).read_text())


def training_for(data,test):
    return [v for v in data['cases'] if v['role']=='train' and v['p']==test['p'] and v['r']==test['r']]


def bootstrap_training(data):
    """Binomial training campaigns independently check delta-method variance."""
    rng=np.random.default_rng(20260923);rows=[];M=data['configuration']['training_shots_per_circuit'];B=8000
    for test in data['cases']:
        if test['role']!='test' or (test['u'],test['tau'],test['noise'])!=(4.,1.,'matched'):continue
        train=training_for(data,test)
        for method in ('raw','number'):
            x=[]
            for point in train:
                v=point['noisy'][method];accepted=rng.binomial(M,min(1.,v['acceptance']),B)
                assert accepted.min()>0
                x.append(rng.binomial(accepted,v['mean'])/accepted)
            x=np.array(x).T;y=np.array([v['ideal'] for v in train]);xm=x.mean(axis=1)
            slope=((x-xm[:,None])*(y-y.mean())).sum(axis=1)/((x-xm[:,None])**2).sum(axis=1)
            predictions=y.mean()+slope*(test['noisy'][method]['mean']-xm)
            result=evaluate_estimator(test,train,'tflo_'+method)
            measured=float(np.var(predictions,ddof=1));ratio=measured/result['training_variance']
            rows.append(dict(p=test['p'],r=test['r'],method=method,repeats=B,seed=20260923,
                analytical_variance=result['training_variance'],binomial_variance=measured,ratio=ratio,
                finite_training_shift=float(predictions.mean()-result['prediction'])))
    return rows


def validate():
    nisq=read('nisq.json');star=read('star_channels.json');ft=read('ft_stage/results.json');local=read('locality.json');growth=read('star_growth/results.json')
    assert len(nisq['cases'])==144 and len(nisq['estimates'])==384
    keys={v['key']:v for v in nisq['cases']};assert len(keys)==144
    for row in nisq['estimates']:
        test=keys[row['key']];new=evaluate_estimator(test,training_for(nisq,test),row['method'])
        assert all(row[k]==v for k,v in new.items())
    for record in (nisq,star,ft,local):
        for name,digest in record['source_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    assert growth['source_sha256']==hashlib.sha256((ROOT/'python/star_growth_probe.py').read_bytes()).hexdigest()
    assert len(star['rows'])==1640 and len(star['campaigns'])==384
    for row in star['rows']:
        assert -.000000000001<=row['q_dephasing']<=.5 and row['mean_seconds']>=0
        assert np.linalg.norm(row['multiplier'])<=1+1e-12
    for row in star['campaigns']:
        assert row['rus_tail_bound']<.011*row['epsilon']
        assert np.isfinite(row['log_sampling_time'])
    assert len(ft['rows'])==16
    for v in ft['rows']:
        assert v['clean_ancillas']==v['batch']-v['batch'].bit_count()
        assert v['T_count']==v['compiled_counts']['actions'].get('qlx_standard_t',0)+v['compiled_counts']['actions'].get('qlx_standard_tdg',0)
        assert v['logical_clocks']>=v['source_only_clocks']
        if v['buffer_policy']=='bounded':assert v['required_buffer_upper_estimate']<=v['ports']
    assert len(local['rows'])==23 and len(local['comparisons'])==33
    for row in local['comparisons']:
        a=next(v for v in local['rows'] if v['patch']==row['patch'] and v['tau']==row['tau'])
        b=next(v for v in local['rows'] if v['patch']==row['finite_reference'] and v['tau']==row['tau'])
        assert row['difference']==abs(a['doublon']-b['doublon'])
        assert all(passed==(row['difference']<=float(budget)) for budget,passed in row['passes'].items())
    assert len(growth['rows'])==20 and sum(v['shots'] for v in growth['rows'])==10000000
    boot=bootstrap_training(nisq)
    assert len(boot)==16 and all(.9<v['ratio']<1.1 for v in boot)
    return dict(status='PASS: source hashes, estimator replay, binomial fit checks, channel and finite-patch accounting',
                nisq_cases=144,held_out_estimates=384,training_bootstraps=16,star_channels=1640,
                star_workloads=384,growth_cases=20,growth_shots=10000000,ft_stages=16,locality_comparisons=33),boot


def savefig(fig,name):
    fig.savefig(ROOT/f'figures/{name}.pdf',bbox_inches='tight')
    fig.savefig(ROOT/f'figures/{name}.png',dpi=180,bbox_inches='tight');plt.close(fig)


def plot_nisq(data):
    fig,axes=plt.subplots(1,3,figsize=(13,3.8),layout='constrained')
    for ax,p in zip(axes[:2],(.001,.0001)):
        for method,color in zip(METHODS,COLORS):
            for noise,style in [('matched','-'),('mismatch','--')]:
                rows=sorted([v for v in data['estimates'] if v['method']==method and v['p']==p and v['u']==4 and v['tau']==1 and v['noise']==noise],key=lambda v:v['r'])
                ax.loglog([v['r'] for v in rows],[abs(v['bias_from_ideal_circuit']) for v in rows],style+'o',color=color,label=METHODS[method] if noise=='matched' else None,ms=4)
        ax.axhline(.005,color='black',ls=':',label='Bias allowance' if p==.001 else None)
        ax.set(title=f'2×2 held-out U=4, τ=1; p={p:g}',xlabel='Second-order steps r',ylabel='Absolute residual mitigation bias')
        ax.set_xticks([4,8,16,32],['4','8','16','32']);ax.grid(alpha=.2)
        ax.xaxis.set_minor_formatter(NullFormatter())
    axes[0].legend(fontsize=7,loc='upper left')
    for method,color in zip(METHODS,COLORS):
        rows=sorted([v for v in data['estimates'] if v['method']==method and v['p']==.001 and v['u']==4 and v['tau']==1 and v['noise']=='matched'],key=lambda v:v['r'])
        axes[2].loglog([v['seconds_including_training'] for v in rows],[abs(v['bias_from_ideal_circuit'])+abs(v['ideal']-v['exact'])+.0025 for v in rows],'-o',color=color,ms=4)
        for v in rows:axes[2].annotate(str(v['r']),(v['seconds_including_training'],abs(v['bias_from_ideal_circuit'])+abs(v['ideal']-v['exact'])+.0025),fontsize=7,xytext=(3,2),textcoords='offset points')
    axes[2].axhline(.01,color='black',ls=':');axes[2].set(title='Matched noise; labels are r',xlabel='One-copy seconds, including training',ylabel='|Mitigation bias| + |Trotter bias| + 0.0025');axes[2].grid(alpha=.2)
    fig.suptitle('Exact simulated channel; dashed lines: 20% p₂ shift + 1% onsite overrotation',fontsize=11)
    savefig(fig,'implementation_nisq')


def plot_star(data,growth):
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for d in (3,5,9,15,25):
        rows=[v for v in data['rows'] if v['p']==.001 and v['d']==d and v['mode']=='raw']
        axes[0,0].loglog([v['theta'] for v in rows],[v['acceptance'] for v in rows],label=f'dᵢ={d}')
    axes[0,0].set(title='Generator-fault projection model',xlabel='Logical angle |θ|',ylabel='Single-projection acceptance');axes[0,0].legend(fontsize=8)
    for mode,label in [('raw','Uncorrected'),('ideal_pcec','Ideal PCEC'),('noisy_pcec','Noisy PCEC'),('switch','Noisy PCEC + switching')]:
        rows=[v for v in data['rows'] if v['p']==.001 and v['d']==9 and v['mode']==mode]
        axes[0,1].loglog([v['theta'] for v in rows],[max(1e-15,v['q_dephasing']) for v in rows],label=label)
    axes[0,1].axhline(4*.001/15,color='black',ls=':',lw=1,label='Original 4p/15 guide')
    axes[0,1].set(title='RUS channel, dᵢ=9, no growth floor',xlabel='Logical angle |θ|',ylabel='q = (1 − Re c)/2');axes[0,1].legend(fontsize=7)
    pairs=[(3,3),(3,5),(3,7),(5,5),(5,7)]
    for basis,offset,color in [('x',-.1,'#2463a8'),('z',.1,'#d18b00')]:
        rows=[next(v for v in growth['rows'] if (v['di'],v['df'],v['p'],v['basis'])==(*pair,.001,basis)) for pair in pairs]
        y=np.array([v['rate'] for v in rows]);ci=np.array([v['confidence_interval'] for v in rows]).T
        axes[1,0].errorbar(np.arange(5)+offset,y,yerr=[y-ci[0],ci[1]-y],fmt='o',capsize=3,color=color,label=f'{basis.upper()} memory')
    axes[1,0].axhline(1e-5,color='black',ls=':',label='Assumed floor in capacity model')
    axes[1,0].set_xticks(range(5),[f'{a}→{b}' for a,b in pairs]);axes[1,0].set_yscale('log')
    axes[1,0].set(title='Circuit-level growth probes at p=0.001',xlabel='Small-patch → target distance',ylabel='Memory failure probability');axes[1,0].legend(fontsize=7)
    for floor in (0,1e-5,1e-4):
        rows=sorted([v for v in data['campaigns'] if v['epsilon']==.001 and v['p']==.001 and v['floor_at_d9']==floor and v['detector_factor']==1 and v['mode']=='switch'],key=lambda v:v['d'])
        axes[1,1].plot([v['d'] for v in rows],[v['log_sampling_time']/np.log(10) for v in rows],'-o',label=f'Floor at dᵢ=9: {floor:g}',ms=4)
    axes[1,1].set(title='Strict 4×4 workload; declared growth scenarios',xlabel='Injection distance dᵢ (data distance 25)',ylabel='log₁₀(serial attempt seconds × exp χ)');axes[1,1].legend(fontsize=7)
    for ax in axes.flat:ax.grid(alpha=.2)
    savefig(fig,'implementation_star')


def plot_ft(data):
    fig,axes=plt.subplots(1,3,figsize=(12,3.7),layout='constrained')
    reference=[v for v in data['rows'] if v['ports']==1 and v['buffer_policy']=='bounded']
    axes[0].plot([v['batch'] for v in reference],[v['T_count'] for v in reference],'-o',label='T count')
    axes[0].plot([v['batch'] for v in reference],[v['routed_CX'] for v in reference],'-s',label='Routed CX count');axes[0].set(ylabel='Gates',xlabel='Hamming-weight batch',title='One 16-site onsite layer');axes[0].legend(fontsize=8)
    for ports,color in [(1,'#2463a8'),(4,'#008c79')]:
        for policy,style in [('bounded','-'),('continuous','--')]:
            rows=[v for v in data['rows'] if v['ports']==ports and v['buffer_policy']==policy]
            axes[1].plot([v['batch'] for v in rows],[1000*v['seconds'] for v in rows],style+'o',color=color,label=f'{ports} ports, {policy}',ms=4)
            axes[2].plot([v['batch'] for v in rows],[v['packed_patches']+v['required_buffer_upper_estimate'] for v in rows],style+'o',color=color,ms=4)
    axes[1].set(title='NN routing + source + feedforward',xlabel='Hamming-weight batch',ylabel='Stage duration (ms)');axes[1].legend(fontsize=8)
    axes[2].set(title='Data, route, and state-buffer patches',xlabel='Hamming-weight batch',ylabel='Logical patches (sources excluded)')
    for ax in axes:ax.set_xticks([1,4,8,16]);ax.grid(alpha=.2)
    savefig(fig,'implementation_ft')


def plot_locality(data):
    fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    for name in ('cross5','2x2','2x3','3x3','4x3'):
        rows=sorted([v for v in data['comparisons'] if v['patch']==name and v['finite_reference']=='3x3'],key=lambda v:v['tau'])
        if rows:axes[0].loglog([v['tau'] for v in rows],[v['difference'] for v in rows],'-o',label=name,ms=4)
        rows=sorted([v for v in data['comparisons'] if v['patch']==name and v['finite_reference']=='4x4'],key=lambda v:v['tau'])
        if rows:axes[1].semilogy([v['tau'] for v in rows],[v['difference'] for v in rows],'-o',label=name,ms=4)
    for ax in axes:
        ax.axhline(.0025,color='black',ls=':',label='0.0025 spatial allowance');ax.axhline(.00125,color='gray',ls=':',label='0.00125 allowance')
        ax.set(xlabel='Evolution time τ',ylabel='|D(patch) − D(finite reference)|');ax.grid(alpha=.2);ax.legend(fontsize=8)
    axes[0].set_title('Nested patches against 3×3');axes[1].set_title('Nested patches against 4×4');axes[1].set_xticks([1,2])
    fig.suptitle('Same coordinate (1,1), Néel state, U=4: finite references retain nearby boundaries',fontsize=11)
    savefig(fig,'implementation_locality')


def summaries(nisq,star,ft,local):
    rows=[]
    for p in (.001,.0001):
        for noise in ('matched','mismatch'):
            for method in METHODS:
                cases=[v for v in nisq['estimates'] if (v['p'],v['noise'],v['method'])==(p,noise,method)]
                rows.append(dict(p=p,noise=noise,method=method,tests=len(cases),bias_passes=sum(v['passes_bias_allowance'] for v in cases),total_passes=sum(v['passes_total_error'] for v in cases),maximum_bias=max(abs(v['bias_from_ideal_circuit']) for v in cases)))
    with (OUT/'nisq_summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for name,data in [('ft_stages',ft['rows']),('star_workloads',star['campaigns']),('locality_comparisons',local['comparisons'])]:
        flat=[{k:v for k,v in row.items() if not isinstance(v,(dict,list))} for row in data]
        with (OUT/f'{name}.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    lines=['# Bounded implementation checks','',
        'These experiments test estimator bias, a specified injection submodel, an explicit FT layer, and finite-patch truncation. The capacity curves remain conditional; these data are not silently substituted for their assumptions.','',
        '## NISQ: trained mitigation with independent interacting tests','',
        'Exact eight-qubit channels include native routing, depolarizing gates, idle noise and SPAM. Six U=0 circuits train each affine correction; U=2,4 circuits at other times are held out. Zero-angle onsite slots retain noise. Each training circuit receives one million attempted shots in the cost model.','',
        '| p | Test noise | Estimator | Bias ≤ 0.005 | Total error plan ≤ 0.01 |','|---|---|---|---:|---:|']
    for v in rows:lines.append(f"| {v['p']:g} | {v['noise']} | {METHODS[v['method']]} | {v['bias_passes']}/{v['tests']} | {v['total_passes']}/{v['tests']} |")
    lines += ['', 'The total plan adds absolute mitigation bias, absolute Trotter bias and a 0.0025 Gaussian statistical allowance; it does not exploit bias cancellation. The binomial training bootstrap checks the delta-method fit variance. It does not establish hardware coverage or transfer to larger lattices.', '',
        '## STAR: channel, acceptance, and measured growth components','',
        'The generator-Pauli fault channel is exact under perfect syndrome projection. Ideal PCEC makes that channel Pauli; one level of noisy correction leaves a nonzero coherent component, which is reported separately. RUS correction angles, their preparation times, and their growth errors are included. Workload truncation tails consume less than 1.1% of epsilon. Detector rejection and high-distance growth remain declared scenarios.', '',
        'Twenty circuit-level X/Z memory probes use 500,000 shots each, ideal small-patch preparation, noisy direct growth, target-distance stabilization, and final readout. At p=0.001, d=5→7 gives 42 and 41 failures: 8.4e-5 and 8.2e-5. These include memory/readout and are not the error of an injected arbitrary-angle state. At p=0.0001, zero failures in those two rows give two-sided 95% upper limits 7.38e-6, not zero error. No high-distance floor is calibrated by these probes.', '',
        '## FT: compiled arithmetic and constructive nearest-neighbor scheduling','',
        'CUDA-Q Logical compiles the explicit four-T clean-AND arithmetic, X-measurement uncomputation and conditional CZs. All two-patch gates are routed to adjacent square-grid patches. For one 4×4 onsite layer, batch 16 reduces 1280 T to 460 T but increases routed CX from 32 to 1524. With one source port and one buffered state, time falls from 46.824 ms to 27.510 ms. With four ports it rises from 12.282 ms to 16.398 ms. Unlimited buffering helps time but can require hundreds of extra patches.', '',
        'This is a valid logical dependency schedule under declared CX, measurement and feedforward macros. Magic-state transport paths are abstracted. It is not a physical lattice-surgery layout, decoded error calibration, optimized router, or full dynamics campaign.', '',
        '## Locality: finite-reference differences','',
        'At tau=1, the 3×3 versus 4×4 doublon difference is 0.04015; 4×3 versus 4×4 is 0.01488. Every tested proper subpatch of 4×4 fails the 0.0025 allowance at tau=1 and tau=2. Similar Trotter coefficients therefore do not certify a small truncation error. These references share coordinate (1,1), retain nearby outer boundaries, and do not establish a bulk limit.', '',
        'See [the manuscript](../../tex/main.pdf), [methods and reruns](../../docs/IMPLEMENTATION_CHECKS.md), and the four `implementation_*.pdf` figures.']
    (OUT/'summary.md').write_text('\n'.join(lines)+'\n')


def run(check=False):
    validation,boot=validate()
    if check:
        assert read('validation.json')==validation
        assert read('training_bootstrap.json')==boot
        print(json.dumps(validation,indent=2));return
    nisq=read('nisq.json');star=read('star_channels.json');growth=read('star_growth/results.json');ft=read('ft_stage/results.json');local=read('locality.json')
    plot_nisq(nisq);plot_star(star,growth);plot_ft(ft);plot_locality(local);summaries(nisq,star,ft,local)
    (OUT/'training_bootstrap.json').write_text(json.dumps(boot,indent=2)+'\n')
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    (OUT/'environment.json').write_text(json.dumps({name:importlib.metadata.version(name) for name in ('numpy','scipy','numba','matplotlib')},indent=2)+'\n')
    print(json.dumps(validation,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args();run(args.check)
