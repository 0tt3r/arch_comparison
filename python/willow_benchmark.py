"""Compare Willow estimates, ideal free dynamics, and sampling designs."""
import argparse
import csv
import gzip
import json
from math import exp

import numpy as np

from model import log_gamma2_depolarizing
from project_paths import ROOT
from willow_archive import OUT, read, save, sha


def load_reduced():
    return json.loads(gzip.decompress((OUT/'reduced.json.gz').read_bytes()))


def combine_paulis(z,zz,z_se=None,zz_se=None):
    """Z has shape (L,L,2); onsite ZZ has shape (L,L).

    Published covariances are absent. The triangle bound on propagated SE
    remains valid for any covariance consistent with the supplied marginals.
    It is not a total-error or Gaussian-coverage certificate.
    """
    L=z.shape[0];assert z.shape==(L,L,2) and zz.shape==(L,L)
    D=(1-z.sum(axis=2)+zz)/4
    signs=np.fromfunction(lambda x,y:(-1.)**(x+y),(L,L))
    M=np.sum(signs*(z[:,:,1]-z[:,:,0]))/(4*L*L)
    result=dict(local_d=D.ravel().tolist(),mean_d=[float(D.mean())],
                staggered_m=[float(M)],local_z=z.ravel().tolist())
    if z_se is not None:
        Dse=(z_se.sum(axis=2)+zz_se)/4
        result['se']=dict(local_d=Dse.ravel().tolist(),mean_d=[float(Dse.mean())],
            staggered_m=[float(z_se.sum()/(4*L*L))],local_z=z_se.ravel().tolist())
    return result


def published_arrays(case):
    L=case['L'];buckets={}
    for v in case['published']:
        key=(v['method'],round(v['tau'],8))
        if key not in buckets:
            buckets[key]={s:np.full(shape,np.nan) for s,shape in
                          [('z',(L,L,2)),('zz',(L,L)),('z_se',(L,L,2)),('zz_se',(L,L))]}
        b=buckets[key];x,y=v['sites'][0]
        if v['kind']=='Z':
            s=0 if v['spins'][0]=='up' else 1
            assert np.isnan(b['z'][x,y,s]),'Duplicate observable'
            b['z'][x,y,s]=v['value'];b['z_se'][x,y,s]=v['se']
        else:
            assert np.isnan(b['zz'][x,y]),'Duplicate onsite ZZ'
            b['zz'][x,y]=v['value'];b['zz_se'][x,y]=v['se']
    assert all(np.isfinite(a).all() for b in buckets.values() for a in b.values())
    return {key:combine_paulis(**b) for key,b in buckets.items()}


def free_dynamics(L,flux,occupations,times):
    """Exact continuous U=0 dynamics for the declared open, product input.

    Landau gauge: vertical links exp(i phi x); local densities are gauge
    invariant for the occupation-basis input. No many-body approximation.
    """
    H=np.zeros((L*L,L*L),complex)
    for x in range(L):
        for y in range(L):
            a=x*L+y
            for xx,yy,value in [(x+1,y,-1.),(x,y+1,-np.exp(1j*flux*x))]:
                if xx<L and yy<L:
                    b=xx*L+yy;H[a,b]=value;H[b,a]=np.conjugate(value)
    E,V=np.linalg.eigh(H)
    initial=np.zeros((L*L,2))
    for v in occupations:
        x,y=v['site'];initial[x*L+y,0 if v['spin']=='up' else 1]=v['n']
    result={}
    for t in times:
        R=(V*np.exp(-1j*t*E))@V.conj().T
        density=(np.abs(R)**2@initial).reshape(L,L,2)
        z=1-2*density
        result[round(t,8)]=combine_paulis(z,z[:,:,0]*z[:,:,1])
    return result


def symmetrize_shot_means(mean,L):
    values=np.array(mean);N=L*L
    d=values[:N].reshape(L,L);z=values[N+2:].reshape(L,L,2)
    return dict(local_d=((d+d.T)/2).ravel().tolist(),mean_d=[values[N]],
                staggered_m=[values[N+1]],local_z=((z+z.transpose(1,0,2))/2).ravel().tolist())


def compare_cases(data):
    metrics=[];audits=[];curves=[]
    for case in data['cases']:
        L=case['L'];key={s:case[s] for s in ('L','U','flux','steps')}
        arrays=published_arrays(case)
        times=[v['tau'] for v in case['shots']]
        if case['U']==0:
            free=free_dynamics(L,0 if case['flux']=='0' else np.pi,case['initial_occupations'],times)
            references={'continuous FLO':free,
                        'published FLO circuit':{t:arrays['exact',t] for t in times}}
            # A second, approximate continuum reference provides an independent
            # gauge/time-convention check, not a replacement for exact FLO.
            refs={t:arrays['TDVP chi=2048',t] for t in times}
            gauge_check=max(abs(free[t]['mean_d'][0]-refs[t]['mean_d'][0]) for t in times if t<=.2)
            short_circuit_check=max(abs(np.array(free[min(times)]['local_z'])-
                arrays['exact',min(times)]['local_z']))
        else:
            references={method:{t:arrays[method,t] for t in times}
                        for method in ('TDVP chi=1024','TDVP chi=2048')}
            gauge_check=None
            short_circuit_check=None
        discrepancies={}
        for filter_name in ('unfiltered','strict','relaxed1','relaxed2'):
            errors=[]
            for v in case['shots']:
                reconstructed=symmetrize_shot_means(v['statistics'][filter_name]['mean'],L)
                errors.extend(np.array(reconstructed['local_z'])-arrays['raw',v['tau']]['local_z'])
                errors.extend(np.array(reconstructed['local_d'])-arrays['raw',v['tau']]['local_d'])
            discrepancies[filter_name]=float(max(abs(x) for x in errors))
        best=min(discrepancies,key=discrepancies.get)
        audits.append(dict(**key,raw_pipeline_max_differences=discrepancies,
            raw_filter_reproduced=best if discrepancies[best]<1e-10 else None,
            free_vs_TDVP_early_mean_difference=gauge_check,
            free_vs_circuit_first_time_max_Z_difference=None if short_circuit_check is None else float(short_circuit_check),
            empty_groups=case['empty_groups'],initial_state_json_audit=case['initial_state_json_audit']))
        for reference,ref in references.items():
            methods=['raw','mitigated']+(['exact'] if case['U']==0 and reference=='continuous FLO' else [])
            for method in methods:
                for observable in ('local_d','mean_d','staggered_m','local_z'):
                    diff=np.array([np.array(arrays[method,t][observable])-ref[t][observable] for t in times])
                    se=np.array([arrays[method,t]['se'][observable] for t in times])
                    metrics.append(dict(**key,method=method,observable=observable,reference=reference,
                        points=diff.size,rmse=float(np.sqrt(np.mean(diff**2))),max_difference=float(np.max(np.abs(diff))),
                        median_reported_se_or_bound=float(np.median(se)),
                        se_kind='reported marginal SE' if observable=='local_z' else 'covariance-agnostic upper bound from reported marginal SEs',
                        fraction_within_01=float(np.mean(np.abs(diff)<=.01)),
                        fraction_within_001=float(np.mean(np.abs(diff)<=.001))))
        for (method,t),v in arrays.items():
            for obs in ('mean_d','staggered_m'):
                curves.append(dict(**key,method=method,tau=t,observable=obs,value=v[obs][0],se_bound=v['se'][obs][0]))
        if case['U']==0:
            for t,v in free.items():
                for obs in ('mean_d','staggered_m'):
                    curves.append(dict(**key,method='continuous FLO',tau=t,observable=obs,value=v[obs][0],se_bound=0.))
    return metrics,audits,curves


def design_forecast(between,within,K,shot_multiplier=1,twirl_multiplier=1):
    """Random-effects planning sensitivity, not a fitted error-mitigation law."""
    latent=np.maximum(0.,np.asarray(between)-np.asarray(within))
    return np.sqrt((latent+np.asarray(within)/shot_multiplier)/(K*twirl_multiplier))


def sampling_table(data,cfg):
    rows=[]
    for case in data['cases']:
        if case['U']!=4 or case['flux']!='0' or case['steps']!=3:continue
        L=case['L'];N=L*L
        shot=next(v for v in case['shots'] if v['tau']==1.)
        for filter_name in ('unfiltered','strict','relaxed1','relaxed2'):
            s=shot['statistics'][filter_name]
            for obs,i in [('central_d',(L//2-1)*L+(L//2-1)),('mean_d',N),('staggered_m',N+1)]:
                row=dict(L=L,U=4,flux='0',steps=3,tau=1.,filter=filter_name,observable=obs,
                    attempted_shots=shot['shots'],acceptance=s['acceptance'],
                    cluster_se=s['cluster_se'][i],estimated_fixed_twirl_floor=s['fixed_twirl_infinite_shots_se'][i])
                for h in (1,10,100):
                    row[f'same_twirls_shots_x{h}_se']=float(design_forecast(s['between_mean_variance'][i],s['within_mean_variance'][i],cfg['twirls'],shot_multiplier=h))
                    row[f'new_twirls_x{h}_se']=float(design_forecast(s['between_mean_variance'][i],s['within_mean_variance'][i],cfg['twirls'],twirl_multiplier=h))
                rows.append(row)
    return rows


def resource_table(data,cfg):
    rows=[]
    for L in cfg['sizes']:
        for r in cfg['steps']:
            cases=[v for v in data['cases'] if v['L']==L and v['steps']==r]
            t=cfg['table_I'][str(L)];G=t['gates_per_step']*r+t['gates_offset']
            rows.append(dict(L=L,steps=r,physical_qubits=t['physical_qubits'],data_qubits=2*L*L,
                gates=G,depth=t['depth_per_step']*r+t['depth_offset'],
                shots_per_time=cfg['twirls']*cfg['shots_per_twirl'],times_per_series=cfg['expected_times'],
                shots_per_series=cfg['twirls']*cfg['shots_per_twirl']*cfg['expected_times'],
                strict_acceptance_median=float(np.median([s['statistics']['strict']['acceptance'] for c in cases for s in c['shots']])),
                gatewise_2q_PEC_variance_at_p001=exp(G*log_gamma2_depolarizing(.001,2)),
                scope='Published gate counts and actual shot census. PEC is a specified hypothetical channel at p=0.001, not a model fitted to Willow.'))
    return rows


def csv_rows(name,rows):
    with (OUT/name).open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def run(check=False):
    data=load_reduced();cfg=read(ROOT/'configs/willow_benchmark.json')
    for path,value in data['inputs'].items():assert sha(ROOT/path)==value,path
    downloads=read(OUT/'downloads.json')
    assert downloads['downloader_sha256']==sha(ROOT/'python/willow_archive.py')
    assert downloads['config_sha256']==sha(ROOT/'configs/willow_benchmark.json')
    assert downloads['index_sha256']==sha(OUT/'archive_index.json')
    initial=next(v for v in downloads['files'] if v['name'].endswith('/initial_states.json'))
    assert initial['sha256']==sha(OUT/'initial_states_supplied.json')
    metrics,audits,curves=compare_cases(data)
    result=dict(comparisons=metrics,raw_audits=audits,curves=curves,
                sampling=sampling_table(data,cfg),resources=resource_table(data,cfg))
    validation=dict(status='PASS: declared archive census, reference separation, and sampling arithmetic',
        cases=len(data['cases']),shot_batches=sum(len(c['shots']) for c in data['cases']),
        attempted_shots=sum(s['shots'] for c in data['cases'] for s in c['shots']),
        raw_pipeline_replayed=sum(v['raw_filter_reproduced'] is not None for v in audits),
        comparisons=len(metrics),interacting_exact_groups_empty=sum('exact' in v['empty_groups'] for v in audits),
        maximum_early_free_reference_difference=max(v['free_vs_TDVP_early_mean_difference'] or 0 for v in audits),
        maximum_first_time_circuit_difference=max(v['free_vs_circuit_first_time_max_Z_difference'] or 0 for v in audits),
        scope='No interacting exact truth is supplied. Published processed estimates are reused. Design extrapolations use 20 observed twirls and are not bias or accuracy certificates.')
    assert validation['cases']==36 and validation['shot_batches']==864 and validation['attempted_shots']==172800000
    assert validation['interacting_exact_groups_empty']==24
    assert validation['raw_pipeline_replayed']==36
    assert validation['maximum_first_time_circuit_difference']<1e-4
    inputs={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'python/willow_benchmark.py',
        ROOT/'python/willow_data.py',ROOT/'python/willow_archive.py',ROOT/'python/phasecraft_data.py',
        ROOT/'configs/willow_benchmark.json',ROOT/cfg['source_pdf'],
        OUT/'reduced.json.gz',OUT/'downloads.json',OUT/'archive_index.json',OUT/'initial_states_supplied.json']}
    if check:
        assert read(OUT/'inputs.json')==inputs
        assert read(OUT/'benchmarks.json')==result
        assert read(OUT/'validation.json')==validation
        print(validation['status'])
    else:
        save(OUT/'inputs.json',inputs);save(OUT/'benchmarks.json',result);save(OUT/'validation.json',validation)
        csv_rows('comparisons.csv',metrics);csv_rows('sampling.csv',result['sampling']);csv_rows('resources.csv',result['resources'])
        print(json.dumps(validation,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true')
    run(p.parse_args().check)
