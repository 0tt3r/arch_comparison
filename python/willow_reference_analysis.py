"""Separate measured experimental discrepancy, circuit error, and reference error."""
import argparse
import csv
import json

import numpy as np

from project_paths import ROOT
from small_dynamics import rectangle
from willow_archive import read, save, sha
from willow_benchmark import load_reduced, published_arrays
from willow_reference import OUT, free_circuit
from willow_sector import sources


def observables(result,L=4,symmetrize=True):
    d=np.zeros((L,L));z=np.zeros((L,L,2))
    for i,(x,y) in enumerate(rectangle(L,L)):
        d[x,y]=result['doublons'][i]
        z[x,y]=[result['z_up'][i],result['z_down'][i]]
    if symmetrize:d=(d+d.T)/2;z=(z+z.transpose(1,0,2))/2
    signs=np.fromfunction(lambda x,y:(-1.)**(x+y),(L,L))
    return dict(local_d=d.ravel(),mean_d=np.array([d.mean()]),local_z=z.ravel(),
                staggered_m=np.array([np.sum(signs*(z[:,:,1]-z[:,:,0]))/(4*L*L)]))


def load_records():
    cfg=read(ROOT/'configs/willow_reference.json');build=read(OUT/'build.json');current=sources()
    assert build['sources']==current
    continuous={};circuits={};all_paths=[];records=[]
    for folder in ('continuous','circuit','direct_check'):
        paths=sorted((OUT/folder).glob('*.json'));all_paths.extend(paths)
        for path in paths:
            v=read(path);c=v['contract'];result=v['result'];records.append(v)
            assert v['status']=='complete' and c['sources']==current
            assert c['binary_sha256']==build['binary_sha256']
            assert (c['L'],c['U'],c['flux'])==(4,4,0)
            assert abs(result['norm_squared']-1)<1e-9
            assert result['dimension']==147232800
            assert v['pbs_jobid']
            prefix=path.name.rsplit('_k',1)[0]
            assert sha(path.parent/(prefix+'.coefficients.txt'))==v['coefficient_sha256']
            assert sha(path.parent/(prefix+'.order.txt'))==v['order_sha256']
            tau=round(result['tau'],8)
            if folder=='continuous':
                assert tau not in continuous and c['mode']=='trajectory'
                continuous[tau]=observables(result)
            elif folder=='circuit':
                key=(tau,c['r']);assert key not in circuits and c['mode']=='circuit'
                circuits[key]=observables(result)
    assert set(continuous)==set(cfg['times'])
    assert set(circuits)=={(t,r) for t in cfg['times'] for r in cfg['steps']}
    direct=[v for v in records if v['contract']['mode']=='exact']
    assert len(direct)==1 and direct[0]['result']['tau']==1.
    direct=observables(direct[0]['result'])
    discrepancy=max(np.max(abs(direct[k]-continuous[1.][k])) for k in direct)
    assert discrepancy<1e-9
    norms=max(abs(v['result']['norm_squared']-1) for v in records)
    tails=max(v['state_tail_bound'] or 0 for v in records)
    assert tails<1e-11
    return continuous,circuits,dict(reference_records=len(records),continuous_times=len(continuous),
        circuit_points=len(circuits),maximum_norm_error=norms,maximum_accumulated_state_tail=tails,
        direct_vs_incremental_max_observable_difference=float(discrepancy)),all_paths


def decompose(experiment,circuit,continuous):
    """Signed identity; experimental discrepancies include sampling and processing."""
    e=np.asarray(experiment)-circuit;t=np.asarray(circuit)-continuous;total=e+t
    return dict(experiment_vs_circuit_rms=float(np.sqrt(np.mean(e*e))),
        circuit_vs_continuous_rms=float(np.sqrt(np.mean(t*t))),
        experiment_vs_continuous_rms=float(np.sqrt(np.mean(total*total))),
        twice_cross_mean=float(2*np.mean(e*t)),
        experiment_vs_continuous_max=float(np.max(abs(total))),
        experiment_vs_circuit_max=float(np.max(abs(e))),
        fraction_total_within_01=float(np.mean(abs(total)<=.01)),
        fraction_total_within_001=float(np.mean(abs(total)<=.001)))


def analyze():
    continuous,circuits,validation,paths=load_records()
    for name in ('smoke','production'):
        accounting=read(OUT/f'{name}_accounting.json')
        assert accounting['job_state']=='F' and accounting['Exit_status']==0
        assert accounting['Resource_List']['ncpus']==128
    data=load_reduced();cfg=read(ROOT/'configs/willow_reference.json')
    metrics=[];curves=[];decomposition=[];reference_comparison=[]
    for r in cfg['steps']:
        case=next(c for c in data['cases'] if (c['L'],c['U'],c['flux'],c['steps'])==(4,4,'0',r))
        pub=published_arrays(case)
        for obs in ('local_d','mean_d','local_z','staggered_m'):
            ideal=np.array([circuits[t,r][obs] for t in cfg['times']])
            target=np.array([continuous[t][obs] for t in cfg['times']])
            for method in ('raw','mitigated'):
                measured=np.array([pub[method,t][obs] for t in cfg['times']])
                m=decompose(measured,ideal,target)
                np.testing.assert_allclose(m['experiment_vs_continuous_rms']**2,
                    m['experiment_vs_circuit_rms']**2+m['circuit_vs_continuous_rms']**2+m['twice_cross_mean'],atol=1e-15)
                metrics.append(dict(steps=r,observable=obs,method=method,points=measured.size,**m))
                if obs=='local_d':
                    for i,t in enumerate(cfg['times']):
                        decomposition.append(dict(tau=t,steps=r,method=method,**decompose(measured[i],ideal[i],target[i])))
            for method in ('TDVP chi=1024','TDVP chi=2048'):
                delta=np.array([pub[method,t][obs] for t in cfg['times']])-target
                reference_comparison.append(dict(steps=r,observable=obs,method=method,
                    continuous_reference_rms=float(np.sqrt(np.mean(delta**2))),
                    continuous_reference_max=float(np.max(abs(delta)))))
        for t in cfg['times']:
            values={'continuous':continuous[t],'ideal circuit':circuits[t,r],
                    **{m:pub[m,t] for m in ('raw','mitigated','TDVP chi=1024','TDVP chi=2048')}}
            for method,v in values.items():
                curves.append(dict(tau=t,steps=r,method=method,mean_d=float(v['mean_d'][0]),
                    central_d=float(v['local_d'][5]),local_d=list(v['local_d']),local_z=list(v['local_z'])))
    # Replay the cheap, independent U=0 ordering check at every saved point.
    audit=read(OUT/'circuit_audit.json');zerr=[];derr=[]
    for c in data['cases']:
        if c['U']!=0:continue
        pub=published_arrays(c);L=c['L']
        for s in c['shots']:
            t=s['tau'];v=free_circuit(L,0 if c['flux']=='0' else np.pi,c['initial_occupations'],t,c['steps'])
            z=np.array(v['local_z']).reshape(L,L,2);d=np.array(v['local_d']).reshape(L,L)
            zerr.append(float(np.max(abs(((z+z.transpose(1,0,2))/2).ravel()-pub['exact',t]['local_z']))))
            derr.append(float(np.max(abs(((d+d.T)/2).ravel()-pub['exact',t]['local_d']))))
    assert len(zerr)==288 and max(zerr)<1e-12 and max(derr)<1e-7
    assert max(zerr)==audit['maximum_symmetrized_Z_difference']
    assert max(derr)==audit['maximum_symmetrized_D_difference']
    validation.update(status='PASS: complete matched-reference census, exact-state controls, direct/incremental comparison, U=0 ordering replay, and signed error accounting',
        free_circuit_reference_points=len(zerr),maximum_free_circuit_Z_difference=max(zerr),
        scope='4x4 U=4 zero-flux center-hole experiment only. Observed discrepancies include sampling and processing; no ensemble bias certificate or 6x6 extrapolation.')
    result=dict(metrics=metrics,curves=curves,decomposition=decomposition,mps_comparison=reference_comparison)
    inputs={str(p.relative_to(ROOT)):sha(p) for p in paths+[OUT/'build.json',OUT/'circuit_audit.json',
        OUT/'smoke_accounting.json',OUT/'production_accounting.json',
        ROOT/'python/willow_reference_analysis.py',ROOT/'python/willow_benchmark.py',ROOT/'refs/2510.26845v1.pdf']}
    return result,validation,inputs


def run(check=False):
    result,validation,inputs=analyze()
    if check:
        assert read(OUT/'results.json')==result
        assert read(OUT/'validation.json')==validation
        assert read(OUT/'inputs.json')==inputs
    else:
        save(OUT/'results.json',result);save(OUT/'validation.json',validation);save(OUT/'inputs.json',inputs)
        for name in ('metrics','decomposition','mps_comparison'):
            with (OUT/f'{name}.csv').open('w',newline='') as out:
                writer=csv.DictWriter(out,fieldnames=list(result[name][0]));writer.writeheader();writer.writerows(result[name])
    print(json.dumps(validation,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true')
    run(p.parse_args().check)
