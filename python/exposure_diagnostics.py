"""Test exposure hypotheses against saved channel and estimator measurements.

Neither ordinate is manufactured by evaluating the assumed resource exponent.
The NISQ variance and bias come from exact noisy density matrices and held-out
fits. STAR points come from generator-fault instrument enumeration. Scope and
failed bias tests are retained in the output.
"""
import hashlib
import json
from math import log

import numpy as np
from scipy.stats import beta

from model import log_gamma2_dephasing,log_gamma2_depolarizing
from project_paths import ROOT

OUT=ROOT/'results/frontier_v1'


def nisq_exposure(case,configuration):
    p=case['p'];g2=case['two_qubit_gates'];g1=case['gate_count']-g2
    idle=8*case['depth']-2*g2-g1
    p2=p*(configuration['mismatch']['two_qubit_noise_multiplier'] if case['noise']=='mismatch' else 1)
    parts=dict(two_qubit=g2*log_gamma2_depolarizing(p2,2),
        one_qubit=g1*log_gamma2_depolarizing(p/10,1),
        idle=idle*log_gamma2_depolarizing(p*configuration['gate_ns']/100,1),
        spam=16*log_gamma2_dephasing(p))
    return dict(exposure=sum(parts.values()),components=parts,idle_locations=idle,
        exposure_scope='Gatewise inverse-norm Pauli reference; coherent mismatch is outside this inverse')


def run():
    nisqpath=ROOT/'results/implementation_v1/nisq.json';starpath=ROOT/'results/implementation_v1/star_channels.json'
    data=json.loads(nisqpath.read_text());star=json.loads(starpath.read_text());cases={v['key']:v for v in data['cases']};rows=[]
    for v in data['estimates']:
        case=cases[v['key']];reference=case['ideal']*(1-case['ideal'])
        F=v['variance_per_attempt']/reference
        idealshots=reference*(data['configuration']['confidence_z']/data['configuration']['statistical_allowance'])**2
        total=v['seconds_including_training']/(idealshots*case['attempt_seconds']) if v['seconds_including_training'] is not None else None
        rows.append(dict(**v,**nisq_exposure(case,data['configuration']),measured_variance_multiplier=F,
            log_measured_variance_multiplier=log(F),total_time_multiplier_with_training=total))
    channels=[dict(**v,normalized_error=v['q_dephasing']/(v['p']*v['theta'])) for v in star['rows']]
    groups=[]
    for p in (.001,.0001):
      for method in ('raw','number','tflo_raw','tflo_number'):
       for noise in ('matched','mismatch'):
        points=[v for v in rows if (v['p'],v['method'],v['noise'])==(p,method,noise)]
        groups.append(dict(p=p,method=method,noise=noise,cases=len(points),bias_pass=sum(v['passes_bias_allowance'] for v in points),
            total_error_pass=sum(v['passes_total_error'] for v in points),
            variance_multiplier_min=min(v['measured_variance_multiplier'] for v in points),
            variance_multiplier_max=max(v['measured_variance_multiplier'] for v in points)))
    names=['python/exposure_diagnostics.py','results/implementation_v1/nisq.json','results/implementation_v1/star_channels.json']
    (OUT/'exposure.json').write_text(json.dumps(dict(nisq=rows,star=channels,summary=groups,
        scope='A test of collapse, not a fitted universal mitigation exponent. Density-matrix data are eight qubits; STAR data exclude generic syndrome/preparation faults.',
        source_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')


def growth_summary():
    files=sorted((OUT/'growth/cases').glob('d*.json'));cases=[json.loads(p.read_text()) for p in files]
    assert len(cases)==24 and sum(v['shots'] for v in cases)==24000000
    for v in cases:
        assert v['status']=='complete' and v['pbs_job']=='7790943.imgt1'
        k,n=v['failures'],v['shots'];upper=float(beta.ppf(1-.05/24,k+1,n-k)) if k<n else 1.
        assert abs(upper-v['simultaneous_upper95'])<1e-15
        # Keys contain decimal points, so append rather than replace a suffix.
        assert hashlib.sha256((OUT/'growth/cases'/(v['contract']['key']+'.stim')).read_bytes()).hexdigest()==v['circuit_sha256']
        for name,digest in v['contract']['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    accounting=(OUT/'growth_production.qstat.txt').read_text()
    assert 'job_state = F' in accounting and 'Exit_status = 0' in accounting
    pairs=[]
    for di,df in ((5,13),(9,13),(13,13),(5,25),(9,25),(13,25)):
      for p in (.001,.0001):
        selected=[v for v in cases if (v['contract']['di'],v['contract']['df'],v['contract']['p'])==(di,df,p)]
        assert len(selected)==2
        q=min(.5,sum(v['simultaneous_upper95'] for v in selected))
        pairs.append(dict(di=di,df=df,p=p,failures={v['contract']['basis']:v['failures'] for v in selected},
            pauli_total_upper=q,chi_per_growth_upper=log_gamma2_dephasing(q),
            scope='Union of X/Z campaign upper bounds, conditional on a common Pauli logical-channel interpretation; ideal encoded input, target memory and readout included, not injected-state error'))
    (OUT/'growth_summary.json').write_text(json.dumps(dict(cases=cases,pairs=pairs,shots=24000000,
        status='PASS: 24 complete cells, source/circuit hashes, simultaneous intervals, PBS F/0',
        accounting_sha256=hashlib.sha256((OUT/'growth_production.qstat.txt').read_bytes()).hexdigest()),indent=2)+'\n')


if __name__=='__main__':
    run();growth_summary()
