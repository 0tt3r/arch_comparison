"""Exact rotation-fault submodel, noisy PCEC controls, and explicit RUS supply.

The projected multirotation channel is exact for independent generator-Pauli
faults and perfect syndrome projection. Other extraction/preparation faults and
growth remain separately declared sensitivities or measured Clifford probes.
"""
from functools import lru_cache
from math import atan, ceil, comb, cos, exp, expm1, log, log1p, pi, sin, tan
import hashlib
import json

import numpy as np

from hubbard_circuits import build_workload
from project_paths import ROOT

OUT=ROOT/'results/implementation_v1'


def wrap(theta):
    angle=(theta+pi/8)%(pi/4)-pi/8
    return 0. if abs(angle)<1e-13 else angle


@lru_cache(maxsize=None)
def preparation(theta,d,p):
    theta=wrap(theta)
    if theta==0:return dict(acceptance=1.,multiplier=1+0j,components=((1.,0.),),alpha=0.,infidelity=0.)
    k=(d+1)//2;alpha=atan(abs(tan(theta))**(1/k));c,s=cos(alpha),sin(alpha)
    # Pair-generator ZZ faults have probability p/15; the remaining one-qubit
    # Z fault has probability (p/10)/3. All other fault locations are excluded.
    q2,q1=p/15,p/30;terms=[];acceptance=0.
    for pairs in range(k):
        for single in (0,1):
            w=pairs+single
            fault=comb(k-1,pairs)*q2**pairs*(1-q2)**(k-1-pairs)*q1**single*(1-q1)**(1-single)
            passed=c**(2*(k-w))*s**(2*w)+c**(2*w)*s**(2*(k-w))
            phi=atan((-1)**w*tan(alpha)**(k-2*w))*(1 if theta>0 else -1)
            terms.append((fault*passed,phi-theta));acceptance+=fault*passed
    components=tuple((weight/acceptance,delta) for weight,delta in terms)
    multiplier=sum(weight*np.exp(-2j*delta) for weight,delta in components)
    return dict(acceptance=acceptance,multiplier=multiplier,components=components,alpha=alpha,infidelity=(1-multiplier.real)/2)


def growth_error(d,p,floor,data_distance):
    # `floor` is anchored at d=9, p=1e-3. Scaling the whole distance law
    # prevents a held d=9 anchor from making lower p worse at smaller d.
    return 0. if d==data_distance or floor==0 or p==0 else min(.25,floor*(p/.001)**5*(d/9)**2*(.01/p)**(-(d-9)/2))


def trial_rounds(theta,d,p,data_distance,detector_factor):
    prep=preparation(theta,d,p)
    acceptance=prep['acceptance']*(1-min(.99,detector_factor*p))**(d*d+6*d)
    rounds=(4+d*((d+1)//2%2==0))/acceptance
    if d<data_distance:rounds+=data_distance+1
    return rounds+9*data_distance


@lru_cache(maxsize=None)
def raw_rus_details(theta,d,p,floor=0.,data_distance=25,detector_factor=1.,cap=32):
    failed=1+0j;total=0j;rounds=0.;attempts=0.;tail=2**(-cap)
    for n in range(cap):
        angle=wrap(theta*2**n)
        if angle==0:
            total+=2**(-n)*failed;tail=0.;break
        c=preparation(angle,d,p)['multiplier']*(1-2*growth_error(d,p,floor,data_distance))
        total+=2**(-n-1)*failed*c;failed*=np.conj(c)
        rounds+=2**(-n)*trial_rounds(angle,d,p,data_distance,detector_factor)
        attempts+=2**(-n)
    return total/(1-tail),rounds,attempts,tail


def raw_rus(theta,d,p,cap=32):
    return raw_rus_details(theta,d,p,cap=cap)[0]


def corrected_trial(theta,d,p,mode,floor=0.,data_distance=25,detector_factor=1.):
    v=preparation(theta,d,p);c=v['multiplier']
    if mode=='raw':return c,0.,0.,0.
    if mode not in ('ideal_pcec','noisy_pcec','switch'):raise ValueError(mode)
    correction=0j;rounds=0.;attempts=0.;tail=0.
    for weight,delta in v['components']:
        if abs(delta)<1e-12:
            correction+=weight;continue
        channel,time,count,truncated=raw_rus_details(-delta,d,p,floor,data_distance,detector_factor)
        correction+=weight*np.exp(2j*delta)*(1 if mode=='ideal_pcec' else channel)
        rounds+=weight*time;attempts+=weight*count;tail+=weight*truncated
    answer=c*correction
    if mode=='switch' and (1-answer.real*(1-2*growth_error(d,p,floor,data_distance)))/2>2*p/15:
        # This original-injection reference includes growth in its error rate.
        return complex(1-4*p/15),-1.,0.,0.
    return answer,rounds,attempts,tail


def rus_channel(theta,d,p,mode='noisy_pcec',floor=0.,data_distance=25,detector_factor=1.,cap=24):
    theta=wrap(theta);failed=1+0j;total=0j;seconds_rounds=0.;switched=0.;extra=0.
    tail=2**(-cap);correction_tail=0.
    for n in range(cap):
        angle=wrap(theta*2**n)
        if angle==0:
            total+=2**(-n)*failed;tail=0.;break
        c,correction_rounds,corrections,ctr_tail=corrected_trial(angle,d,p,mode,floor,data_distance,detector_factor)
        if correction_rounds<0:
            switched+=2**(-n);rounds=4/.9+data_distance+9*data_distance
        else:
            rounds=trial_rounds(angle,d,p,data_distance,detector_factor)+correction_rounds
            c*=1-2*growth_error(d,p,floor,data_distance)
            extra+=2**(-n)*corrections
            correction_tail+=2**(-n)*ctr_tail
        total+=2**(-n-1)*failed*c;failed*=np.conj(c)
        seconds_rounds+=2**(-n)*rounds
    total/=1-tail
    return dict(q_dephasing=(1-total.real)/2,coherent_component=float(total.imag),
                multiplier=[float(total.real),float(total.imag)],
                mean_rounds=seconds_rounds,mean_seconds=seconds_rounds*240e-9,
                original_protocol_attempts=switched,extra_correction_attempts=extra,
                uncorrected_rus_tail=tail+correction_tail,injection_distance=d,data_distance=data_distance)


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    config=json.loads((ROOT/'configs/implementation_checks.json').read_text())['star'];rows=[]
    for p in config['p_values']:
        for theta in np.geomspace(1e-4,pi/8,41):
            for d in (3,5,9,15,25):
                prep=preparation(float(theta),d,p)
                for mode in ('raw','ideal_pcec','noisy_pcec','switch'):
                    result=rus_channel(float(theta),d,p,mode)
                    rows.append(dict(theta=float(theta),p=p,d=d,mode=mode,acceptance=prep['acceptance'],
                                     physical_angle=prep['alpha'],state_infidelity=prep['infidelity'],**result))
    # Actual angle histograms from the selected strict/loose L=4 STAR workloads.
    old=json.loads((ROOT/'results/physics_sweep_v1/costs.json').read_text());campaigns=[]
    for eps in (.01,.001):
        baseline=next(v for v in old if v['arm']=='star' and v['L']==4 and v['p']==.001 and v['absolute_error']==eps and v['scenario']=='held')
        ops,_=build_workload(4,1,baseline['r'],'termwise');histogram={}
        for op in ops:
            if op.name=='rotate':histogram[abs(op.angle)]=histogram.get(abs(op.angle),0)+1
        for p in config['p_values']:
            for floor in config['growth_floor_scenarios']:
                for detector_factor in (1.,4.):
                    for d in config['injection_distances']:
                        for mode in ('noisy_pcec','switch'):
                            cap=max(config['rus_cap'],ceil(log(sum(histogram.values())/(eps*.01),2)))
                            cases=[(count,rus_channel(theta,d,p,mode,floor,config['data_distance'],detector_factor,cap)) for theta,count in histogram.items()]
                            # Inverse-norm dephasing term and uncompensated coherent
                            # channel are kept separate; no pretend Pauli certificate.
                            chi=sum(-2*count*log(max(1e-300,1-2*v['q_dephasing'])) for count,v in cases)
                            # The dephasing inverse also rescales the coherent
                            # residue. Its non-CPTP product requires a norm
                            # product bound, rather than just adding raw biases.
                            coherent_log_norm=sum(count*log1p(abs(v['coherent_component'])/max(1e-300,1-2*v['q_dephasing'])) for count,v in cases)
                            phase_bias_bound=.5*expm1(coherent_log_norm)
                            seconds=sum(count*v['mean_seconds'] for count,v in cases)
                            campaigns.append(dict(epsilon=eps,p=p,r=baseline['r'],d=d,mode=mode,floor_at_d9=floor,growth_floor_anchor_p=.001,
                                detector_factor=detector_factor,rotations=sum(histogram.values()),chi=chi,
                                residual_coherent_bias_bound=phase_bias_bound,rus_tail_bound=sum(count*v['uncorrected_rus_tail'] for count,v in cases),
                                rus_cap=cap,serial_source_and_teleport_seconds=seconds,
                                injection_patch_qubits=2*d*d,data_patch_qubits=2*config['data_distance']**2,
                                log_sampling_time=log(seconds)+chi,
                                scope='Fixed workload; generator-fault channel and declared detector/growth scenarios, not calibrated full injection'))
    record=dict(rows=rows,campaigns=campaigns,configuration=config,
        equations='Perfect projection with independent generator-Pauli faults; exact error-angle mixture. Ideal PCEC multiplies c by conjugate(c); noisy PCEC evaluates one correction level using raw RUS. Diagnostic tail cap 24; workload cap reserves 1 percent of epsilon, correction cap 32. Correction time and growth channel evaluated at each correction angle.',
        source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
                      ('python/star_channel_study.py','configs/implementation_checks.json','refs/PhysRevX.15.021057.pdf','refs/2603.13093v1.pdf')})
    (OUT/'star_channels.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print('STAR:',len(rows),'channel points;',len(campaigns),'workload scenarios')


if __name__=='__main__':run()
