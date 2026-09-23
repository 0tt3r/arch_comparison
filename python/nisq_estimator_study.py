"""Exact eight-qubit bias/variance calibration with held-out interacting tests.

Training uses U=0 circuits with the interaction gate slots retained. A fixed
affine map is fit to exact noisy training expectations; finite training-shot
uncertainty is then charged analytically. Test outputs never enter the fit.
"""
from dataclasses import asdict
import argparse
import hashlib
import json
from math import ceil, sqrt

import numpy as np
from numba import njit

from hubbard_circuits import Operation
from nisq_compilation import gaussian_rzz, schedule
from parallel_nisq import PairEmitter, plaquette_stages
from project_paths import ROOT
from small_dynamics import HubbardSector, rectangle

OUT=ROOT/'results/implementation_v1'


def settings():
    return json.loads((ROOT/'configs/implementation_checks.json').read_text())['nisq']


def native_circuit(tau,steps,u):
    # Routing is independent of U. Density slots survive at U=0 and carry noise.
    sites=rectangle(2,2);order=sorted(range(8),key=lambda m:(m%2,m))
    emitter=PairEmitter(order)
    for stage in plaquette_stages(2,tau,steps,u_over_j=1):
        groups=sorted(stage.groups,key=lambda group:min(emitter.position[m] for b in group for m in b.modes))
        target=[]
        for group in groups:
            target.extend(sorted({m for b in group for m in b.modes},key=emitter.position.get))
        emitter.regroup(target)
        for group in groups:
            for block in group:emitter.block(block)
    occupied={2*i+(x+y)%2 for i,(x,y) in enumerate(sites)}
    gates=[Operation('x',(q,)) for q,m in enumerate(order) if m in occupied]
    for node in emitter.nodes:
        gates.extend(gaussian_rzz(node['pair'],node['matrix']))
        if abs(node['density'])>1e-14:
            gates.append(Operation('rzz',node['pair'],2*u*node['density'],role='onsite'))
    return gates,emitter.position


@njit(cache=True)
def one_qubit(rho,q,U):
    n=rho.shape[0];bit=1<<q
    for i in range(n):
        if i&bit:continue
        j=i^bit
        for k in range(n):
            a,b=rho[i,k],rho[j,k]
            rho[i,k]=U[0,0]*a+U[0,1]*b;rho[j,k]=U[1,0]*a+U[1,1]*b
    for i in range(n):
        if i&bit:continue
        j=i^bit
        for k in range(n):
            a,b=rho[k,i],rho[k,j]
            rho[k,i]=np.conj(U[0,0])*a+np.conj(U[0,1])*b
            rho[k,j]=np.conj(U[1,0])*a+np.conj(U[1,1])*b


@njit(cache=True)
def diagonal_gate(rho,mask,angle):
    n=rho.shape[0]
    for i in range(n):
        a=i&mask;si=1
        while a:
            si=-si;a&=a-1
        for j in range(n):
            a=j&mask;sj=1
            while a:
                sj=-sj;a&=a-1
            rho[i,j]*=np.exp(-.5j*angle*(si-sj))


@njit(cache=True)
def depolarize(rho,masks,p,scratch):
    if p==0:return
    n=rho.shape[0];support=0
    for m in masks:support|=m
    dim=len(masks);strength=p*dim*dim/(dim*dim-1)
    for i in range(n):
        for j in range(n):
            v=0j
            if (i^j)&support==0:
                for m in masks:v+=rho[i^m,j^m]
                v/=dim
            scratch[i,j]=(1-strength)*rho[i,j]+strength*v
    rho[:,:]=scratch


def probabilities(gates,p,noise_multiplier=1.,overrotation=0.):
    c=settings();n=8;rho=np.zeros((256,256),complex);rho[0,0]=1;scratch=np.empty_like(rho)
    sched,events=schedule(gates,n,gate_ns=c['gate_ns'],reset_ns=c['reset_ns'],measurement_ns=c['measurement_ns'])
    layers={}
    for start,gate in events:layers.setdefault(start,[]).append(gate)
    # Same p SPAM convention as the reference density calculation.
    for q in range(n):
        flip=np.arange(256)^(1<<q);rho=(1-p)*rho+p*rho[flip[:,None],flip[None,:]]
    for layer in range(sched['depth_layers']):
        busy=set()
        for gate in layers.get(layer,[]):
            qs=gate.qubits;busy.update(qs);angle=gate.angle*(1+overrotation if gate.role=='onsite' else 1)
            if gate.name in ('rz','rzz'):
                diagonal_gate(rho,sum(1<<q for q in qs),angle)
            else:
                if gate.name=='x':U=np.array([[0,1],[1,0]],complex)
                elif gate.name=='h':U=np.array([[1,1],[1,-1]],complex)/sqrt(2)
                elif gate.name=='rx':U=np.array([[np.cos(angle/2),-1j*np.sin(angle/2)],[-1j*np.sin(angle/2),np.cos(angle/2)]])
                else:raise ValueError(gate.name)
                one_qubit(rho,qs[0],U)
            masks=np.array([sum((1<<q) for k,q in enumerate(qs) if (s>>k)&1) for s in range(2**len(qs))])
            depolarize(rho,masks,p*noise_multiplier if len(qs)==2 else p/10,scratch)
        for q in set(range(n))-busy:
            depolarize(rho,np.array([0,1<<q]),p*c['gate_ns']/100,scratch)
    probs=rho.diagonal().real.copy()
    for q in range(n):probs=(1-p)*probs+p*probs[np.arange(256)^(1<<q)]
    assert abs(probs.sum()-1)<1e-10 and probs.min()>-1e-12
    return probs,sched


def observables(probs,mapping):
    bits=[(np.arange(256)>>q)&1 for q in range(8)]
    up=sum(bits[mapping[m]] for m in range(0,8,2));down=sum(bits[mapping[m]] for m in range(1,8,2))
    i=rectangle(2,2).index((1,1));D=bits[mapping[2*i]]*bits[mapping[2*i+1]]
    out={}
    for key,mask in [('raw',np.ones(256,bool)),('number',(up==2)&(down==2))]:
        acceptance=float(probs[mask].sum());mean=float(np.dot(probs*mask,D)/acceptance)
        out[key]=dict(mean=mean,acceptance=acceptance,variance_per_attempt=mean*(1-mean)/acceptance)
    return out


def affine_fit(training,method):
    x=np.array([v['noisy'][method]['mean'] for v in training]);y=np.array([v['ideal'] for v in training])
    X=np.column_stack((np.ones(len(x)),x));operator=np.linalg.pinv(X);beta=operator@y
    return beta,operator


def evaluate_estimator(test,training,method):
    c=settings();raw_method='number' if method in ('number','tflo_number') else 'raw'
    v=test['noisy'][raw_method];x=v['mean'];training_var=0.
    if method.startswith('tflo'):
        beta,operator=affine_fit(training,raw_method);a,b=beta
        # Delta-method covariance of fitted OLS coefficients for noisy x_i.
        # Includes the residual term, rather than treating labels as a perfect line.
        X=np.array([[1,t['noisy'][raw_method]['mean']] for t in training])
        y=np.array([t['ideal'] for t in training]);residual=y-X@beta
        inv=np.linalg.inv(X.T@X);h=np.array([1.,x]);gradient=[]
        for i in range(len(training)):
            gradient.append(h@inv@(np.array([0.,residual[i]])-b*X[i]))
        training_var=sum(g*g*t['noisy'][raw_method]['variance_per_attempt']/c['training_shots_per_circuit']
                         for g,t in zip(gradient,training))
    else:a,b=0.,1.
    prediction=a+b*x;bias=prediction-test['ideal'];variance=b*b*v['variance_per_attempt']
    allowance=(c['statistical_allowance']/c['confidence_z'])**2-training_var
    shots=ceil(variance/allowance) if allowance>0 else None
    train_attempts=len(training)*c['training_shots_per_circuit'] if method.startswith('tflo') else 0
    train_seconds=sum(t['attempt_seconds']*c['training_shots_per_circuit'] for t in training) if train_attempts else 0.
    return dict(method=method,prediction=prediction,bias_from_ideal_circuit=bias,
        bias_from_exact_dynamics=prediction-test['exact'],variance_per_attempt=variance,
        acceptance=v['acceptance'],fit_intercept=float(a),fit_slope=float(b),training_variance=training_var,
        training_attempts=train_attempts,measurement_attempts=shots,
        seconds_including_training=None if shots is None else shots*test['attempt_seconds']+train_seconds,
        passes_bias_allowance=bool(abs(bias)<=c['bias_allowance']),
        passes_total_error=bool(shots is not None and abs(bias)+abs(test['ideal']-test['exact'])+c['statistical_allowance']<=c['total_error']),
        scope='Exact simulated channel bias; Gaussian/delta-method shot planning, not a hardware or finite-size certificate')


def run():
    OUT.mkdir(parents=True,exist_ok=True);c=settings();path=OUT/'nisq.json'
    sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
             [ROOT/'python/nisq_estimator_study.py',ROOT/'configs/implementation_checks.json']}
    record=json.loads(path.read_text()) if path.exists() else dict(cases=[],source_sha256=sources)
    assert record['source_sha256']==sources,'Use a fresh result file after changing the experiment'
    done={v['key']:v for v in record['cases']}
    for p in c['p_values']:
        for r in c['steps']:
            specifications=[('train',t,0.,'matched') for t in c['training_times']]
            specifications += [('test',t,u,noise) for t in c['test_times'] for u in c['test_interactions'] for noise in ('matched','mismatch')]
            for role,tau,u,noise in specifications:
                key=f'{p}_{r}_{role}_{tau}_{u}_{noise}'
                if key in done:continue
                gates,mapping=native_circuit(tau,r,u)
                changes={} if noise=='matched' else dict(noise_multiplier=c['mismatch']['two_qubit_noise_multiplier'],overrotation=c['mismatch']['onsite_fractional_overrotation'])
                probs,sched=probabilities(gates,p,**changes)
                model=HubbardSector(rectangle(2,2),u=u)
                ideal=model.doublon(model.trotter(tau,r),(1,1));exact=model.doublon(model.exact(tau),(1,1))
                row=dict(key=key,role=role,tau=tau,u=u,p=p,r=r,noise=noise,ideal=ideal,exact=exact,
                    noisy=observables(probs,mapping),gate_count=len(gates),two_qubit_gates=sum(len(g.qubits)==2 for g in gates),
                    depth=sched['depth_layers'],attempt_seconds=(c['reset_ns']+c['measurement_ns']+sched['depth_layers']*c['gate_ns'])*1e-9,
                    circuit_sha256=hashlib.sha256(json.dumps([asdict(g) for g in gates],sort_keys=True).encode()).hexdigest())
                record['cases'].append(row);done[key]=row
                temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(record,indent=2)+'\n');temporary.replace(path)
                print('NISQ',key,'bias',row['noisy']['number']['mean']-ideal,flush=True)
    estimates=[]
    for case in record['cases']:
        if case['role']!='test':continue
        training=[v for v in record['cases'] if v['role']=='train' and v['p']==case['p'] and v['r']==case['r']]
        for method in ('raw','number','tflo_raw','tflo_number'):
            estimates.append({k:case[k] for k in ('key','p','r','tau','u','noise','ideal','exact')}|evaluate_estimator(case,training,method))
    record['estimates']=estimates;record['configuration']=c
    record['scope']='2x2 Neel quench only; exact physical channel, independent U=0 training, held-out U=2,4; zero-angle interaction slots retain noise'
    path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':run()
