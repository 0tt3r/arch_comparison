"""Single-process C6/C7 surface-memory calibration; all output is isolated.

Commands: validate, sample, fit, audit. See results/surface_law_v1/manifest.json.
No multiprocessing, sinter worker pool, or implicit BLAS threading is used.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import socket
import sys
import time

import numpy as np
import pymatching
import scipy
from scipy.optimize import minimize, brentq
from scipy.special import xlogy
from scipy.stats import beta, chi2
import stim

ROOT = Path('/lcrc/project/QC_NA_Sim/architecture_comparison')
OUT = ROOT / 'results/surface_law_v1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def circuit(distance, rounds, p, basis):
    """Compile Stim's hook-safe rotated schedule into dressed CZ macro-layers.

    Every CNOT is H(target)-CZ-H(target), each H has p/10 noise.
    The generator's ancilla H gates are also noisy. Basis changes are fused
    into interaction slots: no additional idle macro-layer or time is assigned.
    Four interaction slots, M slot, R slot; idle p per unused slot. Initial
    product prep and final data measurement carry flip p. MR is split M/R.
    """
    template = stim.Circuit.generated('surface_code:rotated_memory_' + basis,
                                      distance=distance, rounds=rounds)
    coordinates = template.get_final_qubit_coordinates()
    all_qubits = set(coordinates)
    data = {q for q, (x, y) in coordinates.items() if x % 2 and y % 2}
    out = stim.Circuit()

    def noise(gate, targets, probability):
        if p and targets:
            out.append(gate, targets, probability)

    for ins in template.flattened():
        name = ins.name
        ts = [t.value for t in ins.targets_copy() if t.is_qubit_target]
        if name == 'TICK':
            continue  # Ticks below mark precisely the six declared macro-slots.
        if name == 'CX':
            targets = ts[1::2]
            out.append('H', targets)
            noise('DEPOLARIZE1', targets, p / 10)
            out.append('CZ', ts)
            noise('DEPOLARIZE2', ts, p)
            out.append('H', targets)
            noise('DEPOLARIZE1', targets, p / 10)
            noise('DEPOLARIZE1', sorted(all_qubits - set(ts)), p)
            out.append('TICK')
        elif name == 'H':
            out.append(ins)
            noise('DEPOLARIZE1', ts, p / 10)
        elif name in ('R', 'RX'):
            out.append(ins)
            noise('Z_ERROR' if name == 'RX' else 'X_ERROR', ts, p)
        elif name == 'MR':
            noise('X_ERROR', ts, p)
            out.append('M', ts)
            noise('DEPOLARIZE1', sorted(data), p)
            out.append('TICK')
            out.append('R', ts)
            noise('X_ERROR', ts, p)
            noise('DEPOLARIZE1', sorted(data), p)
            out.append('TICK')
        elif name in ('M', 'MX'):
            noise('Z_ERROR' if name == 'MX' else 'X_ERROR', ts, p)
            out.append(ins)
        else:
            out.append(ins)
    return out


def interval(k, n, alpha=.05):
    return [0. if k == 0 else float(beta.ppf(alpha/2, k, n-k+1)),
            1. if k == n else float(beta.ppf(1-alpha/2, k+1, n-k))]


def per_round(q, rounds):
    return float(-.5 * np.expm1(np.log1p(-2 * min(q, .499999999)) / rounds))


def validate():
    rows = []
    for d in (3, 5, 7, 9):
        for basis in ('x', 'z'):
            c = circuit(d, 2*d, .001, basis)
            ideal = circuit(d, 2*d, 0, basis)
            det, obs = ideal.compile_detector_sampler(seed=17).sample(256, separate_observables=True)
            assert not det.any() and not obs.any()
            dem = c.detector_error_model(decompose_errors=True)
            decoder = pymatching.Matching.from_detector_error_model(dem)
            coords = c.get_final_qubit_coordinates()
            data = [q for q, (x, y) in coords.items() if x % 2 and y % 2]
            logical = [q for q in data if coords[q][1 if basis == 'x' else 0] == 1]
            control = stim.Circuit()
            inserted = False
            for ins in ideal:
                if ins.name == ('MX' if basis == 'x' else 'M') and set(t.value for t in ins.targets_copy()) == set(data):
                    control.append('Z_ERROR' if basis == 'x' else 'X_ERROR', logical, 1)
                    inserted = True
                control.append(ins)
            det, obs = control.compile_detector_sampler(seed=18).sample(256, separate_observables=True)
            assert inserted and not det.any() and obs.all()
            assert np.all(decoder.decode_batch(det) != obs)
            # Both representations of batch decoding must agree bit-for-bit.
            det, obs = c.compile_detector_sampler(seed=19).sample(2048, separate_observables=True)
            packed = np.packbits(det, axis=1, bitorder='little')
            assert np.array_equal(decoder.decode_batch(det), decoder.decode_batch(packed, bit_packed_shots=True))
            assert sum(ins.name == 'TICK' for ins in c) == 12*d
            assert sum(len(ins.targets_copy()) // 2 for ins in c if ins.name == 'CZ') == 2*d * (4*d*(d-1))
            distance_found = len(c.shortest_graphlike_error())
            assert distance_found == d, (d, basis, distance_found)
            rows.append(dict(distance=d, basis=basis, graphlike_distance=distance_found,
                             noiseless_shots=256, logical_flip_shots=256,
                             packed_decode_equal=True, macro_layers_per_round=6))
    result = dict(status='pass', rows=rows, stim=stim.__version__, pymatching=pymatching.__version__)
    atomic(OUT/'validation.json', result)
    print(json.dumps(result), flush=True)


def sample(args):
    folder = OUT / args.group
    folder.mkdir(parents=True, exist_ok=True)
    key = f'd{args.distance}_{args.basis}_p{args.p:g}_r{args.rounds}_s{args.seed}'
    path = folder / (key + '.json')
    contract = dict(distance=args.distance, basis=args.basis, p=args.p, rounds=args.rounds,
                    target_shots=args.shots, seed=args.seed, batch=args.batch,
                    program_sha256=sha(__file__), decoder='pymatching MWPM; correlations disabled')
    if path.exists():
        previous = json.loads(path.read_text())
        assert previous['contract'] == contract and previous['status'] == 'complete'
        print('Already complete:', key, flush=True)
        return
    lock = folder / (key + '.lock')
    with lock.open('x') as f:
        json.dump(dict(pid=os.getpid(), host=socket.gethostname(), job=os.environ.get('PBS_JOBID')), f)
    started = time.monotonic()
    try:
        c = circuit(args.distance, args.rounds, args.p, args.basis)
        cpath = folder / (key + '.stim')
        cpath.write_text(str(c) + '\n')
        decoder = pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
        sampler = c.compile_detector_sampler(seed=args.seed)
        shots = failures = 0
        batch_records = []
        while shots < args.shots:
            n = min(args.batch, args.shots-shots)
            det, obs = sampler.sample(n, separate_observables=True, bit_packed=True)
            predicted = decoder.decode_batch(det, bit_packed_shots=True, bit_packed_predictions=True)
            k = int(np.count_nonzero(np.any(predicted != obs, axis=1)))
            failures += k
            shots += n
            batch_records.append([n, k])
            elapsed = time.monotonic()-started
            if shots == args.shots or len(batch_records) % 10 == 0:
                result = dict(status='complete' if shots == args.shots else 'running', contract=contract,
                    shots=shots, failures=failures, interval95=interval(failures, shots),
                    finite_window_per_round=per_round(failures/shots, args.rounds),
                    finite_window_interval95=[per_round(q,args.rounds) for q in interval(failures,shots)],
                    batches=batch_records, elapsed_seconds=elapsed, shots_per_second=shots/elapsed,
                    circuit_sha256=sha(cpath), host=socket.gethostname(), pbs_job=os.environ.get('PBS_JOBID'),
                    stim=stim.__version__, pymatching=pymatching.__version__, numpy=np.__version__,
                    scipy=scipy.__version__, python=sys.version, machine=platform.machine(),
                    peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                atomic(path if shots == args.shots else path.with_suffix('.progress.json'), result)
                print(key, shots, failures, round(elapsed,2), flush=True)
        progress = path.with_suffix('.progress.json')
        if progress.exists():
            progress.unlink()
    finally:
        lock.unlink()


class LawFit:
    """Binomial likelihood with basis amplitudes and optional end corrections.

    q_b(d) = exp(a) * f_b * exp(-h*l), f_x+f_z=2.
    P_b(d,R) = [1-(1-2*q_b(d))**(R+delta_db)]/2.
    delta_db is a fitted effective-round correction, not a negative probability.
    The fit is to complete fixed-shot cells, including zero-failure cells.
    """
    def __init__(self, rows, endpoints=True):
        self.rows = rows
        self.d = np.array([r['contract']['distance'] for r in rows], float)
        self.h = (self.d+1)/2
        self.p = np.array([r['contract']['p'] for r in rows], float)
        self.reference_p = float(max(self.p))
        self.log_p_relative = np.log(self.p/self.reference_p)
        self.r = np.array([r['contract']['rounds'] for r in rows], float)
        self.sign = np.array([1 if r['contract']['basis']=='x' else -1 for r in rows])
        self.n = np.array([r['shots'] for r in rows], float)
        self.k = np.array([r['failures'] for r in rows], float)
        self.saturated = float(-np.sum(xlogy(self.k,self.k/self.n)+xlogy(self.n-self.k,1-self.k/self.n)))
        labels = [(float(p),int(d),int(s)) for p,d,s in zip(self.p,self.d,self.sign)]
        self.labels = sorted(set(labels)) if endpoints else []
        self.index = np.array([self.labels.index(v) for v in labels]) if endpoints else None
        self.bounds = [(-20.,5.),(-1.,10.),(-3.,3.)] + [(-9.5,40.)]*len(self.labels)
        self.x0 = np.array([-3.2, math.log(.007 / self.reference_p), 0.] + [0.]*len(self.labels))
        self.opt = self.optimize()

    def objective(self, x):
        a, l, b = x[:3]
        logistic = 1/(1+np.exp(-self.sign*b))
        logq = a-self.h*(l-self.log_p_relative)+np.log(2*logistic)
        if np.any(logq >= math.log(.49)):
            excess = np.maximum(0,logq-math.log(.49))
            grad = np.zeros(len(x))
            grad[0] = 2e10*np.sum(excess)
            grad[1] = -2e10*np.sum(self.h*excess)
            grad[2] = 2e10*np.sum(excess*self.sign*(1-logistic))
            return 1e12+1e10*np.sum(excess**2),grad
        q = np.exp(logq)
        t = self.r.copy()
        if self.labels:
            t += self.d*x[3:][self.index]
        log_survival = np.log1p(-2*q)
        survival = np.exp(t*log_survival)
        prob = np.maximum(-.5*np.expm1(t*log_survival),1e-300)
        # Remove the constant saturated likelihood to make optimizer tolerance
        # depend on lack of fit, not on an irrelevant O(millions) offset.
        value = np.sum(xlogy(self.k,self.k/(self.n*prob)) +
                       xlogy(self.n-self.k,(self.n-self.k)/(self.n*(1-prob))))
        score = (self.n*prob-self.k)/(prob*(1-prob))
        dq = t*q*survival/(1-2*q)
        grad = np.zeros(len(x))
        grad[0] = np.sum(score*dq)
        grad[1] = -np.sum(score*dq*self.h)
        grad[2] = np.sum(score*dq*self.sign*(1-logistic))
        if self.labels:
            np.add.at(grad[3:],self.index,score*(-.5*survival*log_survival)*self.d)
        return float(value),grad

    def optimize(self, fixed=None, start=None):
        x = self.x0.copy() if start is None else start.copy()
        free = list(range(len(x)))
        if fixed is not None:
            idx,value = fixed
            x[idx] = value
            free.remove(idx)
        def fun(y):
            z=x.copy(); z[free]=y
            value,gradient=self.objective(z)
            return value,gradient[free]
        result=minimize(fun,x[free],jac=True,method='L-BFGS-B',
                        bounds=[self.bounds[i] for i in free],
                        options=dict(maxiter=4000,ftol=1e-13,gtol=1e-6,maxls=50))
        if not result.success:
            # Very flat sparse endpoint directions can end an L-BFGS line
            # search at machine precision. Independently finish with SLSQP.
            fallback=minimize(fun,result.x,jac=True,method='SLSQP',
                bounds=[self.bounds[i] for i in free],options=dict(maxiter=1500,ftol=1e-9))
            if fallback.success and fallback.fun<=result.fun+1e-6:
                result=fallback
        x[free]=result.x
        return dict(x=x,fun=float(result.fun),success=bool(result.success),message=str(result.message))

    def profile_interval(self,index,level=.95):
        center=self.opt['x'][index]
        cutoff=chi2.ppf(level,1)
        def difference(value):
            opt=self.optimize((index,value),self.opt['x'])
            return 2*(opt['fun']-self.opt['fun'])-cutoff
        answers=[];bounded=[]
        for side in (-1,1):
            edge=center
            step=.025
            limit=self.bounds[index][0 if side<0 else 1]
            for _ in range(40):
                edge=max(limit,center-step) if side<0 else min(limit,center+step)
                if difference(edge)>0:
                    answers.append(float(brentq(difference,min(center,edge),max(center,edge),xtol=1e-6)))
                    bounded.append(True)
                    break
                if edge==limit:
                    answers.append(float(limit));bounded.append(False);break
                step*=1.65
        return dict(log_interval=answers,finite=bounded)

    def summary(self, intervals=True):
        a,l,b=self.opt['x'][:3]
        deviance=max(0.,2*self.opt['fun'])
        dof=len(self.rows)-len(self.opt['x'])
        result=dict(A_component_mean=float(np.exp(a)), A_patch_union_bound=float(2*np.exp(a)),
            Lambda=float(np.exp(l)), reference_p=self.reference_p,
            effective_p_threshold=float(self.reference_p*np.exp(l)),
            basis_amplitude_ratio_x_over_z=float(np.exp(b)),
            A_x=float(np.exp(a)*2/(1+np.exp(-b))),A_z=float(np.exp(a)*2/(1+np.exp(b))),
            endpoint_corrections_rounds={f'p{p:g}_d{d}_{"x" if s==1 else "z"}':float(d*v)
                for (p,d,s),v in zip(self.labels,self.opt['x'][3:])},
            endpoint_bounds_hit=any(abs(v+9.5)<.001 or abs(v-40)<.001 for v in self.opt['x'][3:]),
            negative_log_likelihood=self.opt['fun']+self.saturated,deviance=deviance,degrees_of_freedom=dof,
            goodness_of_fit_p_asymptotic=float(chi2.sf(deviance,dof)) if dof>0 else None,
            optimizer_success=self.opt['success'],optimizer_message=self.opt['message'],
            interval_warning='Profile likelihood uses asymptotic chi-square calibration; zero-event endpoint nuisance parameters can hit bounds, especially at p=1e-4.',
            cells=len(self.rows),shots=int(sum(self.n)),failures=int(sum(self.k)))
        if intervals:
            for idx,name in ((0,'A_component_mean'),(1,'Lambda')):
                for level,label in ((.95,'95'),(.975,'97_5')):
                    ci=self.profile_interval(idx,level)
                    result[name+'_profile'+label]=[float(np.exp(v)) for v in ci['log_interval']]
                    result[name+'_profile'+label+'_finite']=ci['finite']
            result['A_patch_union_bound_profile95']=[2*v for v in result['A_component_mean_profile95']]
            result['effective_p_threshold_profile95']=[self.reference_p*v for v in result['Lambda_profile95']]
        return result


def required_distance(a, suppression, volume, budget):
    if suppression <= 1:
        return None
    h=math.log(volume*a/budget)/math.log(suppression)
    # Round upward to an odd d, checking the defining inequality near boundaries.
    d=max(3,2*math.ceil(h-1e-12)-1)
    while math.log(volume*a/budget)-((d+1)/2)*math.log(suppression)>1e-11:
        d+=2
    return d


def load_rows(manifest_name='manifest.json', group='cases'):
    manifest=json.loads((OUT/manifest_name).read_text())
    rows=[]
    for cell in manifest['cases']:
        path=OUT/group/(cell['key']+'.json')
        value=json.loads(path.read_text())
        assert value['status']=='complete' and value['shots']==cell['shots']
        for key in ('distance','basis','p','rounds','seed','batch'):
            assert value['contract'][key]==cell[key]
        assert value['contract']['program_sha256']==manifest['worker_sha256']
        assert sha(path.with_suffix('.stim'))==value['circuit_sha256']
        assert sum(v[0] for v in value['batches'])==value['shots']
        assert sum(v[1] for v in value['batches'])==value['failures']
        assert not path.with_suffix('.lock').exists()
        rows.append(value)
    return manifest,rows


def holdout_analysis(training_fits):
    path=OUT/'holdout_manifest.json'
    if not path.exists():
        return dict(status='not requested')
    manifest=json.loads(path.read_text())
    if any(not (OUT/'holdout'/(c['key']+'.json')).exists() for c in manifest['cases']):
        return dict(status='running; excluded from calibration')
    _,rows=load_rows('holdout_manifest.json','holdout')
    fit_results={}
    for p in (.001,.0001):
        fit_results[str(p)]=LawFit([r for r in rows if r['contract']['p']==p],endpoints=False).summary()
        if p==.001:
            fit_results[str(p)]['distance_5_and_above_fit']=LawFit(
                [r for r in rows if r['contract']['p']==p and r['contract']['distance']>=5],
                endpoints=False).summary()
    cells=[]
    for r in rows:
        c=r['contract'];f=training_fits[str(c['p'])]
        q=f['A_'+c['basis']]*f['Lambda']**(-(c['distance']+1)/2)
        offset=f['endpoint_corrections_rounds'][f"p{c['p']:g}_d{c['distance']}_{c['basis']}"]
        predicted=-.5*math.expm1((c['rounds']+offset)*math.log1p(-2*q))
        cells.append({**{k:c[k] for k in ('p','distance','basis','rounds','seed')},
            **{k:r[k] for k in ('shots','failures','finite_window_per_round','finite_window_interval95')},
            'training_predicted_failures':r['shots']*predicted,
            'plug_in_residual_sigma':(r['failures']-r['shots']*predicted)/math.sqrt(r['shots']*predicted*(1-predicted))})
    return dict(status='complete',role='Independent 100d check; not pooled into training intervals.',
        total_shots=sum(r['shots'] for r in rows),total_failures=sum(r['failures'] for r in rows),
        fits=fit_results,cells=cells,
        caveat='100d fits set endpoint correction to zero. Plug-in residuals omit uncertainty in the training prediction; they are diagnostics, not calibrated predictive tests.')


def endpoint_slopes(rows):
    """Distance-by-distance estimates without imposing the exponential law.

    Four 98.75% Clopper-Pearson intervals give a >=95% simultaneous rectangle
    for both bases at one (p,d). Propagate it monotonically through the parity
    slope and sum. This covers statistical sampling, not nonstationarity.
    """
    output=[]
    for p,d in sorted({(r['contract']['p'],r['contract']['distance']) for r in rows}):
        result=dict(p=p,distance=d,components={})
        for basis in ('x','z'):
            pair=sorted([r for r in rows if (r['contract']['p'],r['contract']['distance'],r['contract']['basis'])==(p,d,basis)],key=lambda r:r['contract']['rounds'])
            assert len(pair)==2
            first,last=pair
            gap=last['contract']['rounds']-first['contract']['rounds']
            def slope(x,y):
                return -.5*math.expm1((math.log1p(-2*y)-math.log1p(-2*x))/gap)
            q=slope(first['failures']/first['shots'],last['failures']/last['shots'])
            a=interval(first['failures'],first['shots'],.0125)
            b=interval(last['failures'],last['shots'],.0125)
            result['components'][basis]=dict(rate=max(0.,q),unconstrained_slope=q,
                interval95_simultaneous_bases=[max(0.,slope(a[1],b[0])),max(0.,slope(a[0],b[1]))],
                fewer_than_30_failures=any(r['failures']<30 for r in pair))
        result['patch_component_sum']=sum(v['rate'] for v in result['components'].values())
        result['patch_component_sum_interval95']=[sum(v['interval95_simultaneous_bases'][i] for v in result['components'].values()) for i in (0,1)]
        result['interpretation']='Finite-length parity slope with common endpoint correction; no exponential distance law imposed. Intervals are simultaneous for the two bases at this p,d, not over the entire grid.'
        output.append(result)
    return output


def bootstrap_gof(repetitions=399):
    """Calibrate lack-of-fit diagnostics when sparse nuisance fits hit bounds."""
    _,rows=load_rows()
    rng=np.random.default_rng(2026092501)
    results={}
    for p in (.001,.0001):
        f=LawFit([r for r in rows if r['contract']['p']==p])
        x=f.opt['x'].copy();a,l,b=x[:3]
        q=np.exp(a-f.h*(l-f.log_p_relative))*2/(1+np.exp(-f.sign*b))
        t=f.r+f.d*x[3:][f.index]
        probabilities=-.5*np.expm1(t*np.log1p(-2*q))
        observed=2*f.opt['fun']
        deviances=[];failed=0
        for rep in range(repetitions):
            f.k=rng.binomial(f.n.astype(np.int64),probabilities).astype(float)
            opt=f.optimize(start=x)
            if not opt['success']:
                failed+=1
            deviances.append(max(0.,2*opt['fun']))
            if (rep+1)%50==0:
                print('bootstrap',p,rep+1,flush=True)
        exceed=sum(v>=observed for v in deviances)
        results[str(p)]=dict(repetitions=repetitions,observed_deviance=observed,
            bootstrap_deviances=deviances,tail_count=exceed,p_value=(1+exceed)/(repetitions+1),
            tail_probability_monte_carlo_interval95=interval(exceed,repetitions),
            optimizer_unsuccessful=failed)
    atomic(OUT/'gof_bootstrap.json',dict(seed=2026092501,results=results,
        manifest_sha256=sha(OUT/'manifest.json'),
        interpretation='Parametric binomial bootstrap under each fitted law; conditions on its nuisance estimates. It diagnoses lack of fit, not extrapolation reliability.'))


def fit():
    manifest,rows=load_rows()
    results={}
    comparisons=[]
    for p in (.001,.0001):
        cells=[r for r in rows if r['contract']['p']==p]
        print('Fitting',p,flush=True)
        primary=LawFit(cells,endpoints=True).summary()
        no_offset=LawFit(cells,endpoints=False).summary()
        lengths={}
        for mult in (10,20):
            selected=[r for r in cells if r['contract']['rounds']==mult*r['contract']['distance']]
            lengths[str(mult)+'d']=LawFit(selected,endpoints=False).summary(intervals=False)
        primary['length_sensitivity_fits']=lengths
        primary['zero_endpoint_fit']=no_offset
        if p==.001:
            primary['distance_5_and_above_fit']=LawFit([r for r in cells if r['contract']['distance']>=5],endpoints=True).summary()
        results[str(p)]=primary
        for volume in (1e8,1e11,1e14):
            for budget in (1e-4,1e-5,1e-6):
                baseline=required_distance(.1,.01/p,volume,budget)
                for convention,factor in (('component_mean',1),('patch_union_bound',2)):
                    d=required_distance(factor*primary['A_component_mean'],primary['Lambda'],volume,budget)
                    a_lo,a_hi=primary['A_component_mean_profile97_5']
                    l_lo,l_hi=primary['Lambda_profile97_5']
                    low=required_distance(factor*a_lo,l_hi,volume,budget)
                    high=required_distance(factor*a_hi,l_lo,volume,budget)
                    comparisons.append(dict(p=p,tiles_times_cycles=volume,budget=budget,convention=convention,
                        assumed_distance=baseline,calibrated_distance=d,distance_change=None if d is None else d-baseline,
                        qubit_ratio=None if d is None else (d/baseline)**2,
                        distance_interval95_bonferroni=[low,high],
                        qubit_ratio_interval95_bonferroni=[None if k is None else (k/baseline)**2 for k in (low,high)],
                        interval_note='Conservative rectangle from two 97.5% profile intervals; asymptotic and conditional on the fitted law.'))
    cell_summaries=[]
    for row in rows:
        c=row['contract']
        cell_summaries.append({**{k:c[k] for k in ('p','distance','basis','rounds','seed')},
            **{k:row[k] for k in ('shots','failures','finite_window_per_round','finite_window_interval95')},
            'fewer_than_30_failures':row['failures']<30})
    joint=LawFit(rows,endpoints=True).summary()
    joint['likelihood_ratio_vs_separate_p_fits']=2*(joint['negative_log_likelihood']-sum(v['negative_log_likelihood'] for v in results.values()))
    joint['likelihood_ratio_degrees_of_freedom']=3
    joint['likelihood_ratio_p_asymptotic']=float(chi2.sf(max(0.,joint['likelihood_ratio_vs_separate_p_fits']),3))
    holdout=holdout_analysis(results)
    for row in comparisons:
        p=str(row['p']);f=results[p]
        alternatives=[f,f['zero_endpoint_fit'],*f['length_sensitivity_fits'].values()]
        if 'distance_5_and_above_fit' in f:
            alternatives.append(f['distance_5_and_above_fit'])
        if holdout['status']=='complete':
            alternatives.append(holdout['fits'][p])
            if 'distance_5_and_above_fit' in holdout['fits'][p]:
                alternatives.append(holdout['fits'][p]['distance_5_and_above_fit'])
        factor=2 if row['convention']=='patch_union_bound' else 1
        ds=[required_distance(factor*v['A_component_mean'],v['Lambda'],row['tiles_times_cycles'],row['budget']) for v in alternatives]
        row['fit_choice_distance_range']=[min(ds),max(ds)]
        row['fit_choice_qubit_ratio_range']=[(min(ds)/row['assumed_distance'])**2,(max(ds)/row['assumed_distance'])**2]
        row['fit_choice_note']='Sensitivity over primary, zero endpoint, 10d, 20d, d>=5 where available, and independent 100d fit when complete; not a confidence interval or rigorous systematic bound.'
        if holdout['status']=='complete':
            long_fit=holdout['fits'][p]
            long_fit=long_fit.get('distance_5_and_above_fit',long_fit)
            long_d=required_distance(factor*long_fit['A_component_mean'],long_fit['Lambda'],row['tiles_times_cycles'],row['budget'])
            row['long_memory_sensitivity_distance']=long_d
            row['long_memory_sensitivity_distance_change']=long_d-row['assumed_distance']
            row['long_memory_sensitivity_qubit_ratio']=(long_d/row['assumed_distance'])**2
            row['long_memory_sensitivity_scope']='Post hoc 100d fit: d>=5 at 1e-3, d=3,5 at 1e-4; endpoint correction neglected.'
    result=dict(schema_version=1,status='complete' if holdout['status']=='complete' else 'training complete; independent validation pending',model='q_b=A_b Lambda^(-(d+1)/2); P_fail(R)=[1-(1-2q_b)^(R+delta_db)]/2',
        fit_scope='Common distance suppression for X/Z, separate basis amplitudes; independent endpoint correction per distance and basis.',
        patch_convention='A_patch_union_bound=A_x+A_z=2*A_component_mean; X/Z memory tests do not resolve Y, so this is a conservative sum, not a measured any-Pauli probability.',
        interval_convention='Cell rates: exact equal-tail 95% Clopper-Pearson. Fits: asymptotic profile-likelihood intervals; model/schedule/extrapolation uncertainty excluded.',
        extrapolation_warning='Distances outside 3..9, surgery, and arbitrarily long memory are not simulated. Effective p_threshold is a low-p fit parameter, not a measured threshold crossing.',
        schedule=dict(cz_layers=4,measurement_ns=100,reset_ns=100,interaction_slot_ns=10,nominal_round_ns=240,
            idle='DEPOLARIZE1(p) once per unused CZ slot and once per data-qubit measurement/reset slot',
            basis_changes='Explicit H gates, DEPOLARIZE1(p/10) each, fused into CZ slots; extra physical duration is unmodeled',
            init_readout='basis-appropriate independent flip p',two_qubit='DEPOLARIZE2(p)'),
        reference_reuse='No existing growth data pooled: old controls use p/10 idle and 80 ns allowance, differing from C6/C7.',
        sources={name:sha(ROOT/name) for name in ('PLAN.md','configs/assumptions.json','python/star_growth_probe.py','python/surface_law_calibration.py')},
        worker_sha256=manifest['worker_sha256'],manifest_sha256=sha(OUT/'manifest.json'),
        total_shots=sum(r['shots'] for r in rows),total_failures=sum(r['failures'] for r in rows),
        cells=cell_summaries,endpoint_slope_estimates=endpoint_slopes(rows),
        fits=results,joint_p_fit=joint,holdout=holdout,resource_comparisons=comparisons,
        shared_law_interpretation='Apply the same law to STAR and FT; it cancels as a differential law assumption, but differing protected volumes and distance rounding can still change their relative costs. It does not cancel against NISQ.',
        sampling_library_references=['https://github.com/quantumlib/Stim','https://pymatching.readthedocs.io/en/stable/'])
    if (OUT/'gof_bootstrap.json').exists():
        result['goodness_of_fit_bootstrap']=json.loads((OUT/'gof_bootstrap.json').read_text())
        assert result['goodness_of_fit_bootstrap']['manifest_sha256']==sha(OUT/'manifest.json')
    atomic(OUT/'fit.json',result)
    if holdout['status']=='complete':
        write_report(result)
    print(json.dumps(results,indent=2),flush=True)


def write_report(result):
    def ci(point, bounds, digits=3):
        return f'{point:.{digits}g} [{bounds[0]:.{digits}g}, {bounds[1]:.{digits}g}]'
    lines=['# Surface-law calibration — 23 September 2026', '',
        '**The reference suppression Λ=10 at p=10⁻³ is optimistic at larger distances.** '
        'The exponential is an approximate resource model: finite-distance lack of fit prevents treating its statistical intervals as prediction guarantees.', '',
        '**Experiment.** 78 million fixed training shots (R=10d,20d), plus 18 million independent R=100d holdout shots; '
        'd=3,5,7,9 at p=10⁻³ and d=3,5,7 at 10⁻⁴ (holdout d=3,5). '
        '[Stim](https://github.com/quantumlib/Stim) 1.16.dev0 and '
        '[PyMatching](https://pymatching.readthedocs.io/en/stable/) 2.3.1, ordinary MWPM, one serial process/node. '
        'Circuit controls and distance checks passed.', '',
        '**C6/C7 implementation.** CZ depolarizing p; every H depolarizing p/10; preparation/readout flips p; '
        'idle depolarizing p per unused CZ slot and per data-qubit measurement/reset slot. '
        'Four CZ slots + 100 ns measurement + 100 ns reset give nominal 240 ns. '
        '**Basis changes are fused into these slots; their additional physical duration is unmodeled.** '
        'Existing growth controls use different idle conventions and were not pooled.', '',
        '**Fit.** q_b(d)=A_b Λ^{−(d+1)/2}; P_b(R)=[1−(1−2q_b)^{R+δ_db}]/2. '
        'X/Z have separate amplitudes, common Λ, and fitted endpoint corrections δ. '
        'AΣ=A_X+A_Z is the conservative Pauli-component sum convention; basis-mean A=AΣ/2. '
        'It is not a directly measured any-Pauli rate.', '',
        '| p | Λ [95% profile interval] | AΣ [95% profile interval] |',
        '|---|---:|---:|']
    for p in ('0.001','0.0001'):
        f=result['fits'][p]
        lines.append(f"| {p} | {ci(f['Lambda'],f['Lambda_profile95'])} | {ci(f['A_patch_union_bound'],f['A_patch_union_bound_profile95'])} |")
    high=result['fits']['0.001'];low=result['fits']['0.0001']
    hold=result['holdout']['fits']
    lines += ['',
        f"Cell intervals are exact Clopper–Pearson; profile intervals are asymptotic and conditional on this law. "
        f"At 10⁻³ the training-law lack-of-fit p={high['goodness_of_fit_p_asymptotic']:.3g}; "
        f"100d also rejects a single d=3–9 law (p={hold['0.001']['goodness_of_fit_p_asymptotic']:.2g}). "
        f"A **post hoc d≥5,100d sensitivity fit** gives Λ={ci(hold['0.001']['distance_5_and_above_fit']['Lambda'],hold['0.001']['distance_5_and_above_fit']['Lambda_profile95'])}, "
        f"AΣ={ci(hold['0.001']['distance_5_and_above_fit']['A_patch_union_bound'],hold['0.001']['distance_5_and_above_fit']['A_patch_union_bound_profile95'])} "
        '(endpoint correction neglected; fit p=0.17). '
        'At 10⁻⁴,d=7 all four cells have zero failures; each million-shot R=140 cell gives q_b<2.64×10⁻⁸ (two-sided 95% interval upper endpoint). '
        'Sparse endpoint parameters hit bounds. pΛ is an effective fit threshold, not a measured threshold crossing.', '',
        '**Conditional resource extrapolation.** Smallest odd d≥3 with VAΣΛ^{−(d+1)/2}≤ε, V=tiles×cycles. '
        'Entries are assumed→fitted d (d² ratio); every entry extrapolates beyond sampled distances. '
        'JSON also gives statistical intervals, fit-choice sensitivity, and basis-mean comparisons.', '',
        '| p | V | ε=10⁻⁴ | ε=10⁻⁵ | ε=10⁻⁶ |',
        '|---|---:|---:|---:|---:|']
    for p in (.001,.0001):
        for volume in (1e8,1e11,1e14):
            selected=[r for r in result['resource_comparisons'] if r['p']==p and r['tiles_times_cycles']==volume and r['convention']=='patch_union_bound']
            cells=[f"{r['assumed_distance']}→{r['calibrated_distance']} ({r['qubit_ratio']:.2f})" for r in selected]
            lines.append(f'| {p:g} | {volume:.0e} | '+ ' | '.join(cells)+' |')
    lines += ['',
        '**Scope.** The same law applies to STAR and FT, so it creates no arm-specific advantage; '
        'different volumes and rounding can still change their relative costs. It does not cancel against NISQ. '
        'Memory only: surgery, correlated/coherent hardware noise, and high-distance performance are unvalidated.', '',
        '**Reproduce:** `~/soft/mambaforge/envs/stim/bin/python python/surface_law_calibration.py fit` then `audit`. '
        '[Evidence and manifests](../results/surface_law_v1/fit.json) retain seeds, circuits, counts, intervals and PBS scripts.']
    (ROOT/'docs/SURFACE_LAW.md').write_text('\n'.join(lines)+'\n')


def audit():
    manifest,rows=load_rows()
    before=json.loads((OUT/'input_snapshot.json').read_text())
    changed=[name for name,digest in before.items() if not (ROOT/name).exists() or sha(ROOT/name)!=digest]
    # The user explicitly has concurrent work in this checkout. Record external
    # changes; only the frozen inputs used by this calibration must be unchanged.
    critical=['PLAN.md','configs/assumptions.json','python/star_growth_probe.py']
    assert not set(critical).intersection(changed),changed
    assert sha(OUT/'worker_snapshot.py')==manifest['worker_sha256']
    assert not list((OUT/'cases').glob('*.lock'))
    assert json.loads((OUT/'validation.json').read_text())['status']=='pass'
    fits=json.loads((OUT/'fit.json').read_text())
    assert fits['holdout']['status']=='complete'
    if 'goodness_of_fit_bootstrap' in fits:
        assert all(r['optimizer_unsuccessful']==0 for r in fits['goodness_of_fit_bootstrap']['results'].values())
    assert fits['total_shots']==sum(r['shots'] for r in rows)
    for p,f in fits['fits'].items():
        assert f['optimizer_success']
        assert f['Lambda_profile95'][0] <= f['Lambda'] <= f['Lambda_profile95'][1]
        assert f['A_component_mean_profile95'][0] <= f['A_component_mean'] <= f['A_component_mean_profile95'][1]
    for row in fits['resource_comparisons']:
        factor=2 if row['convention']=='patch_union_bound' else 1
        f=fits['fits'][str(row['p'])]
        d=row['calibrated_distance']
        logratio=math.log(row['tiles_times_cycles']*factor*f['A_component_mean']/row['budget'])
        assert logratio-((d+1)/2)*math.log(f['Lambda'])<=1e-10
        assert d==3 or logratio-((d-1)/2)*math.log(f['Lambda'])>0
    holdout_count=holdout_shots=0
    if fits['holdout']['status']=='complete':
        _,holdout_rows=load_rows('holdout_manifest.json','holdout')
        holdout_count=len(holdout_rows);holdout_shots=sum(r['shots'] for r in holdout_rows)
        assert holdout_shots==fits['holdout']['total_shots']
    assert not list(OUT.rglob('*.lock'))
    result=dict(status='pass',fixed_shot_cells=len(rows),shots=sum(r['shots'] for r in rows),
        failures=sum(r['failures'] for r in rows),snapshot_files_unchanged=len(before)-len(changed),
        concurrent_changes_observed=changed,critical_calibration_inputs_unchanged=critical,
        task_write_scope=['python/surface_law_calibration.py','results/surface_law_v1/','docs/SURFACE_LAW.md'],
        runtime_seconds_sum=sum(r['elapsed_seconds'] for r in rows),
        circuit_hashes_verified=True,manifest_verified=True,complete_census=True,
        holdout_cells=holdout_count,holdout_shots=holdout_shots,
        active_sampling_locks=len(list(OUT.rglob('*.lock'))),fit_sha256=sha(OUT/'fit.json'))
    atomic(OUT/'audit.json',result)
    print(json.dumps(result),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['validate', 'sample', 'fit', 'audit','bootstrap'])
    parser.add_argument('--distance', type=int, default=3)
    parser.add_argument('--basis', choices=['x','z'], default='x')
    parser.add_argument('--p', type=float, default=.001)
    parser.add_argument('--rounds', type=int, default=30)
    parser.add_argument('--shots', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=20260923)
    parser.add_argument('--batch', type=int, default=4096)
    parser.add_argument('--group', default='pilot')
    parser.add_argument('--bootstraps', type=int, default=399)
    args = parser.parse_args()
    if args.command == 'validate':
        validate()
    elif args.command == 'sample':
        sample(args)
    elif args.command == 'fit':
        fit()
    elif args.command == 'bootstrap':
        bootstrap_gof(args.bootstraps)
    else:
        audit()


if __name__ == '__main__':
    main()
