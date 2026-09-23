"""Read-only adapters for the two supplied workshop models.

Normalize task settings and steps, expose remaining compilation/delivery
differences, and evaluate all injection laws on an identical rotation ledger.
This is a diagnostic, not three independent hardware validations.
"""
import csv
import argparse
import hashlib
import json
from math import ceil, isclose, log, sqrt
from pathlib import Path
import sys
import zipfile
from unittest.mock import patch

from discussion_model import estimate, ideal_shots
from model import log_gamma2_dephasing, workload
from physics_sweep import configuration, planning_coefficient
from project_paths import ROOT

OUT=ROOT/'results/workshop_reconciliation_v1'


def common_contract(tau):
    epsilon=.01
    C=max(planning_coefficient(4,tau,a) for a in ('plaquette','termwise'))
    r=max(ceil(8*tau),ceil(sqrt(C/(.25*epsilon))))
    c={**configuration(.001),'stat_fraction':.4,'logical_fraction':.1,
       'synthesis_fraction':.075,'magic_fraction':.075}
    return dict(L=4,tau=tau,absolute_error=epsilon,steps=r,coefficient=C,Q=1000000,
                ideal_shots=ideal_shots(c,epsilon),hardware=c)


def delivery_row(contract,package,arm,q,seconds,chi,**details):
    Q=contract['Q'];copies=int(Q//q)
    logtime=(log(contract['ideal_shots'])+chi+log(seconds)-log(copies)) if copies else None
    return dict(package=package,arm=arm,tau=contract['tau'],r=contract['steps'],
                absolute_error=contract['absolute_error'],qubits_per_copy=q,seconds_per_attempt=seconds,
                chi=chi,ideal_shots=contract['ideal_shots'],copies_at_Q=copies,Q=Q,
                log10_runtime_at_Q=logtime/log(10) if logtime is not None else None,**details)


def our_rows(contract):
    out=[];c=contract['hardware']
    for arm,scenario,changes in [('nisq','effective',{}),('nisq','two-qubit PEC',{'nisq_model':'pec_two_qubit'}),
                                 ('star','count',{'star_law':'original_count'}),('star','angle',{}),('ft','individual',{}),
                                 ('ft','HWP 8',{'rotation_synthesis':'hamming_weight','hwp_batch_size':8}),('hybrid','angle',{})]:
        candidates=[]
        for factories in c['factory_choices'] if arm in ('ft','hybrid') else [1]:
            trial=estimate(4,contract['tau'],arm,{**c,**changes,'factory_choices':[factories]},.01,steps=contract['steps'])
            copies=contract['Q']//trial.get('qubits_per_copy',contract['Q']+1)
            if copies:candidates.append((trial['seconds_per_attempt']/copies,trial))
        v=min(candidates,key=lambda pair:pair[0])[1]
        row=delivery_row(contract,'ours',arm,v['qubits_per_copy'],v['seconds_per_attempt'],v['chi'],scenario=scenario)
        row.update({k:v.get(k) for k in ('two_qubit_gates','rotations','theta','T_count','distance','hwp_workspace')})
        if arm in ('ft','hybrid'):
            row['rotations']=workload(4,contract['tau'],contract['steps'],'plaquette').rotations
            row['factories']=v['factories']
        out.append(row)
    return out


def martin_rows(contract):
    from crossover.benchmarks.hubbard import HubbardBenchmark
    from crossover.architectures import NISQArchitecture,STARArchitecture,FTQCArchitecture
    b=HubbardBenchmark(eps=.01,t_evolve=contract['tau'])
    native_steps=b.trotter_steps(4)
    b.trotter_steps=lambda L:contract['steps']
    spec=b.circuit(4);out=[]
    physical=NISQArchitecture(p2=.001,p1=.0001,p_meas=.001,t_2q_s=1e-8,t_1q_s=1e-8,t_readout_s=1e-7)
    ex=physical.execute(spec,contract['Q'],.001)
    out.append(delivery_row(contract,'martin','nisq',ex.qubits_used,ex.t_shot_s,-2*log(ex.F),
        scenario='native attenuation rescaling',two_qubit_gates=ex.details['gates_2q'],native_steps=native_steps,
        residual_scope='Native T2 model retained; single-F response assumed'))
    star=STARArchitecture(p=.001,t_cycle_s=2.4e-7,logical_frac=1.)
    ex=star.execute(spec,contract['Q'],.001)
    theta=workload(4,contract['tau'],contract['steps'],'termwise').theta
    for name,chi in [('native count',2*ex.details['rotation_err']),
                     ('canonical count',spec.rotations*log_gamma2_dephasing(4*.001/15)),
                     ('angle',4*1.5*.001*theta+4*1e-5*spec.rotations)]:
        # Native logical faults were already budgeted by execute. Here chi is
        # the injection-only multiplier, using one shared statistical driver.
        out.append(delivery_row(contract,'martin','star',ex.qubits_used/ex.copies,ex.t_shot_s,chi,
            scenario=name,rotations=spec.rotations,theta=theta,distance=ex.details['d'],native_steps=native_steps,
            residual_scope='Native RUS/delivery and logical criterion retained; angle total supplied by common Hamiltonian'))
    ft=FTQCArchitecture(p=.001,t_cycle_s=2.4e-7,synth_frac=.3,logical_frac=.4,T_frac=.3)
    ex=ft.execute(spec,contract['Q'],.0025)
    out.append(delivery_row(contract,'martin','ft',ex.qubits_used/ex.copies,ex.t_shot_s,0,
        scenario='native individual synthesis',rotations=spec.rotations,T_count=ex.details['T_count'],
        distance=ex.details['d'],native_steps=native_steps,
        factory_wait_exposure_ratio=max(ex.details['cycles_logic'],ex.details['cycles_factory'])/ex.details['cycles_logic'],
        residual_scope='Native source menu; stall-time logical undercount reported, not silently corrected'))
    return out


def iqus_rows(contract):
    from fhcost.budget import DEFAULT
    from fhcost import hubbard,nisq,ftqc
    c=DEFAULT.but(p=.001,eps=.1,signal_regime='fixed',s_sig=.1,n_times=1,observable='local',
        n_ancilla=0,tmax_mode='const',tmax_const=contract['tau'],t_round=2.4e-7,dt_gate=1e-8,dt_meas=1e-7,
        frac_trotter=.25,frac_finite=0,frac_syn=.075,frac_logical=.1,frac_magic=.075,frac_mitig=.1,frac_stat=.4)
    native_steps=hubbard.trotter_steps(16,c)
    c=c.but(trotter_mode='fixed_density',steps_per_tau=contract['steps']/contract['tau'])
    out=[]
    # Full circuit: exclude the unbudgeted physical-cone discount. The supplied
    # lightcone_frac config field does not control this function in current code.
    with patch.object(hubbard,'mean_cone_fraction',return_value=1.):
        counts=hubbard.counts(16,c)
        assert isclose(hubbard.eps_absolute(c,16),contract['absolute_error'],rel_tol=1e-14)
        out.append(delivery_row(contract,'iqus','nisq',counts['q_per_copy'],counts['t_circuit'],nisq.log_gamma_sq(16,c),
            scenario='full-circuit two-qubit PEC',two_qubit_gates=counts['g_total'],native_steps=native_steps,
            residual_scope='Compact count/depth proxy; no idle or single-qubit noise'))
        theta=workload(4,contract['tau'],contract['steps'],'termwise').theta
        R=counts['n_rot']
        for name,chi in [('canonical count',R*log_gamma2_dephasing(4*.001/15)),
                         ('angle',4*1.5*.001*theta+4*1e-5*R)]:
            pt=ftqc.star_point(16,c)
            out.append(delivery_row(contract,'iqus','star',pt['phys'],pt['t_shot'],chi,scenario=name,
                rotations=R,theta=theta,distance=pt['d'],native_steps=native_steps,
                residual_scope='Native parallel-step delivery and logical criterion; same R in both injection laws'))
        pt=ftqc.surface_point(16,c)
        out.append(delivery_row(contract,'iqus','ft',pt['phys'],pt['t_shot'],0,scenario='native HWP',
            rotations=R,T_count=pt['n_t'],distance=pt['d'],hwp_workspace=pt['workspace'],native_steps=native_steps,
            residual_scope='Published-code workspace/factory/synthesis issues remain; this is a diagnostic output'))
    return out


def run(use_snapshot=False):
    sys.dont_write_bytecode=True
    OUT.mkdir(parents=True,exist_ok=True)
    source_root=ROOT
    if use_snapshot or not all((ROOT/name).exists() for name in ('crossover_martin','fh_iqus_idears')):
        source_root=ROOT/'build/workshop_sources';source_root.mkdir(parents=True,exist_ok=True)
        with zipfile.ZipFile(OUT/'source_snapshot.zip') as archive:
            for name in archive.namelist():
                path=Path(name)
                if path.is_absolute() or '..' in path.parts:raise ValueError('Invalid source-snapshot path')
            archive.extractall(source_root)
    sys.path[:0]=[str(source_root/'crossover_martin'),str(source_root/'fh_iqus_idears')]
    rows=[];ledgers=[];contracts=[]
    for tau in (1,2):
        contract=common_contract(tau);contracts.append(contract)
        rows.extend(our_rows(contract));rows.extend(martin_rows(contract));rows.extend(iqus_rows(contract))
        w=workload(4,tau,contract['steps'],'termwise')
        for law,chi in [('Martin native: chi = 0.4 p R',.4*.001*w.rotations),
                       ('IQUS native: chi = 1.066 p R',1.066*.001*w.rotations),
                       ('Canonical 2024 inverse',w.rotations*log_gamma2_dephasing(4*.001/15)),
                       ('Angle + floor',4*1.5*.001*w.theta+4*1e-5*w.rotations)]:
            ledgers.append(dict(tau=tau,r=w.steps,rotations=w.rotations,theta=w.theta,law=law,chi=chi))
    sources={}
    for folder,subdir in [('crossover_martin','crossover'),('fh_iqus_idears','fhcost')]:
        for path in sorted((source_root/folder/subdir).rglob('*.py')):
            sources[str(path.relative_to(source_root))]=hashlib.sha256(path.read_bytes()).hexdigest()
    with zipfile.ZipFile(OUT/'source_snapshot.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(sources):
            info=zipfile.ZipInfo(name,date_time=(1980,1,1,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            archive.writestr(info,(source_root/name).read_bytes())
    record=dict(contracts=contracts,rows=rows,identical_star_ledger=ledgers,source_sha256=sources,
                scope='Common task/step/shot convention; native circuit and architecture differences retained. Not an accuracy certificate or a fair ranking of fully specified machines.')
    (OUT/'comparison.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/'comparison.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    text=['# Workshop models at a shared operating point','',
          '4 × 4, tau = 1 and 2, p = 0.001, 0.01 absolute error, 10 ns elementary gates, 100 ns measurement, 240 ns syndrome cycle, one time point, one-week budget. All runtimes use the same Gaussian doublon planning count at D*=0.1 and statistical allowance 0.004. This is a common estimator assumption, not a validation of native mitigation bias.',
          'The common step count covers the larger of our plaquette and termwise planning coefficients using 25% of the error budget. It is transferred to the supplied cost models to isolate counts; their unspecified/different circuits and observables are not thereby accuracy-certified. Their native step estimates are retained in the CSV.',
          'The table preserves native gate counts, register sizes, layouts, delivery, logical criteria, and source menus. Only task settings, steps, elementary clocks, and the final shot driver are normalized. Our source count is selected for throughput at Q=1e6, matching Martin\'s objective; IQUS retains its native point construction. IQUS uses the full circuit rather than a physical-cone discount. Its fixed signal is explicitly set. Adapters do not edit either supplied package.','',
          '| tau | r | Package | Arm / scenario | G2 or R | T / shot | Q / copy | seconds / shot | chi | log10 seconds at Q=1e6 |',
          '|---:|---:|---|---|---:|---:|---:|---:|---:|---:|']
    for v in rows:
        number=v.get('two_qubit_gates') if v['arm']=='nisq' else v.get('rotations')
        nt=v.get('T_count');runtime=v['log10_runtime_at_Q']
        text.append(f'| {v["tau"]} | {v["r"]} | {v["package"]} | {v["arm"]}: {v["scenario"]} | {number if number is not None else "—"} | {round(nt) if nt is not None else "—"} | {v["qubits_per_copy"]:,.0f} | {v["seconds_per_attempt"]:.4g} | {v["chi"]:.3f} | {runtime:.3f} |')
    text+=['','## Identical STAR ledger, different laws','',
           '| tau | rotations | total angle | law | log variance multiplier |','|---:|---:|---:|---|---:|']
    for v in ledgers:text.append(f'| {v["tau"]} | {v["rotations"]} | {v["theta"]:.1f} | {v["law"]} | {v["chi"]:.4f} |')
    text+=['','The two native count-law coefficients are themselves different: 0.4 p R in Martin versus approximately 1.066 p R in IQUS and the canonical Akahoshi inverse. Therefore choosing “count law” alone does not reconcile the implementations.',
           'Native rotations per step and injection schedules also differ. The supplied models are useful alternative scenarios, not independent replications of one numerical result. The shared-ledger table is an algebra check, not an empirical collapse plot.',
           'Martin retains its single-attenuation NISQ estimator and T2 convention. Its FT logical criterion omits factory wait exposure where supply is slower. IQUS retains its unvalidated HWP workspace/synthesis accounting and factory table; its quoted FT runtime is diagnostic. No qLDPC result is imported into the nearest-neighbor hardware comparison.',
           'The IQUS native Trotter calibration belongs to the spin-correlator/dimer task and is not independent calibration of our doublon quench. Its reconciliation preset retains an old four-step workload and does not control all the switches it labels.',
           'Classical ED capacity is not a classical impossibility boundary. No frontier line or Lieb-Robinson certificate is imported without a compatible initial state, boundary-error theorem, and classical-method scope.',
           '','Source hashes are recorded in comparison.json; source_snapshot.zip preserves only the supplied Python modules used for reproduction. Run `python python/workshop_reconciliation.py`; when the supplied folders are absent, or with `--snapshot`, the adapter restores that snapshot under build/. The independent repository inside fh_iqus_idears is left untouched.','']
    (OUT/'summary.md').write_text('\n'.join(text))
    print('Wrote',len(rows),'common-point diagnostics and',len(ledgers),'shared-ledger comparisons')
    check()


def check():
    value=json.loads((OUT/'comparison.json').read_text())
    assert len(value['rows'])==32 and len(value['identical_star_ledger'])==8
    with zipfile.ZipFile(OUT/'source_snapshot.zip') as archive:
        assert set(archive.namelist())==set(value['source_sha256'])
        for name,sha in value['source_sha256'].items():
            assert hashlib.sha256(archive.read(name)).hexdigest()==sha
    for row in value['rows']:
        c=next(v for v in value['contracts'] if v['tau']==row['tau'])
        assert row['r']==c['steps'] and row['absolute_error']==.01
        assert row['ideal_shots']==c['ideal_shots']
        assert row['copies_at_Q']==int(row['Q']//row['qubits_per_copy'])
        actual=(log(c['ideal_shots'])+row['chi']+log(row['seconds_per_attempt'])-log(row['copies_at_Q']))/log(10)
        assert isclose(actual,row['log10_runtime_at_Q'],abs_tol=1e-12)
    print('PASS: shared tasks, shot driver, runtime arithmetic, and captured input sources')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true')
    parser.add_argument('--snapshot',action='store_true',help='Run from the captured source archive even when original folders exist')
    args=parser.parse_args()
    if args.check:check()
    else:run(use_snapshot=args.snapshot)
