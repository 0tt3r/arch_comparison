"""Cultivated-state supply for both FT and hybrid, with explicit handoff risk."""
import json
from math import ceil, log
from pathlib import Path

from compiled_model import compiled_star, encoded
from model import accuracy, finish, log_gamma2_dephasing, star_exposure, workload

ROOT=Path(__file__).resolve().parent


def cultivated_encoded(record,c,arm,source,handoff_error=0.,fixed_slot=False):
    if c['p']!=source['physical_error']:raise ValueError('Cultivation input point is only sourced at p=0.001')
    if arm not in ('ft','hybrid'):raise ValueError(arm)
    w=workload(record['L'],record['tau'],record['steps'],record['algorithm'],c['u_over_j'])
    eps,_=accuracy(c)
    nt=record['synthesis']['T_total'] if arm=='ft' else w.fixed_T
    analog=0 if arm=='ft' else w.rotations
    chi=0 if arm=='ft' else star_exposure(w,c,'star_floor')-nt*log_gamma2_dephasing(4*c['p']/15)
    log_amp=chi/2
    qin=source['state_error']+handoff_error
    volume=source['unshared_fixed_slot_volume_qubit_cycles'] if fixed_slot else source['packed_volume_qubit_cycles']
    slot=source['factory_physical_qubits'];candidates=[]
    for distill in (0,1):
        magic=qin if not distill else 35*qin**3
        if log(nt*magic)+log_amp>log(eps*c['magic_fraction']):continue
        for copies in c['factory_copies']:
            data_tiles=ceil(1.5*(2*w.sites+6))+(c['star_injection_tiles'] if analog else 0)
            # Source faults already enter magic; apply the surface-memory fit
            # only to data, handoff buffers, and any final distillation block.
            protected_tiles=data_tiles+copies*(1+12*distill)
            for d in range(source['reported_end_distance'],c['max_distance']+1,2):
                raw_service=max(volume/slot,source['successful_latency_cycles'])+source['handoff_clocks']*d
                service=(15*raw_service+11*d)/(1-15*qin) if distill else raw_service
                consume=(9*nt+c['star_clocks_per_rotation']*analog+(2*w.rotations if not analog else 0))*d
                supply=source['successful_latency_cycles']+ceil(nt/copies)*service
                seconds=(2*d+max(consume,supply))*c['cycle_seconds']+(nt+analog)*c['feedforward_seconds']
                seconds+=c['reset_seconds']+c['measurement_seconds']
                cycles=ceil(seconds/c['cycle_seconds'])
                log_fault=log(protected_tiles*cycles*c['surface_A'])+(d+1)/2*log(c['p']/c['surface_threshold'])
                if log_fault<=log(eps*c['logical_fraction'])-log_amp:
                    qcopy=2*protected_tiles*d*d+copies*slot
                    candidates.append(finish(w,arm,chi,qcopy,seconds,c,distance=d,T_count=nt,
                        factory_copies=copies,post_cultivation_distillation_levels=distill,
                        cultivated_input_error=source['state_error'],handoff_error=handoff_error,
                        magic_output_error=magic,source_volume_qubit_cycles=volume,
                        source_physical_qubits=slot,source_service_cycles=raw_service,
                        final_state_service_cycles=service,
                        consumption_cycles=consume,supply_cycles=supply,
                        implementation='compiled_T_counts_with_conditional_cultivation_supply'))
                    break
    if not candidates:raise ValueError('No cultivation configuration meets the declared budgets')
    return min(candidates,key=lambda row:row['log10_Q_min'])


def run():
    c=json.loads((ROOT/'assumptions.json').read_text())
    source=json.loads((ROOT/'cultivation_assumptions.json').read_text());rows=[]
    for L in (2,4,6,8):
        record=json.loads((ROOT/f'results/cudaq_v1/plaquette_L{L}_tau1_r4/result.json').read_text())
        for eps in (.001,.01):
            scenario={**c,'relative_error':eps/c['doublon_reference']}
            rows.append({**compiled_star(record,scenario),'source':'STAR reference','absolute_error':eps})
            for arm in ('ft','hybrid'):
                rows.append({**encoded(record,scenario,arm,injection_error=46*c['p']/15),
                             'source':'conventional distillation','absolute_error':eps})
                for error,fixed in [(0.,False),(1e-9,False),(1e-8,False),(0.,True)]:
                    row=cultivated_encoded(record,scenario,arm,source,error,fixed)
                    rows.append({**row,'source':'cultivation fixed slot' if fixed else 'cultivation packed',
                                 'absolute_error':eps})
    out=ROOT/'results/refinement_v3';out.mkdir(exist_ok=True,parents=True)
    (out/'cultivation.json').write_text(json.dumps(rows,indent=2)+'\n')
    for row in rows:
        if row['L']==4 and row.get('handoff_error',0)==0:
            print(row['absolute_error'],row['arm'],row['source'],round(10**row['log10_Q_min']),
                  'distill',row.get('post_cultivation_distillation_levels',row.get('factory_level')),flush=True)


if __name__=='__main__':run()
