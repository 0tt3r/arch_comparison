"""Same-Hamiltonian open-boundary LCU inventories and dynamics query budgets."""
from collections import defaultdict
import hashlib
import json
from math import asin, ceil, exp, lgamma, log, log2, sin, sqrt

import numpy as np
from scipy.special import jv

from hubbard_circuits import build_workload
from project_paths import ROOT

OUT=ROOT/'results/frontier_v1/qubitization'


class OpenHubbard:
    """Particle-hole centered Hubbard; differs by a scalar at half filling."""
    def __init__(self,L,u=4):
        self.L=L;self.u=u;self.terms=defaultdict(float)
        operations,_=build_workload(L,1.,1,'termwise',u)
        for op in operations:
            if op.name!='rotate':continue
            word=['I']*(2*L*L)
            for q,p in zip(op.qubits,op.paulis):word[q]=p
            self.terms[''.join(word)]+=op.angle
    def n_qubits(self):return 2*self.L*self.L
    def n_terms(self,**kwargs):return len(self.terms)
    @property
    def alpha(self):return sum(abs(v) for v in self.terms.values())
    def yield_PauliLCU_Info(self,do_pad=0,return_as='strings'):
        yield from self.terms.items()
        for _ in range(do_pad):yield 'I'*self.n_qubits(),0.


def bessel_degree(x,operator_error):
    """Degree with an absolute Jacobi--Anger tail and analytic remainder.

    |J_n(x)| <= (|x|/2)^n/n! follows from its Poisson integral.
    Beyond K the ratio of these majorants is <= x/(2(K+2)).
    The returned tail is for exp(-i x H/alpha), before LCU normalization.
    """
    x=abs(float(x))
    if x==0:return 0,0.
    K=max(16,ceil(2*x+40))
    def remainder(k):
        a=exp((k+1)*log(x/2)-lgamma(k+2))
        return 2*a/(1-x/(2*(k+2)))
    while remainder(K)>operator_error/100:K+=max(16,ceil(x/4))
    values=2*abs(jv(np.arange(1,K+1),x));tail=np.cumsum(values[::-1])[::-1]
    rest=remainder(K)
    for degree in range(K):
        bound=float(tail[degree]+rest)
        if bound<=operator_error:return degree,bound
    return K,rest


def encoding_inventory(L,split=False):
    # Imported only by the compiler entry point, not by scalar model/tests.
    from pyLIQTR.BlockEncodings.PauliStringLCU import PauliStringLCU
    from pyLIQTR.utils.resource_analysis import estimate_resources
    from pyLIQTR.qubitization.qubitized_gates import QubitizedRotation
    import cirq
    from qualtran.bloqs.multiplexers.select_pauli_lcu import SelectPauliLCU
    model=OpenHubbard(L)
    if split:
        class SplitHubbard(OpenHubbard):
            def __init__(self,L):
                super().__init__(L)
                self.expanded=[(word,np.sign(v)*.5) for word,v in self.terms.items() for _ in range(round(abs(v)*2))]
            def n_terms(self,**kwargs):return len(self.expanded)
            def yield_PauliLCU_Info(self,do_pad=0,return_as='strings'):
                yield from self.expanded
                for _ in range(do_pad):yield 'I'*self.n_qubits(),0.
        model=SplitHubbard(L)
    enc=PauliStringLCU(model,control_val=1)
    if split:
        from qualtran import GateWithRegisters, Signature
        from qualtran.bloqs.state_preparation import PrepareUniformSuperposition
        from pyLIQTR.utils.resource_analysis import pylqt_t_complexity
        class Uniform(GateWithRegisters):
            @property
            def signature(self):return Signature(enc.selection_registers)
            def _t_complexity_(self):return pylqt_t_complexity(PrepareUniformSuperposition(model.n_terms()))
            def decompose_from_registers(self,context,**quregs):
                yield PrepareUniformSuperposition(model.n_terms()).on_registers(target=quregs['selection'])
        enc._prepare_gate=Uniform()
    # pyLIQTR's generic adapter passes magnitudes (and zero padding) into
    # SELECT, which must be unitary. Magnitudes belong only in PREPARE.
    words=list(model.yield_PauliLCU_Info(do_pad=enc.n_pad))
    enc._select_gate=SelectPauliLCU(selection_bitsize=int(log2(enc.n_tot)),target_bitsize=model.n_qubits(),
        select_unitaries=tuple(cirq.DensePauliString(word,coefficient=-1 if coefficient<0 else 1) for word,coefficient in words),control_val=1)
    profile=estimate_resources(enc,profile=True)
    rotation_controls=sum(reg.bitsize for reg in (*enc.selection_registers,*enc.junk_registers))
    reflect=estimate_resources(QubitizedRotation(rotation_controls,control_val=1),profile=True)
    # Reserve k additional clean workspace patches beyond the public signature.
    # Unary iteration needs k-1; comparator/uniform preparation also fits k.
    # The two QSVT branch/phase bits are charged separately.
    return dict(L=L,variant='split_uniform' if split else 'generic_amplitudes',alpha=model.alpha,
                pauli_terms=len(model.terms),lcu_terms=model.n_terms(),encoding=profile,reflection=reflect,
                logical_qubits=model.n_qubits()+2*rotation_controls+3,
                selection_and_junk=rotation_controls,signature_qubits=sum(reg.total_bits() for reg in enc.signature),
                terms_sha256=hashlib.sha256(json.dumps(dict(model.terms),sort_keys=True).encode()).hexdigest())


# P_{2j+1}(a), for the accepted block after j oblivious amplification steps.
AMPLIFICATION=((0,1),(0,3,0,-4),(0,5,0,-20,0,16),(0,7,0,-56,0,112,0,-64))


def amplification(j,base_error):
    coefficients=AMPLIFICATION[j];rho=.25+base_error
    C=sum(n*abs(a)*rho**(n-1) for n,a in enumerate(coefficients) if n)
    a=sin((2*j+1)*asin(.25))
    error=C*base_error
    lower=a-error
    return dict(amplitude=a,lipschitz=C,block_error=error,acceptance_lower=lower**2,
                observable_error=4*error/lower)


def query_budget(inventory,tau,epsilon,simulation_fraction=.25,j=2):
    # Each polynomial is half of cosine/sine; the branch LCU adds another 1/2.
    # Allocate half the simulation budget to the Bessel tail and half to
    # classical phase realization. This is a constructive count, not generated
    # and verified phase angles for every large-degree circuit.
    a=sin((2*j+1)*asin(.25));C=sum(n*abs(v)*.26**(n-1) for n,v in enumerate(AMPLIFICATION[j]) if n)
    delta=min(.001,epsilon*simulation_fraction*a/(16*C))
    degree,tail=bessel_degree(inventory['alpha']*tau,2*delta)
    even=degree-degree%2;odd=degree if degree%2 else degree-1
    queries=even+odd;phases=queries+2
    amp=amplification(j,delta);repeats=2*j+1
    enc,ref=inventory['encoding'],inventory['reflection']
    # Extra acceptance reflections: conservative linear-in-workspace Clifford+T
    # decomposition of a multi-controlled Z, with clean ancillas already reserved.
    k=inventory['selection_and_junk']+2
    direct=repeats*(queries*enc['T']+phases*ref['T'])+2*j*8*k
    rotations=repeats*(queries*enc['Rotations']+phases*ref['Rotations'])
    clifford=repeats*(queries*enc['Clifford']+phases*ref['Clifford']+8)+2*j*40*k
    return dict(degree=degree,cos_degree=even,sin_degree=odd,block_queries=queries*repeats,
                raw_block_queries=queries,amplification_steps=j,polynomial_tail=tail,
                phase_realization_block_allowance=delta/2,base_block_error=delta,
                direct_T=int(direct),rotations=int(rotations),clifford=int(clifford),**amp)


def cost_candidate(inventory,budget,epsilon,c,synthesis_fraction,logical_fraction,cleanup):
    from discussion_model import sampling_cost
    from model import minimum_distance
    from resource_extensions import source_error
    amp=sqrt(budget['acceptance_lower'])
    synthesis_allowance=epsilon*synthesis_fraction
    synthesis_tolerance=synthesis_allowance*amp/(4+synthesis_allowance)
    acceptance=(amp-synthesis_tolerance)**2
    per_rotation=synthesis_tolerance/max(1,budget['rotations'])
    per_T=max(0,ceil(3.02*log2(1/per_rotation)+.77))
    nt=budget['direct_T']+per_T*budget['rotations']
    q0=source_error(c);qout=35*q0**3 if cleanup else q0
    # Conditional bias can be amplified by selection: charge 2q/P, including
    # faults outside the signal subspace, rather than silently treating them as
    # commuting Pauli errors of a deterministic circuit.
    magic=2*nt*qout/(acceptance*epsilon)
    simulation=budget['observable_error']/epsilon
    stat=1-simulation-synthesis_fraction-logical_fraction-magic
    if stat<=0:return None
    width=inventory['logical_qubits'];data=ceil(2*width+sqrt(8*width)+1)+6*c['parallel_lanes']
    choices=[]
    for factories in c['factory_choices']:
        tiles=data+factories*(1+12*cleanup)
        def duration(d):
            source=c['source_cycles_per_T']+2*d
            service=(15*source+11*d)/(1-15*q0) if cleanup else source
            supply=ceil(nt/factories)*service
            # Matched fast-block T consumption. Clifford depth is bracketed
            # by omitting extra serialization or charging the Clifford count
            # plus two Clifford slots per T. Neither is a routed SELECT schedule.
            consume=(nt+(budget['clifford']+2*nt)*c.get('qsp_clifford_serial_fraction',0))*d
            return (2*d+max(supply,consume))*c['cycle_seconds']+nt*c['feedforward_seconds']+c['reset_seconds']+c['measurement_seconds']
        solved=minimum_distance(tiles,duration,log(epsilon*logical_fraction*acceptance/2),c)
        if solved is None:continue
        d=max(solved[0],c['source_end_distance']);seconds=duration(d)
        result=sampling_cost(2*tiles*d*d+factories*c['source_qubits'],seconds,-log(acceptance*(1-epsilon)),{**c,'stat_fraction':stat},epsilon)
        if result['log10_Q_min'] is None:continue
        choices.append(dict(**result,distance=d,factories=factories,cleanup_levels=cleanup,T_count=nt,T_per_rotation=per_T,
            logical_qubits=width,data_tiles=data,attempt_acceptance_lower=acceptance*(1-epsilon),budget_fractions=dict(stat=stat,simulation=simulation,synthesis=synthesis_fraction,logical=logical_fraction,magic=magic)))
    return min(choices,key=lambda v:(v['log10_Q_min'],v['qubits_per_copy']*v['seconds_per_attempt'])) if choices else None


def run():
    from discussion_model import assumptions
    from budget_optimization import coefficient,optimize,settings
    config=json.loads((ROOT/'configs/frontier_study.json').read_text())['qubitization'];OUT.mkdir(parents=True,exist_ok=True)
    if config['u_over_j']!=4:raise ValueError('Equal-weight adapter is specialized to U/J=4')
    inventories=[]
    for L in config['sizes']:
        for split in (False,True):
            row=encoding_inventory(L,split);inventories.append(row);print('inventory',L,split,row['encoding'],flush=True)
    (OUT/'inventories.json').write_text(json.dumps(inventories,indent=2)+'\n')
    rows=[];c={**assumptions(),'data_block':'fast','ft_clocks_per_T':1,'T_consumption_ports':1}
    for L in config['sizes']:
      for tau in config['times']:
       for epsilon in config['errors']:
        for serial in (0.,1.):
         candidates=[]
         for inv in [v for v in inventories if v['L']==L]:
          for simulation in (.05,.15,.3):
           for j in (0,1,2,3):
            budget=query_budget(inv,tau,epsilon,simulation,j)
            for synthesis in settings()['synthesis_fraction_grid']:
             for logical in (.02,.1):
              for cleanup in (0,1):
               result=cost_candidate(inv,budget,epsilon,{**c,'qsp_clifford_serial_fraction':serial},synthesis,logical,cleanup)
               if result:candidates.append(dict(**result,**budget,variant=inv['variant']))
         best=min(candidates,key=lambda v:(v['log10_Q_min'],v['qubits_per_copy']*v['seconds_per_attempt']))
         rows.append(dict(L=L,tau=tau,epsilon=epsilon,clifford_serial_fraction=serial,**best))
        for hwp in (False,True):
            from physics_sweep import planning_coefficient
            local={**c,**({'rotation_synthesis':'hamming_weight','hwp_batch_size':L*L} if hwp else {})}
            # Existing model's calibrated envelope and optimized error allocation.
            C=coefficient(tau,'plaquette') if tau==8 and L>2 else planning_coefficient(L,tau,'plaquette')
            row=optimize(L,tau,'ft',local,epsilon,C_override=C)
            row['coefficient_scope']='3x3 running envelope only' if tau==8 and L>2 else 'extended finite-patch running envelope'
            rows.append(dict(**row,epsilon=epsilon,variant='trotter_hwp' if hwp else 'trotter_individual',clifford_serial_fraction=None))
        print('cost',L,tau,epsilon,flush=True)
    names=['python/qubitization_study.py','configs/frontier_study.json','python/hubbard_circuits.py','python/discussion_model.py','python/budget_optimization.py','python/resource_extensions.py','python/physics_sweep.py','python/model.py','python/dynamics_followup_analysis.py',
        'configs/assumptions.json','configs/discussion_model.json','configs/budget_optimization.json','configs/patch_sensitivity.json',
        'results/discussion_v3/calibration.json','results/discussion_v3/finite_step_guard.json',
        'results/discussion_v3/rectangle_checks/calibration.json','results/dynamics_followup_v1/measurements.json']
    import importlib.metadata
    record=dict(inventories=inventories,rows=rows,configuration=c,search_grid=dict(simulation_nominal=[.05,.15,.3],
        synthesis=settings()['synthesis_fraction_grid'],logical=[.02,.1],amplification_steps=[0,1,2,3],cleanup=[0,1],
        note='Cost scan overrides the two nominal single-query fractions recorded in frontier_study.json; charged simulation bias is the actual normalization bound.'),packages={n:importlib.metadata.version(n) for n in ('pyLIQTR','qualtran','cirq-core')},
        scope='Same OBC centered Hamiltonian. Constructive polynomial/query counts with an explicit phase-realization allowance; large-degree phases and routed SELECT are not compiled. Fast-block consumption sensitivity: no extra Clifford serialization or serial Clifford count plus two slots per T; neither is a routed SELECT schedule. Shared source/error model, conditional fault amplification included.',
        source_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    (OUT/'results.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':run()
