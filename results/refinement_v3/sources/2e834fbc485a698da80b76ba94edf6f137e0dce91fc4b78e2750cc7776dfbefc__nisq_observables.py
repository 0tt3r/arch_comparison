"""Observable-specific PEC exposure and exactly propagated spin parities.

Postselection acceptance is not a certificate of small observable bias.
No state-vector Hubbard dynamics or Trotter calibration is used here.
"""
from math import exp, expm1, log, log1p, pi

from model import log_gamma2_dephasing, log_gamma2_depolarizing


def idle_probability(c):
    return -.75*expm1(c['gate_seconds']/1e-7*log1p(-4*c.get('idle_p_per_100ns',c['p'])/3))


def observable_exposure(events,schedule,observed,c):
    """Exact circuit-dependency cone, conservative about operator commutation."""
    active=set(observed)
    latest=[schedule['depth_layers']]*schedule['width']
    counts={'two_qubit':0,'one_qubit':0,'idle':0}
    for start,gate in sorted(events,key=lambda item:item[0],reverse=True):
        involved=active.intersection(gate.qubits)
        for q in involved:counts['idle']+=latest[q]-start-1
        if involved:
            counts['two_qubit' if len(gate.qubits)==2 else 'one_qubit']+=1
            active.update(gate.qubits)
        for q in gate.qubits:latest[q]=start
    counts['idle']+=sum(latest[q] for q in active)
    factors={'two_qubit':log_gamma2_depolarizing(c['p'],2),
             'one_qubit':log_gamma2_depolarizing(c['p']/10,1),
             'idle':log_gamma2_depolarizing(idle_probability(c),1)}
    pieces={key:count*factors[key] for key,count in counts.items()}
    pieces['spam']=(len(active)+len(observed))*log_gamma2_dephasing(c['p'])
    return {'chi':sum(pieces.values()),'exposure':pieces,'locations':counts,
            'input_support':sorted(active)}


def _backward_pauli(x,z,gate):
    if gate.name in ('x','z'):return
    if gate.name=='h':
        q=gate.qubits[0];x[q],z[q]=z[q],x[q];return
    if gate.name in ('s','sdg'):
        q=gate.qubits[0];z[q]^=x[q];return
    if gate.name=='rz':
        q=gate.qubits[0]
        if not x[q]:return
        quarters=round(gate.angle/(pi/2))
        if abs(gate.angle-quarters*pi/2)>1e-9:
            raise ValueError('Parity does not stay Pauli through this Rz')
        if quarters%2:z[q]^=x[q]
        return
    if gate.name=='rx':
        q=gate.qubits[0]
        if not z[q]:return
        quarters=round(gate.angle/(pi/2))
        if abs(gate.angle-quarters*pi/2)>1e-9:
            raise ValueError('Parity does not stay Pauli through this Rx')
        if quarters%2:x[q]^=z[q]
        return
    if gate.name=='rzz':
        a,b=gate.qubits
        if not (x[a]^x[b]):return
        quarters=round(gate.angle/(pi/2))
        if abs(gate.angle-quarters*pi/2)>1e-9:
            raise ValueError('Parity does not stay Pauli through this RZZ')
        if quarters%2:z[a]^=1;z[b]^=1
        return
    raise ValueError(gate.name)


def parity_acceptance(events,schedule,mode_to_wire,c):
    """Exact two-spin-parity acceptance for this independent Pauli noise model.

    Track all three parity characters through the actual primitive gates.
    If a non-Clifford gate spreads a character beyond one Pauli, reject the
    shortcut rather than assuming a Pauli channel at the block boundary.
    """
    n=schedule['width'];p=c['p'];idle=idle_probability(c)
    eigenvalues=[];log_visibilities=[]
    for mask in (1,2,3):
        x=[0]*n;z=[0]*n
        for mode,wire in mode_to_wire.items():
            if mask&(1<<(mode%2)):z[wire]=1
        log_visibility=sum(z)*log1p(-2*p)  # classical final readout flips
        latest=[schedule['depth_layers']]*n
        for start,gate in sorted(events,key=lambda item:item[0],reverse=True):
            for q in gate.qubits:
                if x[q] or z[q]:
                    log_visibility+=(latest[q]-start-1)*log1p(-4*idle/3)
            k=len(gate.qubits)
            if any(x[q] or z[q] for q in gate.qubits):
                pg=p if k==2 else p/10
                log_visibility+=log1p(-4**k*pg/(4**k-1))
            _backward_pauli(x,z,gate)
            for q in gate.qubits:latest[q]=start
        assert not any(x),'Initial parity is not diagonal in the occupation basis'
        log_visibility+=sum(z[q]*latest[q] for q in range(n))*log1p(-4*idle/3)
        log_visibility+=sum(z)*log1p(-2*p)  # initial preparation bit flips
        log_visibilities.append(log_visibility)
        eigenvalues.append(exp(log_visibility))
    acceptance=(1+sum(eigenvalues))/4
    no_fault_log=(schedule['two_qubit_gates']*log1p(-p)
                  +schedule['one_qubit_gates']*log1p(-p/10)
                  +schedule['idle_qubit_layers']*log1p(-idle)
                  +2*n*log1p(-p))
    assert no_fault_log<=log(acceptance)+1e-12
    return {'acceptance':acceptance,'character_visibilities':eigenvalues,
            'log_character_visibilities':log_visibilities,'log_no_fault_probability':no_fault_log,
            'conditional_bias_upper_bound':-expm1(min(0,no_fault_log-log(acceptance))),
            'bias_status':'worst-case bound; postselection alone is not precision-certified',
            'selection':'both spin parities; not exact spin-resolved particle numbers'}
