"""Explicit alternatives to the frozen discussion model, not new calibrations.

Angles use exp(-i theta P). Hamming-weight phasing follows Campbell App. E.2;
its ancillas are reused between batches and charged throughout each shot.
"""
from math import ceil, log2

from model import log_gamma2_dephasing


def analog_exposure(w, c):
    law = c.get('star_law', 'angle')
    if law == 'angle':
        return 4*c['star_alpha']*c['p']*w.theta, 4*c['star_floor']*w.rotations
    if law == 'original_count':
        return w.rotations*log_gamma2_dephasing(4*c['p']/15), 0.
    raise ValueError(f'Unknown STAR law: {law}')


def source_error(c):
    """Allow lower-p scenarios only when their extrapolation is named."""
    p0 = c['nisq_anchor_p']
    policy = c.get('source_p_policy', 'anchor_only')
    if policy == 'anchor_only':
        if c['p'] != p0:
            raise ValueError('Cultivation reference is anchored only at p=0.001')
        return c['source_error']
    if not 0 < c['p'] <= p0:
        raise ValueError('Source scenarios support only 0 < p <= the anchor')
    if policy == 'held_at_anchor':
        return c['source_error']
    if policy == 'power3':
        return c['source_error']*(c['p']/p0)**3
    raise ValueError(f'Unknown source extrapolation: {policy}')


def plaquette_rotation_stages(w):
    """(equal-angle stage width, repetitions), without pooling across time."""
    if w.algorithm != 'plaquette':
        raise ValueError('HWP implementation covers the plaquette algorithm only')
    L, r = w.L, w.steps
    P, G, B = (L//2)**2, (L//2-1)**2, 2*(L-2)
    groups = [(L*L, r+1), (4*P, r if L == 2 else 2*r), (4*G, r), (4*B, r)]
    groups = [(width, repeat) for width, repeat in groups if width]
    assert sum(width*repeat for width, repeat in groups) == w.rotations
    return groups


def hamming_weight_counts(w, batch_size):
    if not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError('Positive integer HWP batch size required')
    rotations = toffolis = ancillas = 0
    for width, repeat in plaquette_rotation_stages(w):
        full, remainder = divmod(width, batch_size)
        for size, count in ((batch_size, full), (remainder, int(remainder > 0))):
            if not count:
                continue
            alpha = size-size.bit_count()
            rotations += repeat*count*size.bit_length()
            toffolis += repeat*count*alpha
            ancillas = max(ancillas, alpha)
    return dict(synthesized_rotations=rotations, hwp_toffolis=toffolis,
                hwp_workspace=ancillas, hwp_batch_size=batch_size)


def synthesis_plan(w, epsilon, c):
    style = c.get('rotation_synthesis', 'individual')
    details = {}
    if style == 'individual':
        rotations, toffolis = w.rotations, 0
    elif style == 'hamming_weight':
        details = hamming_weight_counts(w, c['hwp_batch_size'])
        rotations, toffolis = details['synthesized_rotations'], details['hwp_toffolis']
    else:
        raise ValueError(style)
    eta = epsilon*c['synthesis_fraction']/(2*rotations)
    t_per_rotation = ceil(3.02*log2(1/eta)+.77)
    return dict(T_count=w.fixed_T+4*toffolis+rotations*t_per_rotation,
                T_per_rotation=t_per_rotation,
                clifford_clock_units=2*w.rotations+2*toffolis, **details)


def compact_physical_counts(L, r, c):
    """Local-encoding count scenario; not a native-grid compilation.

Two spin copies use at most 3N encoding qubits. Weight-three hopping Paulis
cost four CNOTs with a parity ladder; weight-two onsite Paulis cost two.
The two-CNOT alternative assumes fusion. Preparation, packing, routing, and
utilization are declared sensitivities, not properties proved by the mapping.
"""
    N, E = L*L, 2*L*(L-1)
    hopping = 8*E*r
    onsite = 2*r*(N-1)+1
    prep = c.get('local_prep_gates_per_site', 8)*N
    routing = c.get('local_routing_multiplier', 1.)
    gates = ceil(routing*(c.get('local_gates_per_hopping_pauli', 4)*hopping+2*onsite)+prep)
    register = 3*N+2
    packed = ceil((c.get('local_layout_factor', 1.)*register)**.5)**2
    depth = max(4*L, ceil(2*gates/(c.get('local_gate_utilization', .5)*register)))
    return dict(two_qubit_gates=gates, depth=depth, qubits_per_copy=packed,
                encoding_qubits=3*N, preparation_two_qubit_gates=prep)
