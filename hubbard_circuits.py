"""Common OBC Hubbard circuits and an explicit nearest-neighbor JW realization.

No state-vector dynamics or Trotter-error calibration is performed here.
Pauli rotations use exp(-i theta P); native Rz gates use exp(-i phi Z/2).
"""
from collections import Counter
from dataclasses import dataclass
from math import ceil, isclose, pi, sqrt


@dataclass(frozen=True)
class Operation:
    name: str
    qubits: tuple[int, ...]
    angle: float = 0.0
    paulis: str = ""
    role: str = "clifford"


def snake_sites(L):
    return [(x, y) for y in range(L)
            for x in (range(L) if y % 2 == 0 else range(L - 1, -1, -1))]


def mode_map(L):
    return {(x, y, spin): 2 * i + spin
            for i, (x, y) in enumerate(snake_sites(L)) for spin in (0, 1)}


def parity_rotation(i, j, endpoints, theta, role="dynamic"):
    lo, hi = sorted((i, j))
    letters = ["Z"] * (hi - lo + 1)
    letters[i-lo], letters[j-lo] = endpoints
    return Operation("rotate", tuple(range(lo, hi+1)), theta, "".join(letters), role)


def fermionic_hadamard(i, j):
    """F a_i F†=(a_i+a_j)/sqrt(2), F a_j F†=(a_i-a_j)/sqrt(2).

    Two commuting pi/8 Pauli rotations and a Z, including JW parity.
    The ordering of i and j fixes which diagonal entry is negative.
    """
    return [parity_rotation(i, j, "XY", -pi/8, "fixed_T"),
            parity_rotation(i, j, "YX", pi/8, "fixed_T"),
            Operation("z", (j,))]


def plaquette_evolution(modes, duration):
    """Four modes in cyclic order, H = -sum_cycle(c†c+h.c.)."""
    a, b, c, d = modes
    basis = fermionic_hadamard(c, a) + fermionic_hadamard(b, d)
    center = [parity_rotation(b, c, "XX", -duration),
              parity_rotation(b, c, "YY", -duration)]
    # F pairs have disjoint fermionic support and commute, and each F=F†.
    return basis + center + basis


def tiles(L, offset):
    """Disjoint 2x2 tiles; shifted tiles are clipped to the OBC square."""
    result = []
    for y in range(-offset, L, 2):
        for x in range(-offset, L, 2):
            cycle = [(x,y), (x+1,y), (x+1,y+1), (x,y+1)]
            sites = [v for v in cycle if 0 <= v[0] < L and 0 <= v[1] < L]
            if len(sites) > 1:
                result.append(sites)
    return result


def _merge(sequence):
    result = []
    for op in sequence:
        if (result and op.name == result[-1].name == "rotate"
                and op.qubits == result[-1].qubits
                and op.paulis == result[-1].paulis and op.role == result[-1].role):
            previous = result.pop()
            result.append(Operation("rotate", op.qubits, op.angle+previous.angle, op.paulis, op.role))
        else:
            result.append(op)
    return result


def build_workload(L, tau, steps, algorithm="plaquette", u_over_j=4):
    if L < 2 or L % 2 or tau <= 0 or steps < 1:
        raise ValueError("Require even L >= 2, tau > 0, and steps >= 1")
    modes, N, dt = mode_map(L), L**2, tau / steps
    prep = [Operation("x", (modes[x,y,(x+y)%2],)) for x,y in snake_sites(L)]
    onsite = [Operation("rotate", (2*i, 2*i+1), u_over_j/4, "ZZ", "dynamic")
              for i in range(N)]
    if algorithm == "termwise":
        terms = []
        for x,y in snake_sites(L):
            for nx,ny in ((x+1,y), (x,y+1)):
                if nx < L and ny < L:
                    for spin in (0,1):
                        for endpoint in ("XX", "YY"):
                            terms.append(parity_rotation(modes[x,y,spin], modes[nx,ny,spin], endpoint, -.5))
        terms = onsite[:1] + terms + onsite[1:]
        half = [Operation(op.name, op.qubits, op.angle*dt/2, op.paulis, op.role) for op in terms]
        evolution = _merge((half + half[::-1]) * steps)
    elif algorithm == "plaquette":
        # Merge adjacent onsite blocks across steps before expanding their gates.
        stages = []
        for _ in range(steps):
            for name, duration in (("onsite",dt/2),("pink",dt/2),("gold",dt),("pink",dt/2),("onsite",dt/2)):
                if name == "gold" and L == 2:
                    continue
                if stages and stages[-1][0] == name:
                    stages[-1] = (name, stages[-1][1] + duration)
                else:
                    stages.append((name, duration))
        evolution = []
        for name, duration in stages:
            if name == "onsite":
                evolution.extend(Operation(op.name, op.qubits, op.angle*duration, op.paulis, op.role) for op in onsite)
            else:
                for tile in tiles(L, int(name == "gold")):
                    for spin in (0,1):
                        indices = [modes[x,y,spin] for x,y in tile]
                        if len(indices) == 4:
                            evolution.extend(plaquette_evolution(indices, duration))
                        elif len(indices) == 2:
                            evolution.extend(parity_rotation(*indices, endpoint, -duration/2)
                                             for endpoint in ("XX", "YY"))
                        else:
                            raise ValueError("Unexpected clipped plaquette")
    else:
        raise ValueError(algorithm)
    central = (modes[L//2,L//2,0], modes[L//2,L//2,1])
    return prep + evolution, central


def native_path(operations):
    """Parity ladders use only adjacent positions on the qubit path."""
    gates = []
    for op in operations:
        if op.name != "rotate":
            gates.append(op)
            continue
        assert list(op.qubits) == list(range(op.qubits[0], op.qubits[-1]+1))
        for q, p in zip(op.qubits, op.paulis):
            if p == "Y":
                gates.append(Operation("sdg", (q,)))
            if p in "XY":
                gates.append(Operation("h", (q,)))
        ladder = [Operation("cx", (a,b)) for a,b in zip(op.qubits, op.qubits[1:])]
        gates.extend(ladder)
        gates.append(Operation("rz", (op.qubits[-1],), 2*op.angle, role=op.role))
        gates.extend(ladder[::-1])
        for q, p in zip(op.qubits, op.paulis):
            if p in "XY":
                gates.append(Operation("h", (q,)))
            if p == "Y":
                gates.append(Operation("s", (q,)))
    return gates


def native_fermionic_route(operations, width, restore_order=False):
    """Greedy adjacent fSWAP routing of exactly the same fermionic blocks.

    Keep a live mode permutation and reuse it across rotations. Four-CNOT
    fSWAPs are explicit; hopping/fSWAP fusion is not assumed. This is a
    locality-conscious comparator, not an optimized fermionic swap network.
    """
    at=list(range(width))
    position=list(range(width))
    gates=[]
    swaps=0
    def fswap(left):
        nonlocal swaps
        right=left+1
        # SWAP followed by CZ, with CZ = H(target) CX H(target).
        gates.extend([Operation('cx',(left,right)),Operation('cx',(right,left)),
                      Operation('cx',(left,right)),Operation('h',(right,)),
                      Operation('cx',(left,right)),Operation('h',(right,))])
        at[left],at[right]=at[right],at[left]
        position[at[left]],position[at[right]]=left,right
        swaps+=1
    def adjacent(i,j):
        lo,hi=sorted((position[i],position[j]))
        for left in range(hi-1,lo,-1):
            fswap(left)
        return position[i],position[j]
    cursor=0
    while cursor < len(operations):
        op=operations[cursor]
        if op.name != 'rotate':
            gates.append(Operation(op.name,tuple(position[q] for q in op.qubits),op.angle,op.paulis,op.role))
            cursor+=1
        elif op.paulis=='ZZ' and len(op.qubits)==2:
            a,b=adjacent(*op.qubits)
            gates.extend(native_path([Operation('rotate',tuple(sorted((a,b))),op.angle,'ZZ',op.role)]))
            cursor+=1
        elif op.role=='fixed_T':
            second,phase=operations[cursor+1:cursor+3]
            assert second.role=='fixed_T' and phase.name=='z'
            j=phase.qubits[0]
            i=op.qubits[0] if op.qubits[-1]==j else op.qubits[-1]
            a,b=adjacent(i,j)
            gates.extend(native_path(fermionic_hadamard(a,b)))
            cursor+=3
        else:
            second=operations[cursor+1]
            assert op.paulis[0]==op.paulis[-1]
            assert second.paulis[0]==second.paulis[-1]
            assert {op.paulis[0],second.paulis[0]}=={'X','Y'}
            assert second.qubits==op.qubits and second.angle==op.angle
            a,b=adjacent(op.qubits[0],op.qubits[-1])
            gates.extend(native_path([parity_rotation(a,b,'XX',op.angle),
                                      parity_rotation(a,b,'YY',op.angle)]))
            cursor+=2
    if restore_order:
        for i in range(width):
            while position[i]>i:
                fswap(position[i]-1)
    return gates,position,swaps


def path_layout(width):
    side = ceil(sqrt(width))
    return snake_sites(side)[:width], side


def schedule_native(gates, width, gate_ns=10, reset_ns=100, measurement_ns=100):
    """ASAP schedule, preserving all per-qubit dependencies; final common readout.

    All modes are read out, so symmetry information is available. Counting all
    their readout errors for local D is unnecessary; the cost model uses two.
    No crosstalk, concurrent-edge exclusion, or controller bottleneck is assumed.
    """
    positions, side = path_layout(width)
    ready = [0] * width
    occupied = 0
    counts = Counter()
    for gate in gates:
        if len(gate.qubits) == 2:
            a,b = (positions[q] for q in gate.qubits)
            assert abs(a[0]-b[0]) + abs(a[1]-b[1]) == 1
        start = max(ready[q] for q in gate.qubits)
        for q in gate.qubits:
            ready[q] = start + 1
        occupied += len(gate.qubits)
        counts[gate.name] += 1
    depth = max(ready)
    return {"gates": dict(counts), "two_qubit_gates": counts["cx"],
            "one_qubit_gates": sum(counts.values())-counts["cx"],
            "depth_layers": depth, "idle_qubit_layers": width*depth-occupied,
            "gate_ns": gate_ns, "seconds_per_attempt": (reset_ns+depth*gate_ns+measurement_ns)*1e-9,
            "physical_copy_side": side, "physical_copy_qubits": side**2,
            "layout": positions, "schedule": "ASAP on a square-grid Hamiltonian path"}


def angle_summary(ops):
    histogram = Counter((op.role, round(abs(op.angle), 14)) for op in ops if op.name == "rotate")
    return {"dynamic_rotations": sum(n for (role,_), n in histogram.items() if role == "dynamic"),
            "fixed_T": sum(n for (role,_), n in histogram.items() if role == "fixed_T"),
            "theta_dynamic": sum(abs(op.angle) for op in ops if op.name == "rotate" and op.role == "dynamic"),
            "angle_histogram": [{"role": role, "theta": angle, "count": n}
                                for (role,angle),n in sorted(histogram.items())]}
