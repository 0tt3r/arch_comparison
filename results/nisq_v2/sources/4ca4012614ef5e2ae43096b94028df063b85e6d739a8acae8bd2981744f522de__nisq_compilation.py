"""Native RZZ compilation and exact fermionic causal cones.

The logical two-mode blocks are extracted from the v1 workload without changing
their order. All fermionic swaps are adjacent in the live JW order. A separate
Pauli-tree implementation uses extra edges of the physical square grid.
"""
from collections import Counter, deque
from dataclasses import dataclass
from math import atan2, ceil, pi, sqrt

import numpy as np

from hubbard_circuits import Operation, path_layout


@dataclass(frozen=True)
class Block:
    kind: str
    modes: tuple
    theta: float = 0.


def fermionic_blocks(operations):
    prep=[]
    blocks=[]
    cursor=0
    while cursor<len(operations):
        op=operations[cursor]
        if op.name=='x':
            assert not blocks
            prep.extend(op.qubits)
            cursor+=1
        elif op.paulis=='ZZ' and len(op.qubits)==2:
            blocks.append(Block('density',op.qubits,op.angle))
            cursor+=1
        elif op.role=='fixed_T':
            other,phase=operations[cursor+1:cursor+3]
            assert other.role=='fixed_T' and phase.name=='z'
            j=phase.qubits[0]
            i=op.qubits[-1] if op.qubits[0]==j else op.qubits[0]
            blocks.append(Block('hadamard',(i,j)))
            cursor+=3
        else:
            other=operations[cursor+1]
            assert {op.paulis[0],other.paulis[0]}=={'X','Y'}
            assert op.paulis[0]==op.paulis[-1] and other.paulis[0]==other.paulis[-1]
            assert op.qubits==other.qubits and op.angle==other.angle
            blocks.append(Block('hop',(op.qubits[0],op.qubits[-1]),op.angle))
            cursor+=2
    return prep,blocks


def causal_cone(blocks,observed):
    """Conservative exact cone using disjoint-support commutation of even blocks.

    Compressed-mode simulation is appropriate for our initial occupation state
    and diagonal occupation observable. It is not a generic JW qubit partial trace.
    """
    active=set(observed)
    retained=[]
    for block in reversed(blocks):
        if active.intersection(block.modes):
            retained.append(block)
            active.update(block.modes)
    return retained[::-1],sorted(active)


def _odd_matrix(block,position):
    """One-particle matrix in physical basis (|01>, |10>); vacuum stays 1."""
    if block.kind=='hop':
        c,s=np.cos(2*block.theta),-1j*np.sin(2*block.theta)
        return np.array([[c,s],[s,c]])
    if block.kind=='hadamard':
        matrix=np.array([[-1,1],[1,1]],dtype=complex)/sqrt(2)
        return matrix if position[block.modes[0]]<position[block.modes[1]] else matrix[::-1,::-1]
    raise ValueError(block)


def gaussian_rzz(pair,matrix):
    """Any number-conserving two-mode Gaussian gate uses at most two RZZs.

    Write its odd block as exp(i delta) Rz(alpha) Rx(beta) Rz(gamma).
    Opposite single-qubit Z rotations implement odd-sector Rz; XX+YY supplies
    Rx. The equal final Z rotations supply the determinant on |11>.
    """
    a,b=pair
    delta=np.angle(np.linalg.det(matrix))/2
    v=matrix*np.exp(-1j*delta)
    beta=2*atan2(abs(v[0,1]),abs(v[0,0]))
    phase_a=np.angle(v[0,0]) if abs(v[0,0])>1e-12 else 0.
    phase_b=np.angle(v[0,1]) if abs(v[0,1])>1e-12 else -pi/2
    alpha=-phase_a-phase_b-pi/2
    gamma=-phase_a+phase_b+pi/2
    gates=[]
    def rz(q,angle):
        angle=(angle+pi)%(2*pi)-pi
        if abs(angle)>1e-12:
            gates.append(Operation('rz',(q,),angle))
    rz(a,gamma/2);rz(b,-gamma/2)
    if abs(beta)>1e-12:
        # RXX(beta/2), then RYY(beta/2), with RPP(phi)=exp(-i phi P/2).
        # H Sdg H is a single Rx(-pi/2) pulse; absorb the final S into Rz.
        gates.extend([Operation('h',(a,)),Operation('h',(b,)),
                      Operation('rzz',pair,beta/2),
                      Operation('rx',(a,),-pi/2),Operation('rx',(b,),-pi/2),
                      Operation('rzz',pair,beta/2),Operation('h',(a,)),Operation('h',(b,))])
    rz(a,delta+alpha/2+(pi/2 if abs(beta)>1e-12 else 0))
    rz(b,delta-alpha/2+(pi/2 if abs(beta)>1e-12 else 0))
    return gates


def simplify_single_qubit(gates,width):
    """Cancel/merge local gates through disjoint-wire operations, without refits."""
    pending=[[] for _ in range(width)]
    result=[]
    inverse={'h':'h','x':'x','z':'z','s':'sdg','sdg':'s'}
    def flush(q):
        result.extend(pending[q]);pending[q]=[]
    for gate in gates:
        if gate.name in ('s','sdg','z'):
            gate=Operation('rz',gate.qubits,{'s':pi/2,'sdg':-pi/2,'z':pi}[gate.name])
        if len(gate.qubits)==2:
            for q in gate.qubits:flush(q)
            result.append(gate)
        else:
            q=gate.qubits[0]
            stack=pending[q]
            if stack and inverse.get(gate.name)==stack[-1].name:
                stack.pop()
            elif stack and gate.name==stack[-1].name and gate.name in ('rz','rx'):
                angle=(stack.pop().angle+gate.angle+pi)%(2*pi)-pi
                if abs(angle)>1e-12:stack.append(Operation(gate.name,(q,),angle))
            else:
                stack.append(gate)
    for q in range(width):flush(q)
    return result


def fused_route(blocks,modes,occupied,lookahead=8):
    """Route the unchanged ordered blocks and fuse all available pair operations.

    Adjacent even fermionic operations on disjoint modes commute, permitting
    pair fusion through disjoint operations. Lookahead chooses which endpoint
    moves; it never reorders noncommuting target blocks.
    """
    at=list(modes)
    position={mode:i for i,mode in enumerate(at)}
    nodes=[]
    last=[None]*len(modes)
    swaps=0
    fused=0
    def append(pair,matrix=None,density=0.):
        nonlocal fused
        pair=tuple(sorted(pair))
        previous=last[pair[0]]
        if previous is not None and previous==last[pair[1]]:
            node=nodes[previous]
            if matrix is not None:node['matrix']=matrix@node['matrix']
            node['density']+=density
            fused+=1
        else:
            last[pair[0]]=last[pair[1]]=len(nodes)
            nodes.append({'pair':pair,'matrix':np.eye(2,dtype=complex) if matrix is None else matrix,
                          'density':density})
    def swap(left):
        nonlocal swaps
        append((left,left+1),np.array([[0,1],[1,0]],dtype=complex))
        at[left],at[left+1]=at[left+1],at[left]
        position[at[left]],position[at[left+1]]=left,left+1
        swaps+=1
    for index,block in enumerate(blocks):
        i,j=block.modes
        lo,hi=sorted((position[i],position[j]))
        if hi-lo>1:
            # Endpoint and balanced meeting positions, with deterministic ties.
            candidates=[]
            for meeting in sorted({lo,hi-1,(lo+hi-1)//2}):
                sequence=list(range(lo,meeting))+list(range(hi-1,meeting,-1))
                trial=at.copy()
                for left in sequence:trial[left],trial[left+1]=trial[left+1],trial[left]
                loc={mode:q for q,mode in enumerate(trial)}
                score=sum((lookahead-k)*max(0,abs(loc[next_block.modes[0]]-loc[next_block.modes[1]])-1)
                          for k,next_block in enumerate(blocks[index+1:index+1+lookahead]))
                candidates.append((score,meeting,sequence))
            for left in min(candidates)[2]:swap(left)
        pair=tuple(sorted((position[i],position[j])))
        if block.kind=='density':append(pair,density=block.theta)
        else:append(pair,_odd_matrix(block,position))
    gates=[Operation('x',(i,)) for i,mode in enumerate(modes) if mode in occupied]
    for node in nodes:
        gates.extend(gaussian_rzz(node['pair'],node['matrix']))
        if abs(node['density'])>1e-12:
            gates.append(Operation('rzz',node['pair'],2*node['density']))
    return simplify_single_qubit(gates,len(modes)),position,{'fermionic_swaps':swaps,'pair_fusions':fused}


def _cx_rzz(a,b):
    return [Operation('h',(b,)),Operation('rzz',(a,b),pi/2),
            Operation('rz',(a,),-pi/2),Operation('rz',(b,),-pi/2),Operation('h',(b,))]


def grid_parity_route(operations,width):
    """Two-root parity trees on the induced 2D support, then a native RZZ.

    This genuinely uses off-path grid edges. It permutes no fermionic modes;
    the complete JW parity support is retained.
    """
    positions,_=path_layout(width)
    gates=[]
    off_path=0
    for op in operations:
        if op.name!='rotate':
            gates.append(op);continue
        support=set(op.qubits)
        root=min(support,key=lambda q:sum(abs(positions[q][0]-positions[j][0])+abs(positions[q][1]-positions[j][1]) for j in support))
        parent={root:None};depth={root:0};queue=deque([root])
        while queue:
            q=queue.popleft()
            for j in sorted(support-parent.keys()):
                if sum(abs(a-b) for a,b in zip(positions[q],positions[j]))==1:
                    parent[j]=q;depth[j]=depth[q]+1;queue.append(j)
        assert len(parent)==len(support)
        # Every string has at least two sites. Pick a root-child edge whose
        # removal best balances the two trees (gate count is independent).
        children=[q for q in support if parent[q]==root]
        def subtree_size(child):
            count=0
            for q in support:
                while q is not None and q!=child:q=parent[q]
                count+=q==child
            return count
        child=min(children,key=lambda q:abs(2*subtree_size(q)-len(support)))
        for q,p in zip(op.qubits,op.paulis):
            if p=='Y':gates.append(Operation('sdg',(q,)))
            if p in 'XY':gates.append(Operation('h',(q,)))
        ladder=[(q,parent[q]) for q in sorted(support,key=lambda q:(-depth[q],q)) if q not in (root,child)]
        for a,b in ladder:gates.extend(_cx_rzz(a,b))
        gates.append(Operation('rzz',(root,child),2*op.angle))
        for a,b in reversed(ladder):gates.extend(_cx_rzz(a,b))
        for q,p in zip(op.qubits,op.paulis):
            if p in 'XY':gates.append(Operation('h',(q,)))
            if p=='Y':gates.append(Operation('s',(q,)))
        off_path+=sum(abs(a-b)>1 for a,b in ladder)*2+(abs(root-child)>1)
    return simplify_single_qubit(gates,width),{'off_path_entanglers_before_local_simplification':off_path}


def schedule(gates,width,gate_ns=10,reset_ns=100,measurement_ns=100):
    positions,side=path_layout(width)
    ready=[0]*width
    counts=Counter();events=[]
    for gate in gates:
        if len(gate.qubits)==2:
            a,b=(positions[q] for q in gate.qubits)
            assert abs(a[0]-b[0])+abs(a[1]-b[1])==1
        start=max(ready[q] for q in gate.qubits)
        for q in gate.qubits:ready[q]=start+1
        events.append((start,gate))
        counts[gate.name]+=1
    depth=max(ready,default=0)
    two=sum(len(g.qubits)==2 for g in gates)
    return {'width':width,'gates':dict(counts),'two_qubit_gates':two,
            'one_qubit_gates':len(gates)-two,'depth_layers':depth,
            'idle_qubit_layers':width*depth-sum(len(g.qubits) for g in gates),
            'seconds_per_attempt':(reset_ns+gate_ns*depth+measurement_ns)*1e-9,
            'physical_copy_qubits':side*side,'layout':positions},events


def cudaq_operations(gates):
    """Keep native RZZ explicit as a weight-two Pauli rotation in CUDA-Q P0."""
    return [Operation('rotate',g.qubits,g.angle/2,{'rzz':'ZZ','rx':'X'}[g.name],'native')
            if g.name in ('rzz','rx') else g for g in gates]
