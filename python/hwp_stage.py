"""Explicit equal-angle onsite stage with temporary logical AND arithmetic.

The arithmetic uses b-popcount(b) clean work bits. Every AND is an exact
four-T isometry on a zero target; X measurement and a conditional CZ uncompute
it. Routing below is a constructive logical-patch circuit, not a decoded
physical lattice-surgery simulation.
"""
from dataclasses import asdict
from math import ceil, pi, sqrt

from hubbard_circuits import Operation


def weight_network(inputs,ancillas):
    columns={0:list(inputs)};compute=[];undo=[];used=0;weight_bits=[]
    column=0
    while columns.get(column):
        active=columns[column]
        while len(active)>1:
            count=3 if len(active)>=3 else 2
            operands=[active.pop(0) for _ in range(count)]
            a,b=operands[:2];t=ancillas[used];used+=1
            if count==2:
                forward=[Operation('and',(a,b,t)),Operation('cx',(a,b))]
                backward=[Operation('cx',(a,b)),Operation('unand',(a,b,t))]
                result=b
            else:
                c=operands[2]
                forward=[Operation('cx',(a,b)),Operation('cx',(a,c)),Operation('and',(b,c,t)),
                         Operation('cx',(a,t)),Operation('cx',(b,c)),Operation('cx',(a,c))]
                backward=[Operation('cx',(a,c)),Operation('cx',(b,c)),Operation('cx',(a,t)),
                          Operation('unand',(b,c,t)),Operation('cx',(a,c)),Operation('cx',(a,b))]
                result=c
            compute.extend(forward);undo=backward+undo
            active.append(result);columns.setdefault(column+1,[]).append(t)
        weight_bits.append((column,active[0]));column+=1
    assert used==len(inputs)-len(inputs).bit_count()
    return compute,weight_bits,undo


def and_gates(a,b,t):
    # Relative-phase Toffoli followed by S-dagger is exact when target is zero.
    return [Operation('h',(t,)),Operation('t',(t,)),Operation('cx',(b,t)),
            Operation('tdg',(t,)),Operation('cx',(a,t)),Operation('t',(t,)),
            Operation('cx',(b,t)),Operation('tdg',(t,)),Operation('h',(t,)),Operation('sdg',(t,))]


def stage_operations(sites,batch,theta):
    workspace=batch-batch.bit_count();ancillas=list(range(2*sites,2*sites+workspace));ops=[]
    # exp(-i theta Z_up Z_down) conjugates a Z rotation by a parity CNOT.
    for i in range(sites):ops.append(Operation('cx',(2*i,2*i+1)))
    for start in range(0,sites,batch):
        targets=[2*i+1 for i in range(start,min(start+batch,sites))]
        if len(targets)==1:
            ops.append(Operation('rz',(targets[0],),2*theta));continue
        forward,bits,reverse=weight_network(targets,ancillas)
        for op in forward:
            ops.extend(and_gates(*op.qubits) if op.name=='and' else [op])
        ops.extend(Operation('rz',(q,),2*theta*2**weight) for weight,q in bits)
        ops.extend(reverse)
    for i in range(sites):ops.append(Operation('cx',(2*i,2*i+1)))
    return ops,2*sites+workspace


def nearest_neighbor_route(ops,width):
    """Greedy token routing on a square patch grid, with explicit SWAPs.

    Each remote two-patch action routes one endpoint along a shortest path.
    Resulting dependencies are scheduled ASAP. Conditional CZs conservatively
    occupy their edge on both outcomes; every uncompute's feedforward is charged.
    """
    side=ceil(sqrt(width));at=list(range(side*side));position={q:q for q in at}
    gates=[]
    def emit(name,qubits,angle=0.,dependency=None):
        gates.append(dict(name=name,qubits=list(qubits),angle=angle,dependency=dependency))
    def pair(a,b,name,angle=0.,dependency=None):
        while abs(position[a]//side-position[b]//side)+abs(position[a]%side-position[b]%side)>1:
            here=position[a];target=position[b]
            nxt=here+((side if target//side>here//side else -side) if here//side!=target//side else (1 if target%side>here%side else -1))
            other=at[nxt]
            for x,y in ((here,nxt),(nxt,here),(here,nxt)):emit('cx',(x,y))
            at[here],at[nxt]=at[nxt],at[here];position[a]=nxt;position[other]=here
        emit(name,(position[a],position[b]),angle,dependency)
    for op in ops:
        if op.name=='unand':
            a,b,t=op.qubits
            marker=len(gates);emit('mx_reset',(position[t],))
            pair(a,b,'conditional_cz',dependency=marker)
        elif len(op.qubits)==2:pair(*op.qubits,op.name,op.angle)
        else:emit(op.name,(position[op.qubits[0]],),op.angle)
    return gates,dict(width=side*side,side=side,final_positions=position)


def schedule_patch_gates(gates,width,T_per_angle,ports,source_interval=150,distance=25,buffer_policy='bounded'):
    """Constructive logical-clock schedule with shared supply and dependencies.

    Adjacent CX/CZ takes three d-round clocks; H/S one; X-measure/reset two.
    Each synthesized T consumes one clock, a shared port, and a source state.
    Bounded mode holds one state per port and pauses production until consumed.
    Continuous mode reports the storage required by uninterrupted production.
    """
    if buffer_policy not in ('bounded','continuous'):raise ValueError(buffer_policy)
    ready=[0.]*width;port_ready=[0.]*ports;finish=[];events=[];consumed=[];source=source_interval/distance
    source_ready=[source]*ports
    for i,g in enumerate(gates):
        start=max(ready[q] for q in g['qubits'])
        if g['dependency'] is not None:start=max(start,finish[g['dependency']]+1)
        synthesis=T_per_angle[float(g['angle']).hex()] if g['name']=='rz' else dict(T=int(g['name'] in ('t','tdg')),clifford=0)
        nt=synthesis['T']
        if nt:
            time=start
            for _ in range(nt):
                lane=min(range(ports),key=lambda j:max(port_ready[j],source_ready[j]) if buffer_policy=='bounded' else port_ready[j])
                supply=source_ready[lane] if buffer_policy=='bounded' else (len(consumed)//ports+1)*source
                time=max(time,port_ready[lane],supply)
                consumed.append(time);port_ready[lane]=time+1;source_ready[lane]=time+source;time+=1
            end=time+synthesis['clifford']
        else:end=start+({'cx':3,'conditional_cz':3,'mx_reset':2}.get(g['name'],1))
        for q in g['qubits']:ready[q]=end
        finish.append(end);events.append(dict(index=i,start=start,end=end,**g))
    makespan=max(ready);peak=0
    for j,t in enumerate(sorted(consumed)):
        produced=ports*int((t+1e-12)//source)
        peak=max(peak,min(len(consumed),produced)-j)
    if buffer_policy=='bounded':peak=min(ports,len(consumed))
    return dict(logical_clocks=makespan,seconds=makespan*distance*240e-9,
                T_consumed=len(consumed),required_buffer_upper_estimate=peak,
                source_only_clocks=ceil(len(consumed)/ports)*source,buffer_policy=buffer_policy,events=events)
