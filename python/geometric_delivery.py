"""Constructive nearest-neighbor magic transport for the compiled onsite stage.

A state walks from a perimeter output port via SWAPs, is consumed next to its
recipient, and the SWAPs are undone to restore the data mapping. All traversed
patches are reserved until the return. This deliberately conservative circuit
upper scenario replaces free transport; it is not an optimal surgery layout.
"""
from collections import deque
import gzip
import hashlib
import json
from math import ceil,sqrt

from hwp_stage import stage_operations,nearest_neighbor_route
from project_paths import ROOT

OUT=ROOT/'results/frontier_v1/delivery'


def path_to_neighbor(start,target,side):
    todo=deque([(start,[start])]);seen={start}
    while todo:
        node,path=todo.popleft()
        if abs(node[0]-target[0])+abs(node[1]-target[1])==1:return path
        for step in ((1,0),(-1,0),(0,1),(0,-1)):
            other=(node[0]+step[0],node[1]+step[1])
            if other!=target and other not in seen and all(0<=x<side for x in other):
                seen.add(other);todo.append((other,path+[other]))
    raise ValueError('No delivery path')


def footprint(path):
    nodes={(2*x,2*y) for x,y in path}
    nodes.update((a[0]+b[0],a[1]+b[1]) for a,b in zip(path,path[1:]))
    return nodes


def schedule(gates,data_side,angle_cost,ports,distance=25,source_rounds=150):
    side=data_side+2;fine=2*side-1
    centers=[(0,side//2),(side-1,side//2),(side//2,0),(side//2,side-1)][:ports]
    ready={};source_ready=[source_rounds/distance]*ports;events=[];finish=[];T=0;swaps=0
    def location(q):return (q//data_side+1,q%data_side+1)
    def time_for(nodes):return max((ready.get(v,0.) for v in nodes),default=0.)
    def reserve(nodes,t):
        for v in nodes:ready[v]=t
    for index,g in enumerate(gates):
        targets=[location(q) for q in g['qubits']];nodes=footprint(targets)
        start=time_for(nodes)
        if g['dependency'] is not None:start=max(start,finish[g['dependency']]+1)
        synth=angle_cost[float(g['angle']).hex()] if g['name']=='rz' else dict(T=int(g['name'] in ('t','tdg')),clifford=0)
        t=start
        for _ in range(synth['T']):
            choices=[]
            for port,center in enumerate(centers):
                path=path_to_neighbor(center,targets[0],side);occupied=footprint(path+[targets[0]])
                begin=max(t,source_ready[port],time_for(occupied));walk=len(path)-1
                # Each SWAP = three adjacent CX (nine logical clocks).
                # Injection = CX, Z measurement, conditional S (3+1+1).
                end=begin+18*walk+5
                choices.append((end,port,begin,path,occupied,walk))
            end,port,begin,path,occupied,walk=min(choices,key=lambda v:(v[0],v[1]))
            reserve(occupied,end);source_ready[port]=end+source_rounds/distance
            events.append(dict(kind='T_delivery',gate=index,port=port,start=begin,end=end,
                path=[list(v) for v in path],target=list(targets[0]),reserved=sorted(map(list,occupied)),swaps=2*walk))
            T+=1;swaps+=2*walk;t=end
        if synth['T']:
            end=t+synth['clifford'];occupied=footprint(targets)
            if end>t:events.append(dict(kind='synthesis_Cliffords',gate=index,start=t,end=end,reserved=sorted(map(list,occupied))))
        else:
            end=t+({'cx':3,'conditional_cz':3,'mx_reset':2}.get(g['name'],1))
            occupied=nodes
            events.append(dict(kind=g['name'],gate=index,start=t,end=end,reserved=sorted(map(list,occupied))))
        reserve(occupied,end);finish.append(end)
    clocks=max(ready.values())
    return dict(logical_clocks=clocks,seconds=clocks*distance*240e-9,T_count=T,delivery_SWAPS=swaps,
        delivery_CX=3*swaps+T,array_tiles=fine*fine,buffer_patches=ports,
        physical_qubits=2*distance**2*fine*fine+1024*ports,events=events,
        source_only_clocks=ceil(T/ports)*source_rounds/distance)


def validate_events(events):
    by_node={}
    for event in events:
        assert event['end']>=event['start']
        for node in event['reserved']:by_node.setdefault(tuple(node),[]).append((event['start'],event['end']))
        if event['kind']=='T_delivery':
            path=event['path']+[event['target']]
            assert len(set(map(tuple,path)))==len(path)
            assert all(abs(a[0]-b[0])+abs(a[1]-b[1])==1 for a,b in zip(path,path[1:]))
            assert event['swaps']==2*(len(path)-2)
    for intervals in by_node.values():
        ordered=sorted(intervals)
        assert all(a[1]<=b[0]+1e-10 for a,b in zip(ordered,ordered[1:]))


def run():
    oldpath=ROOT/'results/implementation_v1/ft_stage/results.json';old=json.loads(oldpath.read_text());c=old['configuration'];rows=[]
    OUT.mkdir(parents=True,exist_ok=True)
    for ref in old['rows']:
        if ref['buffer_policy']!='bounded':continue
        ops,width=stage_operations(c['L']**2,ref['batch'],c['theta']);gates,layout=nearest_neighbor_route(ops,width)
        result=schedule(gates,layout['side'],ref['rotation_classes'],ref['ports'],c['code_distance'],c['source_interval_cycles'])
        validate_events(result['events']);assert result['T_count']==ref['T_count']
        with gzip.open(OUT/f"batch{ref['batch']}_ports{ref['ports']}.json.gz",'wt') as f:json.dump(result.pop('events'),f)
        rows.append(dict(batch=ref['batch'],ports=ref['ports'],abstract_delivery_seconds=ref['seconds'],
            space_time_qubit_seconds=result['physical_qubits']*result['seconds'],**result))
    names=['python/geometric_delivery.py','python/hwp_stage.py','results/implementation_v1/ft_stage/results.json']
    (OUT/'results.json').write_text(json.dumps(dict(rows=rows,
        scope='Explicit logical SWAP delivery and restoration, perimeter ports, bounded buffers, one empty mediator tile per adjacent patch edge; square footprint includes all gaps and corners. 1024 physical qubits per reference source outside array. No decoded physical surgery, no claim of optimal routing or complete dynamics timing.',
        source_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
    for v in rows:print(v['batch'],v['ports'],v['seconds'],v['delivery_SWAPS'],v['array_tiles'])


if __name__=='__main__':run()
