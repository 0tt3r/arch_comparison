"""Parallel regrouping networks for the unchanged plaquette product formula."""
from dataclasses import dataclass
from itertools import permutations
from math import pi

import numpy as np

from hubbard_circuits import Operation, mode_map, snake_sites, tiles
from nisq_compilation import Block, _odd_matrix, gaussian_rzz, simplify_single_qubit


@dataclass
class Stage:
    name: str
    groups: list


def plaquette_stages(L, tau, steps, u_over_j=4):
    """Same stage order and boundary merging as build_workload()."""
    if L < 2 or L % 2 or tau <= 0 or steps < 1:
        raise ValueError('Require even L >= 2, tau > 0, steps >= 1')
    durations=[]
    dt=tau/steps
    for _ in range(steps):
        for name,t in [('onsite',dt/2),('pink',dt/2),('gold',dt),('pink',dt/2),('onsite',dt/2)]:
            if name=='gold' and L==2:continue
            if durations and durations[-1][0]==name:
                durations[-1]=(name,durations[-1][1]+t)
            else:durations.append((name,t))
    modes=mode_map(L);stages=[]
    for name,t in durations:
        if name=='onsite':
            groups=[[Block('density',(2*i,2*i+1),u_over_j*t/4)] for i in range(L*L)]
        else:
            groups=[]
            for tile in tiles(L,int(name=='gold')):
                for spin in (0,1):
                    indices=[modes[x,y,spin] for x,y in tile]
                    if len(indices)==2:groups.append([Block('hop',tuple(sorted(indices)),-t/2)])
                    else:
                        a,b,c,d=indices
                        basis=[Block('hadamard',(c,a)),Block('hadamard',(b,d))]
                        groups.append(basis+[Block('hop',tuple(sorted((b,c))),-t)]+basis)
        stages.append(Stage(name,groups))
    return stages


class PairEmitter:
    """Persistent fermionic order, with pair fusion through disjoint operations."""
    def __init__(self,modes):
        self.at=list(modes);self.position={m:i for i,m in enumerate(modes)}
        self.nodes=[];self.last=[None]*len(modes);self.swaps=0;self.fusions=0

    def append(self,pair,matrix=None,density=0.):
        pair=tuple(sorted(pair));previous=self.last[pair[0]]
        if previous is not None and previous==self.last[pair[1]]:
            node=self.nodes[previous]
            if matrix is not None:node['matrix']=matrix@node['matrix']
            node['density']+=density;self.fusions+=1
        else:
            self.last[pair[0]]=self.last[pair[1]]=len(self.nodes)
            self.nodes.append({'pair':pair,'matrix':np.eye(2,dtype=complex) if matrix is None else matrix,
                               'density':density})

    def swap(self,left):
        self.append((left,left+1),np.array([[0,1],[1,0]],dtype=complex))
        self.at[left],self.at[left+1]=self.at[left+1],self.at[left]
        self.position[self.at[left]]=left;self.position[self.at[left+1]]=left+1
        self.swaps+=1

    def regroup(self,target):
        """Odd-even transposition sorting; every swap is an adjacent fermionic swap."""
        rank={m:i for i,m in enumerate(target)};layers=0
        while self.at!=target:
            changed=False
            for parity in (0,1):
                pairs=[i for i in range(parity,len(target)-1,2)
                       if rank[self.at[i]]>rank[self.at[i+1]]]
                for i in pairs:self.swap(i)
                if pairs:layers+=1;changed=True
            if not changed:raise AssertionError('Sorting network stalled')
        return layers

    def block(self,block):
        a,b=sorted(self.position[m] for m in block.modes)
        # Within a four-mode tile this takes at most two swaps; all other tiles
        # operate on disjoint physical intervals and can run simultaneously.
        for left in range(b-1,a,-1):self.swap(left)
        pair=tuple(self.position[m] for m in block.modes)
        if block.kind=='density':self.append(pair,density=block.theta)
        else:self.append(pair,_odd_matrix(block,self.position))

    def native(self,occupied,initial_order):
        gates=[Operation('x',(i,)) for i,m in enumerate(initial_order) if m in occupied]
        for node in self.nodes:
            gates.extend(gaussian_rzz(node['pair'],node['matrix']))
            if abs(node['density'])>1e-12:
                gates.append(Operation('rzz',node['pair'],2*node['density']))
        return simplify_single_qubit(gates,len(initial_order))


def _group_order(groups,position,policy):
    supports=[sorted({m for b in group for m in b.modes}) for group in groups]
    key=(lambda s:min(position[m] for m in s)) if policy=='left' else (lambda s:sum(position[m] for m in s)/len(s))
    return sorted(zip(groups,supports),key=lambda item:(key(item[1]),item[1]))


def observable_stages(stages,observed):
    """Exact for an initial occupation state and a final occupation observable.

    The initial onsite stage is a scalar on the initial state; the final onsite
    stage commutes with the measured observable. Then use even-block causality.
    """
    remaining=stages[1:] if stages and stages[0].name=='onsite' else stages[:]
    if remaining and remaining[-1].name=='onsite':remaining=remaining[:-1]
    active=set(observed);output=[]
    for stage in reversed(remaining):
        groups=[]
        for group in reversed(stage.groups):
            kept=[]
            for block in reversed(group):
                if active.intersection(block.modes):
                    kept.append(block);active.update(block.modes)
            if kept:groups.append(kept[::-1])
        if groups:output.append(Stage(stage.name,groups[::-1]))
    return output[::-1],sorted(active)


def parallel_route(stages,modes,occupied,group_order='mean',optimize_inside=False):
    emitter=PairEmitter(modes);history=[]
    for stage in stages:
        before=emitter.swaps
        groups=_group_order(stage.groups,emitter.position,group_order)
        used={m for _,support in groups for m in support}
        # Unused modes are singleton groups, keeping the route a permutation.
        groups += [([], [m]) for m in modes if m not in used]
        groups.sort(key=lambda item:(sum(emitter.position[m] for m in item[1])/len(item[1])
                                    if group_order=='mean' else min(emitter.position[m] for m in item[1])))
        target=[]
        for blocks,support in groups:
            if optimize_inside and blocks:
                offset=len(target)
                def score(order):
                    trial=list(order);cost=0
                    for block in blocks:
                        lo,hi=sorted(trial.index(m) for m in block.modes)
                        cost+=hi-lo-1
                        for i in range(hi-1,lo,-1):trial[i],trial[i+1]=trial[i+1],trial[i]
                    motion=sum(abs(emitter.position[m]-(offset+i)) for i,m in enumerate(order))
                    return motion+2*cost,order
                target.extend(min(permutations(support),key=score))
            else:target.extend(sorted(support,key=emitter.position.get))
        layers=emitter.regroup(target)
        routed=emitter.swaps-before
        for blocks,_ in groups:
            for block in blocks:emitter.block(block)
        history.append({'stage':stage.name,'regroup_swap_layers':layers,'regroup_swaps':routed,
                        'internal_swaps':emitter.swaps-before-routed,'parallel_groups':len(stage.groups)})
    return emitter.native(occupied,modes),emitter.position,{
        'initial_mode_order':list(modes),'fermionic_swaps':emitter.swaps,'pair_fusions':emitter.fusions,
        'stage_history':history,'group_order':group_order,'optimize_inside':optimize_inside,
        'routing':'parallel odd-even fermionic regrouping on the nearest-neighbor snake'}
