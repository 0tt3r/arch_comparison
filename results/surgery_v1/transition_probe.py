"""Circuit-level rotated-patch surgery calibration; isolated surgery_v1 outputs.

Ideal external Bell references and terminal Pauli measurements are a channel
test harness, not physical operations credited to the surgery implementation.
All syndrome extraction uses surface_law_v1's dressed-CZ noise convention.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import socket
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pymatching
import scipy
from scipy.optimize import minimize, brentq
from scipy.stats import beta, chi2
import stim

ROOT = Path('/lcrc/project/QC_NA_Sim/architecture_comparison')
sys.path.insert(0,str(ROOT/'python'))
OUT = ROOT / 'results/surgery_v1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def interval(k, n, alpha=.05):
    return [0. if not k else float(beta.ppf(alpha/2, k, n-k+1)),
            1. if k == n else float(beta.ppf(1-alpha/2, k+1, n-k))]


def bits(n):
    while n:
        t = n & -n
        yield t.bit_length()-1
        n ^= t


def xor_masks(items):
    out = 0
    for item in items:
        out ^= item
    return out


class Knowledge:
    """CSS stabilizers with symbolic measurement signs, excluding logical input.

    Updating on an anticommuting measurement preserves all commuting products.
    A redundant measurement supplies one local detector. Pruning at a split
    retains only the code stabilizers, not the just-measured logical parity.
    """
    def __init__(self):
        self.rows = {'X': [], 'Z': []}

    def express(self, basis, support, coefficients=False):
        pivots = {}
        for index, (s, r) in enumerate(self.rows[basis]):
            coeff = 1 << index
            while s:
                k = s.bit_length()-1
                if k not in pivots:
                    pivots[k] = (s, r, coeff)
                    break
                t, u, v = pivots[k]
                s ^= t
                r ^= u
                coeff ^= v
        r = coeff = 0
        while support:
            k = support.bit_length()-1
            if k not in pivots:
                return None
            t, u, v = pivots[k]
            support ^= t
            r ^= u
            coeff ^= v
        return (r, coeff) if coefficients else r

    def measure(self, basis, support, record):
        known = self.express(basis, support, coefficients=True)
        if known is not None:
            old, coeff = known
            j = next((i for i in bits(coeff) if self.rows[basis][i][0] == support),
                     next(bits(coeff)))
            self.rows[basis][j] = (support, record)
            return old ^ record
        other = self.rows['Z' if basis == 'X' else 'X']
        anti = [i for i, (s, _) in enumerate(other) if (s & support).bit_count() % 2]
        if anti:
            s, r = other[anti[0]]
            for j in anti[1:]:
                other[j] = (other[j][0] ^ s, other[j][1] ^ r)
            other.pop(anti[0])
        self.rows[basis].append((support, record))
        return None


class Builder:
    def __init__(self, p):
        self.p = p
        self.c = stim.Circuit()
        self.ids = {}
        self.xy_by_q=[]
        self.know = Knowledge()
        self.detector_records = []
        self.rounds = 0
        self.active_tile_rounds = 0.
        self.extra_idle_tile_rounds = 0.
        self.observables = []
        self.observable_prefix = None
        self.geometry_checks = 0

    def q(self, xy):
        if xy not in self.ids:
            self.ids[xy] = len(self.ids)
            self.xy_by_q.append(xy)
            self.c.append('QUBIT_COORDS', [self.ids[xy]], xy)
        return self.ids[xy]

    def rectangle(self, x0, y0, w, h):
        data = {self.q((x0+2*x+1, y0+2*y+1)) for x in range(w) for y in range(h)}
        checks = []
        for x in range(w+1):
            for y in range(h+1):
                xy = (x0+2*x, y0+2*y)
                basis = 'X' if ((xy[0]+xy[1])//2) % 2 else 'Z'
                offsets = [(1,1),(-1,1),(1,-1),(-1,-1)] if basis == 'X' else [(1,1),(1,-1),(-1,1),(-1,-1)]
                neighbors = []
                for dx, dy in offsets:
                    pos = (xy[0]+dx, xy[1]+dy)
                    neighbors.append(self.ids.get(pos) if x0 < pos[0] < x0+2*w and y0 < pos[1] < y0+2*h else None)
                weight = sum(q is not None for q in neighbors)
                if weight == 4 or (weight == 2 and ((basis == 'X' and y in (0,h)) or (basis == 'Z' and x in (0,w)))):
                    checks.append((self.q(xy), basis, neighbors))
        for _, b, ns in checks:
            for _, a, ms in checks:
                assert a == b or len(set(ns)-{None} & (set(ms)-{None})) % 2 == 0
        return dict(data=data, checks=checks, rect=(x0,y0,w,h))

    def noise(self, gate, qs, p=None):
        qs = list(qs)
        if self.p and qs:
            self.c.append(gate, qs, self.p if p is None else p)

    def reset(self, qs, basis, noisy=True):
        qs = sorted(qs)
        self.c.append('RX' if basis == 'X' else 'R', qs)
        if noisy:
            self.noise('Z_ERROR' if basis == 'X' else 'X_ERROR', qs)

    def init_patch(self, patch, basis='Z', ref=False):
        self.reset(patch['data'], basis)
        self.reset([a for a,_,_ in patch['checks']], 'Z')
        for _, b, ns in patch['checks']:
            if b == basis:
                self.know.rows[b].append((sum(1 << q for q in ns if q is not None), 0))
        if ref:
            assert basis == 'Z'
            r = self.q((-10-len(self.ids), -10))
            self.reset([r], 'X', noisy=False)
            for q in bits(self.logical(patch, 'X')):
                self.c.append('CX', [r, q])
            patch['ref'] = r

    def logical(self, patch, basis):
        x0, y0, _, _ = patch['rect']
        dim = 0 if basis == 'X' else 1
        edge = (x0 if basis == 'X' else y0)+1
        return sum(1 << q for xy,q in self.ids.items() if q in patch['data'] and xy[dim] == edge)

    def detector(self, record):
        if record:
            n = self.c.num_measurements
            self.c.append('DETECTOR', [stim.target_rec(i-n) for i in bits(record)], [0, 0, self.rounds])
            self.detector_records.append(record)

    def extract(self, patches, nrounds, d):
        data = set().union(*(p['data'] for p in patches))
        checks = sorted((s for patch in patches for s in patch['checks']), key=lambda x:x[0])
        anc = {a for a,_,_ in checks}
        assert not (data & anc) and len(anc) == len(checks)
        active = data | anc
        xanc = [a for a,b,_ in checks if b == 'X']
        for _ in range(nrounds):
            data_mask=sum(1<<q for q in data)
            inactive={axis:[row for row in self.know.rows[axis] if not(row[0]&data_mask)] for axis in ('X','Z')}
            self.c.append('H', xanc)
            self.noise('DEPOLARIZE1', xanc, self.p/10)
            for slot in range(4):
                pairs = [(a, ns[slot]) if b == 'X' else (ns[slot], a) for a,b,ns in checks if ns[slot] is not None]
                flat = [q for pair in pairs for q in pair]
                assert len(set(flat)) == len(flat), ('gate collision', slot)
                for a,b in pairs:
                    xy = self.xy_by_q[a]
                    uv = self.xy_by_q[b]
                    assert abs(xy[0]-uv[0]) == abs(xy[1]-uv[1]) == 1
                self.geometry_checks += len(pairs)
                targets = [b for _,b in pairs]
                self.c.append('H', targets)
                self.noise('DEPOLARIZE1', targets, self.p/10)
                self.c.append('CZ', flat)
                self.noise('DEPOLARIZE2', flat)
                self.c.append('H', targets)
                self.noise('DEPOLARIZE1', targets, self.p/10)
                self.noise('DEPOLARIZE1', sorted(active-set(flat)))
                self.c.append('TICK')
            self.c.append('H', xanc)
            self.noise('DEPOLARIZE1', xanc, self.p/10)
            ordered = [a for a,_,_ in checks]
            start = self.c.num_measurements
            self.noise('X_ERROR', ordered)
            self.c.append('M', ordered)
            self.noise('DEPOLARIZE1', sorted(data))
            self.c.append('TICK')
            self.reset(ordered, 'Z')
            self.noise('DEPOLARIZE1', sorted(data))
            self.c.append('TICK')
            prior=[self.know.express(axis,sum(1<<q for q in ns if q is not None)) for _,axis,ns in checks]
            for j,(_,b,ns) in enumerate(checks):
                r = self.know.measure(b, sum(1<<q for q in ns if q is not None), 1 << (start+j))
                if prior[j] is not None: r=prior[j] ^ (1<<(start+j))
                if r is not None:
                    self.detector(r)
            # Rebase onto the just-measured *local* checks. Keeping a generic
            # elimination basis would turn next-round time differences into
            # artificial long detector hyperedges after a deformation.
            current=Knowledge()
            for j,(_,axis,ns) in enumerate(checks):
                current.rows[axis].append((sum(1<<q for q in ns if q is not None),1<<(start+j)))
            for axis in ('X','Z'): current.rows[axis].extend(inactive[axis])
            self.know=current
            self.rounds += 1
            self.active_tile_rounds += len(data)/(d*d)

    def merge(self, a, b, spectators, nrounds, d, basis):
        x,y,w,h = a['rect']
        u,v,s,t = b['rect']
        merged = self.rectangle(min(x,u), min(y,v), (max(x+2*w,u+2*s)-min(x,u))//2,
                                (max(y+2*h,v+2*t)-min(y,v))//2)
        bridge = merged['data'] - a['data'] - b['data']
        prep = 'X' if basis == 'Z' else 'Z'
        self.reset(bridge, prep)
        old_anc = {q for patch in [a,b] for q,_,_ in patch['checks']}
        self.reset({q for q,_,_ in merged['checks']} - old_anc, 'Z')
        for q in bridge:
            self.know.rows[prep].append((1 << q, 0))
        self.extract([merged]+spectators, nrounds, d)
        # Split via destructive bridge readout in its preparation basis.
        bridge = sorted(bridge)
        start = self.c.num_measurements
        local_before={a:(axis,sum(1<<q for q in ns if q is not None),self.know.express(axis,sum(1<<q for q in ns if q is not None))) for patch in [merged]+spectators for a,axis,ns in patch['checks']}
        self.noise('Z_ERROR' if prep == 'X' else 'X_ERROR', bridge)
        self.c.append('MX' if prep == 'X' else 'M', bridge)
        retained=set().union(*(patch['data'] for patch in [a,b]+spectators))
        self.noise('DEPOLARIZE1',sorted(retained))
        self.extra_idle_tile_rounds += len(retained)/(d*d)*100/240
        self.c.append('TICK')
        bridge_mask=sum(1<<q for q in bridge)
        readout={q:1<<(start+j) for j,q in enumerate(bridge)}
        for axis,support,rec in local_before.values():
            if axis==prep and support & bridge_mask == support:
                self.detector(rec ^ xor_masks(readout[q] for q in bits(support)))
        for j,q in enumerate(bridge):
            self.know.measure(prep, 1<<q, 1<<(start+j))
        # Retain only local split-code checks; the measured logical parity is
        # deliberately not a detector (including it would erase time-like faults).
        keep = Knowledge()
        for patch in [a,b]+spectators:
            for a_check,axis,ns in patch['checks']:
                support = sum(1<<q for q in ns if q is not None)
                r = self.know.express(axis, support)
                if a_check in local_before:
                    old_axis,old_support,old_rec=local_before[a_check]
                    difference=old_support ^ support
                    if not difference: r=old_rec
                    elif axis==prep and difference & bridge_mask == difference:
                        r=old_rec ^ xor_masks(readout[q] for q in bits(difference))
                if r is not None:
                    keep.rows[axis].append((support, r))
        self.know = keep
        self.reset(old_anc - {q for q,_,_ in merged['checks']}, 'Z')
        return bridge

    def ideal_close(self, patches):
        """Perfect local syndrome boundary; never closes a logical relation."""
        for patch in patches:
            for _,axis,ns in patch['checks']:
                support = sum(1<<q for q in ns if q is not None)
                r = self.know.express(axis, support)
                j = self.c.num_measurements
                self.mpp(axis, support)
                if r is not None:
                    self.detector(r ^ (1<<j))

    def mpp(self, basis, support):
        targets = []
        for q in bits(support):
            if targets:
                targets.append(stim.target_combiner())
            targets.append((stim.target_x if basis == 'X' else stim.target_z)(q))
        self.c.append('MPP', targets)

    def observable(self, name, basis, support):
        pauli = stim.PauliString(self.c.num_qubits)
        for q in bits(support):
            pauli[q] = basis
        if self.observable_prefix is None:
            self.observable_prefix = self.c.without_noise()
        prefix = self.observable_prefix
        solution = prefix.solve_flow_measurements([stim.Flow(output=pauli)])[0]
        assert solution is not None, (name, 'not an output stabilizer')
        assert prefix.has_flow(stim.Flow(output=pauli, measurements=solution), unsigned=True)
        # Remove detector-equivalent history to limit record supports.
        record = sum(1<<i for i in solution)
        for det in reversed(self.detector_records):
            if (record ^ det).bit_count() < record.bit_count():
                record ^= det
        j = self.c.num_measurements
        self.mpp(basis, support)
        record ^= 1<<j
        self.c.append('OBSERVABLE_INCLUDE', [stim.target_rec(i-self.c.num_measurements) for i in bits(record)], len(self.observables))
        self.observables.append(name)


def build(kind, d, p, rm=None, length=0, basis='z'):
    if kind=='cnot_core':
        rm=d if rm is None else rm
        original,meta=build('cnot',d,p,rm,0,basis)
        core=stim.Circuit();ticks=0
        start=6*(d+1); end=start+12*rm+3
        for ins in original:
            if ins.name not in ('DEPOLARIZE1','DEPOLARIZE2','X_ERROR','Z_ERROR') or start<=ticks<end:
                core.append(ins)
            if ins.name=='TICK': ticks+=1
        volume=6*rm+2*rm/d
        meta={**meta,'rounds':2*rm,'pre_rounds':0,'post_rounds':0,'probe_rounds':2*d+2*rm,
              'volume':volume,'nominal_volume':6*rm,'volume_with_split_idle':volume+8*100/240,
              'harness':meta['harness']+'; pre/post storage is noiseless to isolate the two-merge core'}
        return core,meta
    if kind == 'memory':
        from surface_law_calibration import circuit
        c = circuit(d, rm, p, basis)
        return c, dict(observables=['memory'], rounds=rm, volume=rm, nominal_volume=rm,
                       harness='original noisy product-state preparation and destructive readout')
    rm = d if rm is None else rm
    b = Builder(p)
    if kind in ('idle_control','cnot_memory'):
        n=3 if kind=='cnot_memory' else 2+length
        patches=[b.rectangle(2*i*(d+1),0,d,d) for i in range(n)]
        for patch in patches: b.init_patch(patch,ref=True)
        b.c=b.c.without_noise(); b.p=0; b.extract(patches,1,d)
        b.p=p; b.rounds=0; b.active_tile_rounds=0.
        b.extract(patches if kind=='cnot_memory' else patches[:2],d,d)
        b.extract(patches,rm,d)
        if kind=='cnot_memory':
            for _ in range(2):
                b.noise('DEPOLARIZE1',sorted(set().union(*(q['data'] for q in patches))))
                b.extra_idle_tile_rounds += 3*100/240; b.c.append('TICK')
                if _==0: b.extract(patches,rm,d)
        b.noise('DEPOLARIZE1',sorted(set().union(*(q['data'] for q in patches[:2]))))
        b.extra_idle_tile_rounds += 2*100/240
        b.c.append('TICK')
        b.extract(patches[:2],d,d)
        b.ideal_close(patches)
        for j,patch in enumerate(patches):
            for axis in ('Z','X'):
                b.observable(f'patch{j}_{axis}',axis,b.logical(patch,axis)|(1<<patch['ref']))
        return b.c,dict(observables=b.observables,rounds=b.rounds,volume=b.active_tile_rounds,
                       nominal_volume=b.active_tile_rounds,volume_with_split_idle=b.active_tile_rounds+b.extra_idle_tile_rounds,
                       harness='matched ideal encoded Bell input/output boundary; independent patches',
                       pre_rounds=d,post_rounds=d,split_readout_extra_ns=300 if kind=='cnot_memory' else 100,syndrome_round_ns=240,
                       geometry_checks=b.geometry_checks)
    a = b.rectangle(0,0,d,d)
    gap = max(1, length*d)
    z = b.rectangle(0,2*(d+gap),d,d)
    if kind == 'parity':
        for patch in [a,z]:
            b.init_patch(patch, ref=True)
        b.c=b.c.without_noise()
        b.p=0
        b.extract([a,z],1,d)
        b.p=p
        b.rounds=0; b.active_tile_rounds=0.
        b.extract([a,z], d, d)
        b.merge(a,z,[],rm,d,'Z')
        b.extract([a,z], d, d)
        b.ideal_close([a,z])
        refs = (1<<a['ref']) | (1<<z['ref'])
        b.observable('parity', 'Z', refs)
        b.observable('data_a_bit', 'Z', b.logical(a,'Z') | (1<<a['ref']))
        b.observable('data_b_bit', 'Z', b.logical(z,'Z') | (1<<z['ref']))
        b.observable('data_phase', 'X', b.logical(a,'X') | b.logical(z,'X') | refs)
        nominal = 4*d+(2+length)*rm
    elif kind == 'cnot':
        # C above A, with T directly right of A: ZZ then XX.
        control, ancilla = a,z
        target = b.rectangle(2*(d+1),2*(d+1),d,d)
        b.init_patch(control, ref=True)
        b.init_patch(ancilla, basis='X')
        b.init_patch(target, ref=True)
        b.c=b.c.without_noise()
        b.p=0
        b.extract([control,ancilla,target],1,d)
        b.p=p
        b.rounds=0; b.active_tile_rounds=0.
        b.extract([control,ancilla,target],d,d)
        b.merge(control,ancilla,[target],rm,d,'Z')
        b.merge(ancilla,target,[control],rm,d,'X')
        # Destructive logical Z readout of the ancilla; frame corrections are
        # derived from exact stabilizer flows rather than guessed record signs.
        qs=sorted(ancilla['data']); start=b.c.num_measurements
        b.noise('X_ERROR',qs); b.c.append('M',qs)
        b.noise('DEPOLARIZE1',sorted(control['data']|target['data']))
        b.extra_idle_tile_rounds += 2*100/240
        b.c.append('TICK')
        for j,q in enumerate(qs):
            rec=b.know.measure('Z',1<<q,1<<(start+j))
            if rec is not None: b.detector(rec)
        keep=Knowledge()
        for patch in [control,target]:
            for _,axis,ns in patch['checks']:
                support=sum(1<<q for q in ns if q is not None)
                rec=b.know.express(axis,support)
                if rec is not None: keep.rows[axis].append((support,rec))
        b.know=keep
        b.extract([control,target],d,d)
        b.ideal_close([control,target])
        b.observable('control_z','Z', b.logical(control,'Z') | (1<<control['ref']))
        b.observable('control_x','X', b.logical(control,'X') | b.logical(target,'X') | (1<<control['ref']))
        b.observable('target_z','Z', b.logical(control,'Z') | b.logical(target,'Z') | (1<<target['ref']))
        b.observable('target_x','X', b.logical(target,'X') | (1<<target['ref']))
        nominal=5*d+6*rm
    else:
        raise ValueError(kind)
    return b.c, dict(observables=b.observables, rounds=b.rounds, volume=b.active_tile_rounds,
                     volume_with_split_idle=b.active_tile_rounds+b.extra_idle_tile_rounds,
                     nominal_volume=nominal, gap_data_columns=gap, geometry_checks=b.geometry_checks,
                     harness='ideal encoded Bell inputs and ancilla; terminal commuting Pauli readout; noisy pre/post storage and bridge preparation/readout',
                     split_readout_extra_ns=100*(1 if kind=='parity' else 3),
                     syndrome_round_ns=240, pre_rounds=d, post_rounds=d)


def inspect_case(kind,d,p,rm=None,length=0,basis='z',shots=1024):
    t=time.monotonic()
    c,meta=build(kind,d,p,rm,length,basis)
    ideal=c.without_noise()
    det,obs=ideal.compile_detector_sampler(seed=917).sample(128,separate_observables=True)
    assert not det.any() and not obs.any(), ('nondeterministic',kind,d,length,rm)
    dem=c.detector_error_model(decompose_errors=True)
    decoder=pymatching.Matching.from_detector_error_model(dem)
    distance=len(c.shortest_graphlike_error())
    det,obs=c.compile_detector_sampler(seed=918).sample(shots,separate_observables=True,bit_packed=True)
    pred=decoder.decode_batch(det,bit_packed_shots=True,bit_packed_predictions=True)
    return dict(kind=kind,d=d,p=p,rm=rm,length=length,meta=meta,seconds=time.monotonic()-t,
                detectors=c.num_detectors,qubits=c.num_qubits,graphlike_distance=distance,
                pilot_shots=shots,failures=int(np.count_nonzero(np.any(obs!=pred,axis=1))))


def cases():
    rows=[]
    for p in (.001,.0001):
        for d in ((3,5,7,9) if p==.001 else (3,5,7)):
            shots=(1000000 if d<9 else 4000000) if p==.001 else (4000000 if d==3 else 8000000)
            for length in (0,1,3):
                for rm in (d-2,d,d+2):
                    for kind in ('parity','idle_control'):
                        rows.append(dict(kind=kind,d=d,p=p,rm=rm,length=length,basis='z',shots=shots))
            if d<9:
                rows.append(dict(kind='cnot',d=d,p=p,rm=d,length=0,basis='z',shots=shots*2))
            for axis in ('x','z'):
                rows.append(dict(kind='memory',d=d,p=p,rm=10*d,length=0,basis=axis,shots=shots))
    for p in (.001,.0001):
        for d in (3,5,7):
            shots=2000000 if p==.001 else (8000000 if d==3 else 16000000)
            for kind in ('cnot_memory','cnot_core'):
                rows.append(dict(kind=kind,d=d,p=p,rm=d,length=0,basis='z',shots=shots))
    for j,row in enumerate(rows):
        row['id']=f"{row['kind']}_d{row['d']}_r{row['rm']}_l{row['length']}_p{row['p']:g}_{row['basis']}"
        row['seed']=302609240000+j*10000019
        row['batch']=25000
    return rows


def prepare():
    path=OUT/'manifest.json'
    assert not path.exists(), 'Preserve the existing manifest; use run/audit.'
    (OUT/'input_sources').mkdir(parents=True,exist_ok=True)
    for rel in ('python/surgery_calibration.py','python/surface_law_calibration.py','results/surface_law_v1/fit.json'):
        (OUT/'input_sources'/Path(rel).name).write_bytes((ROOT/rel).read_bytes())
    manifest=dict(schema=1,cases=cases(),source_sha256=sha(__file__),
        memory_source_sha256=sha(ROOT/'python/surface_law_calibration.py'),
        memory_fit_sha256=sha(ROOT/'results/surface_law_v1/fit.json'),
        python=sys.version,stim=stim.__version__,pymatching=pymatching.__version__,numpy=np.__version__,scipy=scipy.__version__,
        noise='Exact surface_law_v1: Pauli DEPOLARIZE2(p), H DEPOLARIZE1(p/10), prep/readout flip p; idle p per CZ and measurement/reset macro-slot',
        decoder='ordinary PyMatching MWPM; correlations disabled; same matching predicts all observable bits',
        encoded_boundary='Ideal Bell reference preparation and ideal joint terminal Pauli readout; original noisy SPAM memory sampled separately',
        geometry='ZZ vertical, XX horizontal; four hook-safe CZ slots. l=0 minimal one-data-column seam, l=1/3 gap of l*d data columns.',
        volume='Area-equivalent data tiles integrated over active extraction rounds; extra 100 ns split idles reported separately. l=0 adds R_m/d tile-rounds to nominal two-patch volume.',
        lifetime='R_pre=R_post=d. Bridges live only during merge. CNOT ancilla lives through both merges.',
        analysis='Cell intervals exact Clopper-Pearson; ratios conditional on fixed saved memory fit. Fits use binomial likelihood and likelihood-ratio profile intervals. <30 failures flagged.',
        compute_budget='One backfill allocation, select=1:ncpus=128, <=3 h; fixed-shot cells, no failure-count stopping.',
        shots=sum(x['shots'] for x in cases()))
    atomic(path,manifest)
    print('Prepared',len(manifest['cases']),'cells,',manifest['shots'],'shots',flush=True)


def schedule_check():
    """Match every noisy gate location and every CZ slot to the memory helper."""
    from collections import Counter
    from surface_law_calibration import circuit
    evidence=[]
    for d in (3,5,7,9):
        b=Builder(.001);patch=b.rectangle(0,0,d,d);b.init_patch(patch);b.extract([patch],3,d)
        original=circuit(d,3,.001,'z'); oc=original.get_final_qubit_coordinates()
        prefix=stim.Circuit()
        for ins in original:
            ts=[t.value for t in ins.targets_copy() if t.is_qubit_target]
            if ins.name=='X_ERROR' and set(ts)=={q for q,xy in oc.items() if xy[0]%2 and xy[1]%2} and prefix.num_measurements:
                break
            prefix.append(ins)
        def inventory(c):
            coords=c.get_final_qubit_coordinates(); counts=Counter(); slots=[]
            for ins in c:
                if ins.name in ('QUBIT_COORDS','DETECTOR','OBSERVABLE_INCLUDE','TICK','SHIFT_COORDS'): continue
                qs=[t.value for t in ins.targets_copy()]
                if ins.name in ('CZ','DEPOLARIZE2'):
                    pairs=[tuple(sorted((tuple(coords[a]),tuple(coords[z])))) for a,z in zip(qs[::2],qs[1::2])]
                    for pair in pairs: counts[(ins.name,tuple(ins.gate_args_copy()),pair)]+=1
                    if ins.name=='CZ': slots.append(sorted(pairs))
                else:
                    for q in qs: counts[(ins.name,tuple(ins.gate_args_copy()),tuple(coords[q]))]+=1
            return counts,slots
        assert inventory(prefix)==inventory(b.c),(d,'memory extraction schedule differs')
        evidence.append(dict(d=d,rounds=3,noise_locations_equal=True,CZ_slots_equal=True))
    atomic(OUT/'schedule_validation.json',dict(status='pass',checks=evidence))


def filter_observables(c,keep):
    out=stim.Circuit()
    for ins in c:
        if ins.name=='OBSERVABLE_INCLUDE':
            old=int(ins.gate_args_copy()[0])
            if old in keep: out.append('OBSERVABLE_INCLUDE',ins.targets_copy(),keep.index(old))
        else: out.append(ins)
    return out


def physical_logical_flip_check(c,meta,d):
    if meta.get('harness','').startswith('original'): return True
    coords=c.get_final_qubit_coordinates()
    # X string on the first patch's left edge, just before ideal terminal readout.
    support=[q for q,xy in coords.items() if xy[0]==1 and 0<xy[1]<2*d and xy[1]%2==1]
    out=stim.Circuit(); inserted=False
    for ins in c.without_noise():
        if ins.name=='MPP' and not inserted:
            out.append('X_ERROR',support,1); inserted=True
        out.append(ins)
    det,obs=out.compile_detector_sampler(seed=7301).sample(64,separate_observables=True)
    assert inserted and not det.any() and np.all(np.any(obs,axis=1)), 'logical flip control failed'
    return True


def single_fault_check(c,decoder):
    entries=[]
    for ins in c.detector_error_model(decompose_errors=False).flattened():
        if ins.type=='error':
            entries.append(([t.val for t in ins.targets_copy() if t.is_relative_detector_id()],
                            [t.val for t in ins.targets_copy() if t.is_logical_observable_id()]))
    bad=0
    for j in range(0,len(entries),2048):
        es=entries[j:j+2048]
        det=np.zeros((len(es),c.num_detectors),dtype=np.uint8)
        obs=np.zeros((len(es),c.num_observables),dtype=np.uint8)
        for i,(ds,ls) in enumerate(es):det[i,ds]=1;obs[i,ls]=1
        bad+=int(np.count_nonzero(np.any(decoder.decode_batch(det)!=obs,axis=1)))
    return dict(mechanisms=len(entries),misdecoded=bad,
                scope='all single undecomposed detector-error-model mechanisms; not an exhaustive multi-fault distance proof')


def sample_case(task):
    row,group,deadline=task
    started=time.monotonic()
    manifest=json.loads((OUT/'manifest.json').read_text())
    assert sha(__file__)==manifest['source_sha256'],'source changed during campaign'
    folder=OUT/group; folder.mkdir(parents=True,exist_ok=True)
    path=folder/(row['id']+'.json')
    if path.exists():
        old=json.loads(path.read_text())
        assert old['case']==row and old['source_sha256']==manifest['source_sha256']
        if old['status']=='complete': return row['id'],'already complete'
    lock=folder/(row['id']+'.lock')
    with lock.open('x') as f: json.dump(dict(host=socket.gethostname(),pid=os.getpid(),job=os.getenv('PBS_JOBID')),f)
    try:
        c,meta=build(row['kind'],row['d'],row['p'],row['rm'],row['length'],row['basis'])
        ideal=c.without_noise()
        det,obs=ideal.compile_detector_sampler(seed=917).sample(128,separate_observables=True)
        assert not det.any() and not obs.any()
        dem=c.detector_error_model(decompose_errors=True)
        decoder=pymatching.Matching.from_detector_error_model(dem)
        distance=len(c.shortest_graphlike_error())
        expected=min(row['d'],row['rm']) if row['kind']=='parity' else row['d']
        assert distance==expected,(row['id'],distance,expected)
        physical_logical_flip_check(c,meta,row['d'])
        single=single_fault_check(c,decoder)
        assert expected==1 or single['misdecoded']==0,(row['id'],single)
        circuit_path=folder/(row['id']+'.stim')
        circuit_path.write_text(str(c)+'\n')
        dem_path=folder/(row['id']+'.dem'); dem_path.write_text(str(dem)+'\n')
        hist=np.zeros(1<<c.num_observables,dtype=np.int64); batch_rows=[]; shots=0
        progress=folder/(row['id']+'.progress.json')
        if progress.exists():
            saved=json.loads(progress.read_text())
            assert saved['case']==row and saved['circuit_sha256']==sha(circuit_path)
            hist=np.array(saved['histogram'],dtype=np.int64); shots=saved['shots']; batch_rows=saved['batches']
        for j in range(len(batch_rows),math.ceil(row['shots']/row['batch'])):
            if time.time()>deadline: break
            n=min(row['batch'],row['shots']-shots)
            sampler=c.compile_detector_sampler(seed=row['seed']+j*100003)
            det,obs=sampler.sample(n,separate_observables=True,bit_packed=True)
            pred=decoder.decode_batch(det,bit_packed_shots=True,bit_packed_predictions=True)
            packed=obs^pred
            codes=np.zeros(n,dtype=np.uint32)
            for k in range(packed.shape[1]): codes |= packed[:,k].astype(np.uint32)<<(8*k)
            counts=np.bincount(codes,minlength=len(hist)); hist+=counts; shots+=n
            batch_rows.append(dict(shots=n,failures=int(n-counts[0]),histogram_sha256=hashlib.sha256(counts.tobytes()).hexdigest()))
            if len(batch_rows)%20==0 or shots==row['shots']:
                result=dict(status='complete' if shots==row['shots'] else 'running',case=row,meta=meta,shots=shots,
                    histogram=hist.tolist(),batches=batch_rows,source_sha256=manifest['source_sha256'],
                    circuit_sha256=sha(circuit_path),dem_sha256=sha(dem_path),graphlike_distance=distance,
                    single_fault_check=single,
                    noiseless_check=True,logical_flip_check=True,detectors=c.num_detectors,
                    elapsed_seconds=time.monotonic()-started,host=socket.gethostname(),job=os.getenv('PBS_JOBID'))
                atomic(path if shots==row['shots'] else progress,result)
        if shots<row['shots']:
            raise RuntimeError(f"Compute cutoff: incomplete fixed-shot cell {row['id']} ({shots}/{row['shots']})")
        return row['id'],dict(shots=shots,failures=int(shots-hist[0]),seconds=time.monotonic()-started)
    finally:
        lock.unlink()


def run(workers,smoke=False,hours=2.65):
    manifest=json.loads((OUT/'manifest.json').read_text())
    rows=manifest['cases']
    if smoke:
        schedule_check()
        wanted=[('parity',3,1,0,.001),('parity',7,7,3,.001),('cnot',5,5,0,.001),('idle_control',5,5,1,.001),('cnot_memory',3,3,0,.001),('cnot_core',3,3,0,.001)]
        rows=[{**r,'shots':25000} for r in rows if (r['kind'],r['d'],r['rm'],r['length'],r['p']) in wanted]
        assert len(rows)==6
    deadline=time.time()+hours*3600
    rows=sorted(rows,key=lambda r:hashlib.sha256(r['id'].encode()).hexdigest())
    failures=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(sample_case,(r,'smoke' if smoke else 'cases',deadline)):r['id'] for r in rows}
        for f in as_completed(futures):
            try: print(f.result(),flush=True)
            except Exception as e:
                failures.append(dict(case=futures[f],error=repr(e))); print('FAILED',failures[-1],flush=True)
    atomic(OUT/('smoke_status.json' if smoke else 'run_status.json'),dict(cases=len(rows),failures=failures,status='pass' if not failures else 'incomplete',job=os.getenv('PBS_JOBID')))
    assert not failures,failures


def scalar_scale_fit(rows):
    n=np.array([r['shots'] for r in rows]); k=np.array([r['failures'] for r in rows]); q0=np.array([r['prediction'] for r in rows])
    from scipy.special import xlogy
    def nll(x):
        q=-np.expm1(-math.exp(float(x))*q0)
        return float(-np.sum(xlogy(k,q)+xlogy(n-k,1-q)))
    optimum=minimize(lambda x:nll(x[0]),[0.],method='Nelder-Mead')
    x=float(optimum.x[0]); best=nll(x); target=chi2.ppf(.95,1)/2
    lo=brentq(lambda y:nll(y)-best-target,x-30,x) if sum(k)>0 else None
    hi=brentq(lambda y:nll(y)-best-target,x,x+20)
    sat=float(-np.sum(xlogy(k,k/n)+xlogy(n-k,1-k/n)))
    dev=max(0,2*(best-sat)); df=len(rows)-1
    return dict(multiplier=math.exp(x) if sum(k)>0 else 0.,profile95=[0. if lo is None else math.exp(lo),math.exp(hi)],
                deviance=dev,df=df,gof_p=float(chi2.sf(dev,df)) if df else None,
                law='q = 1-exp(-multiplier * V * p_L(d))',failures=int(sum(k)),shots=int(sum(n)),
                warning='descriptive common multiplier; lack of fit means use per-cell ratios')


def time_fit(rows):
    """Exact scalar profile: solve each binomial intercept before the slope."""
    from scipy.special import xlogy
    from scipy.optimize import minimize_scalar
    keys=sorted(set((r['d'],r['length']) for r in rows))
    groups=[]
    for key in keys:
        rr=[r for r in rows if (r['d'],r['length'])==key]
        groups.append((key,np.array([(r['rm']+1)/2 for r in rr]),
                       np.array([r['failures'] for r in rr]),np.array([r['shots'] for r in rr])))
    failures=sum(r['failures'] for r in rows)
    if failures<30: return dict(status='insufficient_failures',failures=failures)
    def profile(t, details=False):
        total=0.; amplitudes={}
        for key,h,k,n in groups:
            if not k.sum():
                amplitudes[f'd{key[0]}_l{key[1]}']=0.
                continue
            def score(a):
                q=np.exp(a-h*t)
                return float(np.sum((n*q-k)/(1-q)))
            a=brentq(score,-120.,float(min(h*t)-1e-10),xtol=1e-11)
            q=np.exp(a-h*t)
            total-=float(np.sum(xlogy(k,q)+xlogy(n-k,1-q)))
            amplitudes[f'd{key[0]}_l{key[1]}']=float(np.exp(a))
        return (total,amplitudes) if details else total
    opt=minimize_scalar(profile,bounds=(0.,15.),method='bounded',options={'xatol':1e-9})
    t=float(opt.x); best,amplitudes=profile(t,True); target=chi2.ppf(.95,1)/2
    objective=lambda z:profile(z)-best-target
    low=brentq(objective,0.,t) if objective(0.)>0 else 0.
    high=brentq(objective,t,15.) if objective(15.)>0 else None
    sat=0.
    for _,h,k,n in groups: sat-=float(np.sum(xlogy(k,k/n)+xlogy(n-k,1-k/n)))
    dev=max(0.,2*(best-sat)); df=len(rows)-sum(bool(k.sum()) for _,h,k,n in groups)-1
    return dict(status='fit',Lambda_t=float(np.exp(t)),Lambda_t_profile95=[float(np.exp(low)),None if high is None else float(np.exp(high))],
                B_by_d_length=amplitudes,deviance=dev,df=df,gof_p=float(chi2.sf(dev,df)) if df>0 else None,
                optimizer_success=bool(opt.success),failures=failures,
                law='q_time_only(d,l,R_m) = B_d,l * Lambda_t^(-(R_m+1)/2)',
                warning='outcome-only errors can include spatial mechanisms; profile intervals and chi-square GOF are asymptotic, conditional on this form; zero-count amplitudes are boundary estimates')


def analyze():
    manifest=json.loads((OUT/'manifest.json').read_text()); law=json.loads((OUT/'input_sources/fit.json').read_text())['fits']
    raw=[]; complete=[]
    for case in manifest['cases']:
        path=OUT/'cases'/(case['id']+'.json')
        if not path.exists(): continue
        row=json.loads(path.read_text())
        assert row['status']=='complete' and row['shots']==case['shots']
        complete.append(row); hist=np.array(row['histogram']); labels=np.arange(len(hist)); n=row['shots']
        rates={'total':int(n-hist[0])}
        for j,name in enumerate(row['meta']['observables']): rates[name]=int(hist[(labels&(1<<j))!=0].sum())
        if case['kind']=='parity':
            rates['space']=int(hist[(labels>>1)!=0].sum())
            rates['time_only']=int(hist[labels==1].sum())
        fit=law[str(case['p'])]; pl=fit['A_patch_union_bound']*fit['Lambda']**(-(case['d']+1)/2)
        for metric,k in rates.items():
            lo,hi=interval(k,n); pred=row['meta']['volume']*pl
            raw.append(dict(**case,metric=metric,failures=k,rate=k/n,lo95=lo,hi95=hi,fewer_than_30=k<30,
                volume=row['meta']['volume'],nominal_volume=row['meta']['nominal_volume'],
                volume_with_split_idle=row['meta'].get('volume_with_split_idle',row['meta']['volume']),
                prediction=pred,nominal_prediction=row['meta']['nominal_volume']*pl,nominal_ratio=(k/n)/(row['meta']['nominal_volume']*pl),ratio=(k/n)/pred,ratio_lo95=lo/pred,ratio_hi95=hi/pred,
                circuit_sha256=row['circuit_sha256'],graphlike_distance=row['graphlike_distance']))
    assert len(complete)==len(manifest['cases']),(len(complete),len(manifest['cases']))
    path=OUT/'raw.csv'
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(raw[0])); writer.writeheader(); writer.writerows(raw)
    temporal={str(p):time_fit([r for r in raw if r['kind']=='parity' and r['metric']=='time_only' and r['p']==p]) for p in (.001,.0001)}
    local_temporal={}
    for p in (.001,.0001):
        for d in ((3,5,7,9) if p==.001 else (3,5,7)):
            for l in (0,1,3):
                local_temporal[f'p{p:g}_d{d}_l{l}']=time_fit([r for r in raw if r['kind']=='parity' and r['metric']=='time_only' and r['p']==p and r['d']==d and r['length']==l])
    scales={}
    for p in (.001,.0001):
        for kind in ('parity','cnot','idle_control','cnot_memory','cnot_core'):
            for l in ((0,) if kind.startswith('cnot') else (0,1,3)):
                rs=[r for r in raw if r['p']==p and r['kind']==kind and r['length']==l and r['rm']==r['d'] and r['metric']=='total']
                scales[f'{kind}_l{l}_p{p:g}']=scalar_scale_fit(rs)
    crossover=[]
    for p in (.001,.0001):
        for d in ((3,5,7,9) if p==.001 else (3,5,7)):
            for l in (0,1,3):
                rs=[r for r in raw if r['p']==p and r['d']==d and r['length']==l and r['kind']=='parity']
                by={(r['rm'],r['metric']):r for r in rs}; point=[]; separated=[]
                for rm in (d-2,d,d+2):
                    t=by[(rm,'time_only')]; s=by[(rm,'space')]
                    if t['rate']<s['rate']: point.append(rm)
                    if t['hi95']<s['lo95']: separated.append(rm)
                crossover.append(dict(p=p,d=d,length=l,smallest_tested_point_estimate=min(point) if point else None,
                                      smallest_tested_separated_intervals=min(separated) if separated else None,
                                      note='pointwise 95% intervals; no extrapolation below/above tested R_m'))
    # Original SPAM controls: compare against the saved *basis-specific*,
    # endpoint-corrected probability, not an any-Pauli union bound.
    controls=[]
    for r in raw:
        if r['kind']!='memory' or r['metric']!='total': continue
        f=law[str(r['p'])]; q=f['A_'+r['basis']]*f['Lambda']**(-(r['d']+1)/2)
        endpoints=f['endpoint_corrections_rounds']
        delta=endpoints[f"p{r['p']:g}_d{r['d']}_{r['basis']}"]
        predicted=.5*(1-(1-2*q)**(r['rm']+delta))
        old=[]
        for src in (ROOT/'results/surface_law_v1/cases').glob('*.json'):
            v=json.loads(src.read_text()); cc=v.get('contract',{})
            if cc.get('p')==r['p'] and cc.get('distance')==r['d'] and cc.get('basis')==r['basis'] and cc.get('rounds')==r['rm'] and v.get('status')=='complete': old.append(v)
        old_k=sum(v['failures'] for v in old);old_n=sum(v['shots'] for v in old)
        old_ci=interval(old_k,old_n) if old_n else None
        controls.append(dict(id=r['id'],p=r['p'],d=r['d'],basis=r['basis'],measured=r['rate'],interval95=[r['lo95'],r['hi95']],
                             q_per_round=q,endpoint_rounds=delta,rounds=r['rm'],fit_probability=predicted,
                             fit_point_in_new_interval=r['lo95']<=predicted<=r['hi95'],original_shots=old_n,original_failures=old_k,
                             original_interval95=old_ci,intervals_overlap=None if old_ci is None else max(old_ci[0],r['lo95'])<=min(old_ci[1],r['hi95'])))
    result=dict(status='complete',cells=len(complete),shots=sum(r['shots'] for r in complete),
                temporal=temporal,temporal_by_geometry=local_temporal,multipliers=scales,time_below_space=crossover,memory_controls=controls,
                ratio_intervals='Numerator-only exact binomial intervals, conditional on saved memory-fit point estimates; not calibration/systematic confidence bands.',
                scope='Specific isolated four-slot rotated-patch constructions with ideal encoded input/output probes; no congestion, twists/Y surgery, correlated noise, or d>9.',
                source_sha256=sha(__file__),raw_sha256=sha(path))
    atomic(OUT/'fit.json',result)
    manifest['artifacts']={r['case']['id']:dict(circuit_sha256=r['circuit_sha256'],dem_sha256=r['dem_sha256'],shots=r['shots'],seed=r['case']['seed']) for r in complete}
    atomic(OUT/'manifest.json',manifest)
    print('Analysis:',len(complete),'fixed-shot cells;',result['shots'],'shots',flush=True)


def rebuild_check(case):
    path=OUT/'cases'/(case['id']+'.stim')
    c,meta=build(case['kind'],case['d'],case['p'],case['rm'],case['length'],case['basis'])
    assert str(c)+'\n'==path.read_text()
    ideal=c.without_noise();det,obs=ideal.compile_detector_sampler(seed=917).sample(128,separate_observables=True)
    assert not det.any() and not obs.any()
    return case['id']


def audit(replay=False):
    m=json.loads((OUT/'manifest.json').read_text())
    assert sha(__file__)==m['source_sha256']
    assert sha(ROOT/'python/surface_law_calibration.py')==m['memory_source_sha256']
    rows=[]
    for case in m['cases']:
        path=OUT/'cases'/(case['id']+'.json'); r=json.loads(path.read_text())
        assert r['status']=='complete' and r['case']==case and r['shots']==case['shots']
        assert sum(r['histogram'])==r['shots']
        assert sum(b['shots'] for b in r['batches'])==r['shots']
        assert sum(b['failures'] for b in r['batches'])==r['shots']-r['histogram'][0]
        cp=path.with_suffix('.stim'); dp=path.with_suffix('.dem')
        assert sha(cp)==r['circuit_sha256'] and sha(dp)==r['dem_sha256']
        assert r['noiseless_check'] and r['logical_flip_check']
        assert r['graphlike_distance']==1 or r['single_fault_check']['misdecoded']==0
        rows.append(case['id'])
    if replay:
        assert os.getenv('PBS_JOBID'),'Full parallel rebuild belongs on a compute node; ordinary audit is safe on login.'
        with ProcessPoolExecutor(max_workers=32) as pool:
            assert len(list(pool.map(rebuild_check,m['cases'])))==len(rows)
    fit=json.loads((OUT/'fit.json').read_text());assert fit['raw_sha256']==sha(OUT/'raw.csv')
    atomic(OUT/'audit.json',dict(status='pass',cells=len(rows),shots=sum(c['shots'] for c in m['cases']),
        checks=['fixed-shot manifest census','batch and observable histograms','source, circuit, DEM and CSV hashes','per-case zero-noise and inserted logical-flip controls','per-case shortest graphlike fault distance','all single undecomposed DEM mechanisms decoded correctly except intentional R_m=1 time boundary'],
        rebuilt_all_circuits=replay,job=os.getenv('PBS_JOBID')))
    print('AUDIT PASS',len(rows),'cells',flush=True)


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=list(csv.DictReader((OUT/'raw.csv').open()))
    fig,axs=plt.subplots(2,3,figsize=(12,7),sharex=True)
    colors=['#2487a0','#d78325','#775ca6']
    for i,p in enumerate((.001,.0001)):
        for j,kind in enumerate(('parity','cnot','idle_control')):
            ax=axs[i,j]
            for k,l in enumerate((0,) if kind=='cnot' else (0,1,3)):
                rs=sorted([r for r in rows if float(r['p'])==p and r['kind']==kind and int(r['length'])==l and r['rm']==r['d'] and r['metric']=='total'],key=lambda r:int(r['d']))
                ds=np.array([int(r['d']) for r in rs]); ys=np.array([float(r['rate']) for r in rs]);lo=np.array([float(r['lo95']) for r in rs]);hi=np.array([float(r['hi95']) for r in rs])
                for z,r in enumerate(rs):
                    marker='o' if int(r['failures'])>=30 else ('v' if int(r['failures'])==0 else 's')
                    yy=ys[z] if ys[z]>0 else hi[z]
                    ax.errorbar(ds[z]+.04*k,yy,yerr=[[max(0,yy-lo[z])],[max(0,hi[z]-yy)]],fmt=marker,color=colors[k],ms=4,capsize=2,mfc=colors[k] if int(r['failures'])>=30 else 'white')
                ax.plot(ds,[float(r['prediction']) for r in rs],'--',color=colors[k],label=('memory law' if kind=='cnot' else f'l={l}: memory law'))
            ax.set_yscale('log');ax.set_xticks([3,5,7,9]);ax.grid(alpha=.2)
            ax.set_title(f'{kind.replace("_"," ")} | p={p:g}')
            if j==0: ax.set_ylabel('Decoded gadget / window failure')
            if i==1: ax.set_xlabel('Code distance d')
            ax.legend(fontsize=7)
    fig.suptitle('Isolated lattice surgery, R_m=d: measured channel failure versus V p_L(d)')
    fig.text(.5,.01,'95% Clopper–Pearson intervals; hollow markers: <30 failures; triangles: zero-failure upper bounds. Ideal encoded input/output probes.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.035,1,.95]);fig.savefig(ROOT/'figures/surgery_law.png',dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['pilot','prepare','smoke','run','analyze','audit','plot'])
    parser.add_argument('--kind',default='parity',choices=['parity','cnot','memory','idle_control','cnot_memory','cnot_core'])
    parser.add_argument('--d',type=int,default=3)
    parser.add_argument('--p',type=float,default=.001)
    parser.add_argument('--rm',type=int)
    parser.add_argument('--length',type=int,default=0)
    parser.add_argument('--shots',type=int,default=1024)
    parser.add_argument('--workers',type=int,default=112)
    parser.add_argument('--hours',type=float,default=2.65)
    parser.add_argument('--replay',action='store_true')
    args=parser.parse_args()
    if args.command=='prepare': return prepare()
    if args.command in ('run','smoke'): return run(args.workers,args.command=='smoke',args.hours)
    if args.command=='analyze': return analyze()
    if args.command=='audit': return audit(args.replay)
    if args.command=='plot': return plot()
    row=inspect_case(args.kind,args.d,args.p,args.rm,args.length,shots=args.shots)
    atomic(OUT/'pilot'/f'{args.kind}_d{args.d}_r{args.rm}_l{args.length}.json',row)
    print(json.dumps(row),flush=True)


if __name__=='__main__':
    main()
