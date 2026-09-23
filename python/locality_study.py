"""Nested exact patches and directly measured differences from finite references."""
import hashlib
import json

from project_paths import ROOT
from small_dynamics import HubbardSector,rectangle

OUT=ROOT/'results/implementation_v1'


def run():
    c=json.loads((ROOT/'configs/implementation_checks.json').read_text())['locality']
    observed=(1,1)
    patches={'cross5':[(1,0),(0,1),(1,1),(2,1),(1,2)],
             '2x2':rectangle(2,2),'2x3':rectangle(2,3),'3x3':rectangle(3,3)}
    rows=[];sources=['python/locality_study.py','python/small_dynamics.py','configs/implementation_checks.json']
    for name,sites in patches.items():
        model=HubbardSector(sites)
        for tau in c['times']:
            psi=model.exact(tau)
            rows.append(dict(patch=name,tau=tau,sites=sites,observed_site=observed,
                             doublon=model.doublon(psi,observed),norm_error=float(abs(sum(abs(psi)**2)-1))))
    # Reuse completed large exact calculations; no new large state simulation.
    for folder in ('dynamics_4x4','dynamics_followup_v1'):
        for path in sorted((ROOT/'results'/folder).glob('*/*_exact_*.json')):
            data=json.loads(path.read_text());contract=data.get('contract',{})
            if contract.get('algorithm')!='exact' or contract['tau'] not in c['times']:continue
            nx,ny=contract['nx'],contract['ny'];name=f'{nx}x{ny}'
            if name not in ('4x3','4x4') or any(v['patch']==name and v['tau']==contract['tau'] for v in rows):continue
            sites=rectangle(nx,ny);D=data['result']['doublons'][sites.index(observed)]
            rows.append(dict(patch=name,tau=contract['tau'],sites=sites,observed_site=observed,
                             doublon=D,norm_error=abs(data['result']['norm_squared']-1),source=str(path.relative_to(ROOT))))
            sources.append(str(path.relative_to(ROOT)))
    comparisons=[]
    for larger in rows:
        for smaller in rows:
            if smaller['tau']!=larger['tau'] or not set(map(tuple,smaller['sites']))<set(map(tuple,larger['sites'])):continue
            delta=abs(smaller['doublon']-larger['doublon'])
            comparisons.append(dict(tau=smaller['tau'],patch=smaller['patch'],finite_reference=larger['patch'],
                retained_sites=len(smaller['sites']),difference=delta,
                passes={str(eps):delta<=eps for eps in c['spatial_allowances']},
                scope='Measured truncation error relative to this finite reference; outer-boundary error unbounded'))
    record=dict(rows=rows,comparisons=comparisons,configuration=c,
        source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources},
        scope='Common Neel product state, OBC U/J=4, observed coordinate (1,1). Even reference sizes retain a nearby boundary; no bulk saturation claim.')
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'locality.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Locality:',len(rows),'exact expectations and',len(comparisons),'nested comparisons')


if __name__=='__main__':run()
