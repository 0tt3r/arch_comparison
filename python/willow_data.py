"""Reduce Willow shots and published Z/ZZ estimates without refitting mitigation."""
import argparse
import gzip
import io
import json
from pathlib import Path

import h5py
import numpy as np

from phasecraft_data import DataOnlyReader, columns
from project_paths import ROOT
from willow_archive import CACHE, OUT, read, save, sha


def nonempty_columns(group):
    shape = group['axis1'].attrs.get('shape')
    if shape is not None:
        shape = DataOnlyReader(io.BytesIO(shape)).load()
        if np.prod(shape)==0:
            return None
    return columns(group)


def labels_for(L):
    return [(x,y,s) for x in range(L) for y in range(L) for s in ('up','down')]


def selected_columns(c):
    """Keep Z and same-site, opposite-spin ZZ; derive doublons from these."""
    selected=[]
    for i, kind in enumerate(c['obs_type']):
        if kind=='Z' or (kind=='ZZ' and c['site_coords'][i][0]==c['site_coords'][i][1]
                        and len(set(c['spins'][i]))==2):
            selected.append(i)
    return selected


def published_rows(path):
    rows=[];empty=[]
    with h5py.File(path) as h:
        for group in h:
            c=nonempty_columns(h[group])
            if c is None:
                empty.append(group)
                continue
            for i in selected_columns(c):
                rows.append(dict(method=str(c['method'][i]),group=group,
                    tau=float(c['ev_time'][i]),kind=str(c['obs_type'][i]),
                    sites=[list(v) for v in c['site_coords'][i]],spins=list(c['spins'][i]),
                    value=float(c['obs_val'][i]),se=float(c['obs_err'][i])))
    return rows,empty


def zero_time_state(rows,L):
    zero=[v for v in rows if v['method']=='TDVP chi=2048' and v['tau']==0 and v['kind']=='Z']
    state={(*v['sites'][0],v['spins'][0]):round((1-v['value'])/2) for v in zero}
    assert set(state)==set(labels_for(L)),'Reference input must cover every physical mode'
    assert all(abs(v['value'])>1-1e-10 for v in zero),'Expected a product input state'
    assert sum(state.values())==L*L-1
    return state


def initial_state_issues(state_json,L):
    state=state_json['center_hole'][f'{L}x{L}']
    coords=state['ups']+state['downs']+state['holons']+state['doublons']
    coords += [xy for pair in state['singlet_pairs'] for xy in pair]
    return dict(out_of_bounds=sorted({tuple(xy) for xy in coords if any(t<0 or t>=L for t in xy)}),
                unexpected_singlet_pairs=len(state['singlet_pairs']),
                expected_center_hole_product_state=not state['singlet_pairs'] and len(state['holons'])==1)


def shot_values(bits,labels,L):
    up=[labels.index(f'({x},{y},up)') for x in range(L) for y in range(L)]
    down=[labels.index(f'({x},{y},down)') for x in range(L) for y in range(L)]
    u,d=bits[:,up].astype(float),bits[:,down].astype(float)
    D=u*d
    phase=np.array([(-1)**(x+y) for x in range(L) for y in range(L)])
    # Columns: all local doublons, spatial mean, staggered magnetization,
    # then single-mode Z in coordinate/spin order. Ancillas are excluded.
    z=np.stack((1-2*u,1-2*d),axis=2).reshape(len(bits),2*L*L)
    return np.column_stack((D,D.mean(1),(u-d)@phase/(2*L*L),z)),u.sum(1),d.sum(1)


def grouped_statistics(groups):
    """Variance across independent randomized circuits, preserving shot grouping.

    The cluster SE includes finite-shot fluctuations; do not add them twice.
    Between-circuit and within-circuit terms for design extrapolation are
    separately estimated with a nonnegative random-effects component.
    """
    K=len(groups)
    if K<2 or min(map(len,groups))<2:
        raise ValueError('At least two circuits and two accepted shots per circuit required')
    means=np.array([v.mean(0) for v in groups])
    variances=np.array([v.var(0,ddof=1) for v in groups])
    n=np.array([len(v) for v in groups])
    within=(variances/n[:,None]).mean(0)
    between=means.var(0,ddof=1)
    latent=np.maximum(0.,between-within)
    return dict(mean=means.mean(0).tolist(),cluster_se=np.sqrt(between/K).tolist(),
        within_mean_variance=within.tolist(),between_mean_variance=between.tolist(),
        estimated_twirl_variance=latent.tolist(),
        fixed_twirl_infinite_shots_se=np.sqrt(latent/K).tolist(),accepted_per_twirl=n.tolist())


def reduce_case(folder,L,U,flux,r,cfg,records):
    for filename in ('shots.h5','exp_vals.h5'):
        name=str((folder/filename).relative_to(CACHE))
        assert sha(folder/filename)==records[name]['sha256'],'Cached input changed'
    rows,empty=published_rows(folder/'exp_vals.h5')
    state=zero_time_state(rows,L)
    nup=sum(v for (*xy,s),v in state.items() if s=='up')
    ndown=sum(v for (*xy,s),v in state.items() if s=='down')
    shots=[]
    with h5py.File(folder/'shots.h5') as h:
        labels=list(h['metadata'].attrs['qubit_order'])
        assert h['metadata'].attrs['system_dims']==f'{L}x{L}'
        assert h['metadata'].attrs['onsite_term']==U and h['metadata'].attrs['n_layers']==r
        assert len(h['data'])==cfg['expected_times']
        for key in sorted(h['data'],key=float):
            ds=h['data'][key];shape=tuple(map(int,ds.attrs['original_shape']))
            assert shape==(cfg['twirls']*cfg['shots_per_twirl'],cfg['table_I'][str(L)]['physical_qubits'])
            bits=np.unpackbits(ds[:],axis=1)[:,:shape[1]]
            values,up,down=shot_values(bits,labels,L)
            distance=np.abs(up-nup)+np.abs(down-ndown)
            item=dict(tau=float(key),shots=shape[0],qubits=shape[1],statistics={})
            for name,limit in [('unfiltered',2*L*L),('strict',0),('relaxed1',1),('relaxed2',2)]:
                groups=[]
                for k in range(cfg['twirls']):
                    a=k*cfg['shots_per_twirl'];b=a+cfg['shots_per_twirl']
                    groups.append(values[a:b][distance[a:b]<=limit])
                stats=grouped_statistics(groups)
                stats['acceptance']=sum(map(len,groups))/shape[0]
                item['statistics'][name]=stats
            shots.append(item)
    return dict(L=L,U=U,flux=flux,steps=r,empty_groups=empty,
        initial_occupations=[dict(site=[x,y],spin=s,n=v) for (x,y,s),v in sorted(state.items())],
        observable_order='L² local doublons (x-major), mean doublon, staggered magnetization, then 2L² Z (x,y,up/down)',
        shots=shots,published=rows)


def reduce_all():
    cfg=read(ROOT/'configs/willow_benchmark.json')
    downloads=read(OUT/'downloads.json')
    assert downloads['config_sha256']==sha(ROOT/'configs/willow_benchmark.json')
    records={v['name']:v for v in downloads['files']}
    states=read(CACHE/'programmable_simulation_fermi_hubbard/initial_states.json')
    cases=[]
    for L in cfg['sizes']:
        for flux in cfg['fluxes']:
            for U in cfg['interactions']:
                for r in cfg['steps']:
                    folder=CACHE/f'programmable_simulation_fermi_hubbard/center_hole/{L}x{L}/mag_field_{flux}/U_{U}/trotter_layers_{r}'
                    result=reduce_case(folder,L,U,flux,r,cfg,records)
                    result['initial_state_json_audit']=initial_state_issues(states,L)
                    cases.append(result)
                    print(f'Reduced L={L} flux={flux} U={U} r={r}',flush=True)
    result=dict(cases=cases,inputs={str(path.relative_to(ROOT)):sha(path) for path in
        [Path(__file__),ROOT/'python/phasecraft_data.py',ROOT/'configs/willow_benchmark.json',OUT/'downloads.json']},
        source=cfg['record'],license=cfg['license'],
        scope='Published processed estimates are not refitted. U=0 exact groups depend on step count and are circuit references. Interacting exact groups are empty. TDVP is approximate. Randomized-compilation blocks are retained.')
    path=OUT/'reduced.json.gz'
    compressed=gzip.compress(json.dumps(result,separators=(',',':'),allow_nan=False).encode(),mtime=0)
    temporary=path.with_suffix('.tmp');temporary.write_bytes(compressed);temporary.replace(path)
    print('Saved',len(cases),'cases;',len(compressed),'bytes')


if __name__=='__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    reduce_all()
