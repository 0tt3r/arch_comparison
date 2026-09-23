"""Reduce the pinned experimental archive to observable and shot summaries.

No mitigation is refitted. The published U=0 result is a training consistency
check, and U=4 tensor-network comparisons are not exact-truth error estimates.
"""
import argparse
import hashlib
import io
import json
import pickle
import re
from pathlib import Path
import urllib.request
import zipfile

import h5py
import numpy as np

from project_paths import ROOT

OUT = ROOT / 'results/published_benchmarks_v1'
CACHE = ROOT / 'build/fermionic_dynamics'
URL = 'https://zenodo.org/api/records/17799843/files/fermionic_dynamics.zip/content'
MD5 = 'abc5b5ab33d6d41989578dd6bc505618'


class OpaqueUncertainty:
    """Ignore redundant pickled uncertainty objects; use obs_val/obs_err columns."""

    def __setstate__(self, state):
        pass


class DataOnlyReader(pickle.Unpickler):
    def find_class(self, module, name):
        if module == 'numpy' and name in ('ndarray', 'dtype'):
            return getattr(np, name)
        if module in ('numpy.core.multiarray', 'numpy._core.multiarray') and name == '_reconstruct':
            return np.core.multiarray._reconstruct
        if module == 'uncertainties.core' and name in ('Variable', 'AffineScalarFunc', 'LinearCombination'):
            return OpaqueUncertainty
        raise ValueError(f'Unsupported external serialized type: {module}.{name}')


def columns(group):
    """Read this fixed-format dataframe without executing arbitrary pickle globals."""
    result = {}
    for i in range(int(group.attrs['nblocks'])):
        names = [v.decode() for v in group[f'block{i}_items'][:]]
        values = group[f'block{i}_values'][:]
        if values.dtype.kind == 'O':
            values = DataOnlyReader(io.BytesIO(bytes(values[0]))).load()
        assert values.ndim == 2 and values.shape[1] == len(names)
        result.update(zip(names, values.T))
    return result


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download():
    archive = ROOT / 'build/phasecraft_17799843.zip'
    if not archive.exists():
        archive.parent.mkdir(exist_ok=True)
        tmp = archive.with_suffix('.download')
        with urllib.request.urlopen(URL, timeout=60) as response, tmp.open('wb') as out:
            while block := response.read(1024*1024):
                out.write(block)
        tmp.replace(archive)
    assert hashlib.md5(archive.read_bytes()).hexdigest() == MD5, 'Archive checksum mismatch'
    with zipfile.ZipFile(archive) as z:
        for U in (0, 4):
            for name in ('shots.h5', 'exp_vals.h5'):
                member = f'fermionic_dynamics/U_{U}/{name}'
                path = ROOT / 'build' / member
                if not path.exists():
                    path.parent.mkdir(parents=True, exist_ok=True)
                    # Fixed member names, never paths supplied by the archive.
                    with z.open(member) as source, path.open('wb') as target:
                        while block := source.read(1024*1024):
                            target.write(block)
    return archive


def shot_observables(bits, labels):
    up = [labels.index(f'({x},{y},up)') for x in range(7) for y in range(4)]
    down = [labels.index(f'({x},{y},down)') for x in range(7) for y in range(4)]
    doublons = bits[:, up]*bits[:, down]
    accepted = (bits[:, up].sum(axis=1) == 14) & (bits[:, down].sum(axis=1) == 14)
    return doublons, accepted


def reduce_data():
    archive = download()
    sources = {str(archive.relative_to(ROOT)): sha(archive)}
    observable_rows, shot_rows, methods = [], [], {}
    for U in (0, 4):
        path = CACHE / f'U_{U}/shots.h5'
        sources[str(path.relative_to(ROOT))] = sha(path)
        with h5py.File(path) as h:
            labels = list(h['metadata'].attrs['qubit_order'])
            for time in sorted(h['data'], key=float):
                ds = h['data'][time]
                shape = ds.attrs['original_shape']
                bits = np.unpackbits(ds[:], axis=-1)[:, :shape[1]]
                assert tuple(bits.shape) == tuple(shape) and bits.shape[1] == 56
                dsamples, accepted = shot_observables(bits, labels)
                sample = dsamples.mean(axis=1)
                shot_rows.append(dict(U=U, tau=float(time), shots=len(bits), qubits=56,
                    accepted=int(accepted.sum()), number_acceptance=float(accepted.mean()),
                    mean_doublon=float(sample.mean()), variance_per_shot=float(sample.var(ddof=1)),
                    iid_standard_error=float(sample.std(ddof=1)/np.sqrt(len(bits))),
                    mean_local_bernoulli_variance=float(np.mean(dsamples.mean(axis=0)*(1-dsamples.mean(axis=0))))))
        path = CACHE / f'U_{U}/exp_vals.h5'
        sources[str(path.relative_to(ROOT))] = sha(path)
        with h5py.File(path) as h:
            methods[U] = {}
            for key in h:
                c = columns(h[key])
                methods[U][key] = sorted(set(map(str, c['method'])))
                # Keep the global mean and all single-site doublons, not an
                # optimistically selected time/site subset.
                selected = np.flatnonzero(np.isin(c['obs_type'], ['doublon_sum', 'Doublon']))
                for i in selected:
                    value, error = float(c['obs_val'][i]), float(c['obs_err'][i])
                    if not np.isfinite(value):
                        continue
                    summed = c['obs_type'][i] == 'doublon_sum'
                    coords = c['site_coords'][i]
                    bond = float(c['max_bond_dimension'][i])
                    method = str(c['method'][i])
                    # Some composite observables omit the numeric metadata,
                    # while their method label still gives the bond dimension.
                    match = re.fullmatch(r'TDVP chi=(\d+)', method)
                    if not np.isfinite(bond) and match:
                        bond = int(match[1])
                    observable_rows.append(dict(U=U, tau=float(c['ev_time'][i]),
                        method=method, group=key,
                        observable='mean_doublon' if summed else 'local_doublon',
                        site=None if summed else str(coords),
                        value=value/(28 if summed else 1),
                        reported_se=error/(28 if summed else 1) if np.isfinite(error) else None,
                        bond_dimension=int(bond) if np.isfinite(bond) else None))
    OUT.mkdir(parents=True, exist_ok=True)
    result = dict(source_url='https://zenodo.org/records/17799843', archive_md5=MD5,
        sources=sources, reducer_sha256=sha(Path(__file__)), methods=methods,
        shots=shot_rows, observables=observable_rows,
        scope='Published processed estimates, not a pipeline reimplementation. U=0 was used for training; U=4 references are finite-bond MPS. Shot arrays do not expose twirling identifiers; iid shot SE is diagnostic only.')
    (OUT/'phasecraft_reduced.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('Reduced', len(shot_rows), 'shot batches and', len(observable_rows), 'published observable rows')


if __name__ == '__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    reduce_data()
