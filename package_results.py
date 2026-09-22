"""Pack or restore full v1 MLIR dumps while retaining original artifact hashes.

Small angle-synthesis fixtures remain plain text for direct test access.
Packing retains the local originals; restoration never replaces different data.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent / 'results' / 'cudaq_v1'


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def pack():
    manifest_bytes = (ROOT / 'manifest.json').read_bytes()
    manifest = json.loads(manifest_bytes)
    # The compiler inventories its output directory recursively. If it ran
    # after an earlier packaging pass, discard those stale packaging entries;
    # their current hashes belong in storage.json, not the raw-artifact manifest.
    expected = {name: value for name, value in manifest['artifacts'].items()
                if name != 'storage.json' and not name.endswith('.mlir.gz')}
    if expected != manifest['artifacts']:
        manifest['artifacts'] = expected
        manifest_bytes = (json.dumps(manifest, indent=2) + '\n').encode()
        (ROOT / 'manifest.json').write_bytes(manifest_bytes)
    files = {}
    for path in sorted(ROOT.rglob('*.mlir')):
        if path.name.startswith('angle_'):
            continue
        name = str(path.relative_to(ROOT))
        data = path.read_bytes()
        digest = sha256(data)
        if expected.get(name) != digest:
            raise ValueError(f'Original artifact does not match its manifest: {name}')
        packed = gzip.compress(data, compresslevel=6, mtime=0)
        destination = path.with_suffix('.mlir.gz')
        destination.write_bytes(packed)
        files[name] = {'stored_as': str(destination.relative_to(ROOT)),
                       'sha256': digest, 'stored_sha256': sha256(packed),
                       'bytes': len(data), 'stored_bytes': len(packed)}
    if not files:
        raise ValueError('No full MLIR dumps found; use --restore on a fresh checkout')
    record = {'original_manifest_sha256': sha256(manifest_bytes), 'files': files}
    (ROOT / 'storage.json').write_text(json.dumps(record, indent=2) + '\n')
    print(f"Packed {len(files)} dumps: "
          f"{sum(r['bytes'] for r in files.values()) / 1e6:.1f} MB -> "
          f"{sum(r['stored_bytes'] for r in files.values()) / 1e6:.1f} MB")


def restore():
    record = json.loads((ROOT / 'storage.json').read_text())
    for name, entry in record['files'].items():
        packed = (ROOT / entry['stored_as']).read_bytes()
        if sha256(packed) != entry['stored_sha256']:
            raise ValueError(f'Compressed artifact checksum failed: {name}')
        data = gzip.decompress(packed)
        if sha256(data) != entry['sha256']:
            raise ValueError(f'Original artifact checksum failed: {name}')
        path = ROOT / name
        if path.exists():
            if sha256(path.read_bytes()) != entry['sha256']:
                raise ValueError(f'Refusing to replace different local data: {name}')
        else:
            path.write_bytes(data)
    print(f"Verified/restored {len(record['files'])} original dumps")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore', action='store_true')
    args = parser.parse_args()
    restore() if args.restore else pack()
