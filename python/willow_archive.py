"""Fetch selected members of the pinned ZIP via checked HTTP byte ranges.

Only about 1 GB of the 8.8 GB archive is needed for the declared comparison.
Each extracted member is checked against its ZIP size and CRC and hashed.
The full published archive MD5 is recorded, not claimed verified by this path.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import struct
import time
import urllib.error
import urllib.request
import zipfile
import zlib

from project_paths import ROOT

OUT = ROOT/'results/willow_v1'
CACHE = ROOT/'build/willow'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


class RangeReader(io.RawIOBase):
    def __init__(self, url, size):
        self.url, self.size, self.pos = url, size, 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = offset + (0 if whence==0 else self.pos if whence==1 else self.size)
        if self.pos < 0:
            raise ValueError('Negative archive position')
        return self.pos

    def read(self, count=-1):
        count = min(self.size-self.pos, count if count>=0 else self.size-self.pos)
        if count <= 0:
            return b''
        expected = f'bytes {self.pos}-{self.pos+count-1}/{self.size}'
        for attempt in range(3):
            try:
                request = urllib.request.Request(self.url, headers={
                    'Range':f'bytes={self.pos}-{self.pos+count-1}'})
                with urllib.request.urlopen(request, timeout=45) as response:
                    if response.status != 206 or response.headers.get('Content-Range') != expected:
                        raise ValueError('Server did not return the requested archive range')
                    data = response.read(count+1)
                if len(data) != count:
                    raise ValueError('Truncated archive range')
                self.pos += count
                return data
            except (urllib.error.URLError, TimeoutError):
                if attempt==2:
                    raise
                time.sleep(2)


def archive_index(reader):
    with zipfile.ZipFile(reader) as archive:
        return [dict(name=v.filename, size=v.file_size, compressed=v.compress_size,
                     offset=v.header_offset, crc=v.CRC, compression=v.compress_type)
                for v in archive.infolist()]


def selected_members(cfg, index):
    base = 'programmable_simulation_fermi_hubbard'
    names = {f'{base}/initial_states.json'}
    for L in cfg['sizes']:
        for flux in cfg['fluxes']:
            for U in cfg['interactions']:
                for r in cfg['steps']:
                    prefix = f'{base}/{cfg["initial_state"]}/{L}x{L}/mag_field_{flux}/U_{U}/trotter_layers_{r}'
                    names.update(f'{prefix}/{file}' for file in ('shots.h5','exp_vals.h5'))
    selected = [v for v in index if v['name'] in names]
    if len(selected) != len(names):
        raise ValueError('Archive does not contain the complete declared workload')
    return selected


def unpack_member(reader, entry):
    reader.seek(entry['offset'])
    header = reader.read(30)
    if header[:4] != b'PK\x03\x04':
        raise ValueError('Invalid member header')
    n, x = struct.unpack('<HH', header[26:30])
    name = reader.read(n).decode('utf-8')
    if name != entry['name']:
        raise ValueError('Member name differs from index')
    reader.seek(x, 1)
    compressed = reader.read(entry['compressed'])
    if entry['compression']==8:
        data = zlib.decompress(compressed, -15)
    elif entry['compression']==0:
        data = compressed
    else:
        raise ValueError('Unsupported ZIP compression')
    if len(data) != entry['size'] or zlib.crc32(data) != entry['crc']:
        raise ValueError('Member size or CRC mismatch')
    return data


def fetch():
    cfg = read(ROOT/'configs/willow_benchmark.json')
    reader = RangeReader(cfg['url'], cfg['archive_bytes'])
    index_path = OUT/'archive_index.json'
    if index_path.exists():
        index = read(index_path)
    else:
        index = archive_index(reader)
        save(index_path, index)
    selected = selected_members(cfg, index)
    records = []
    for i, entry in enumerate(selected):
        path = CACHE/entry['name']
        if path.exists():
            data = path.read_bytes()
            valid = len(data)==entry['size'] and zlib.crc32(data)==entry['crc']
        else:
            valid = False
        if not valid:
            data = unpack_member(reader, entry)
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix('.download')
            temporary.write_bytes(data)
            temporary.replace(path)
        records.append(dict(entry, sha256=hashlib.sha256(data).hexdigest()))
        print(f'{i+1}/{len(selected)} {entry["name"]}', flush=True)
    save(OUT/'downloads.json', dict(url=cfg['url'], archive_bytes=cfg['archive_bytes'],
        archive_md5_verified=False, published_archive_md5=cfg['published_archive_md5'],
        selected_compressed_bytes=sum(v['compressed'] for v in selected),
        index_sha256=sha(index_path), config_sha256=sha(ROOT/'configs/willow_benchmark.json'),
        downloader_sha256=sha(Path(__file__)), files=records,
        verification='Selected member sizes and CRCs from pinned HTTPS archive; SHA-256 of every extracted file. Full archive MD5 not checked.'))


if __name__ == '__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    fetch()
