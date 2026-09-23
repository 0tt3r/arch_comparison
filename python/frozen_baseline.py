"""Check original source hashes without rewriting the frozen model record.

The archive is a copy of the five changed Python files from commit 1503a89.
This checks provenance only: callers separately replay the saved costs with
the current model. Non-source inputs must always match their live files.
"""
import hashlib
import zipfile

from project_paths import ROOT


def verify_sources(record):
    archived = 0
    for name, expected in record.items():
        current = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        if current == expected:
            continue
        assert name.startswith('python/') and name.endswith('.py'), name
        with zipfile.ZipFile(ROOT/'results/physics_sweep_v1/baseline_sources.zip') as archive:
            assert hashlib.sha256(archive.read(name)).hexdigest() == expected, name
        archived += 1
    if archived:
        print(f'Verified {archived} original Python sources from the frozen baseline archive')
