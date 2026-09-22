"""Exercise file lookup and compressed-artifact restoration after relocation."""
import gzip
import hashlib
import json

from plot_capacity import costs


def test_current_cost_inputs_do_not_depend_on_working_directory(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    rows=costs()
    assert len(rows)==32
    assert {r['arm'] for r in rows}=={'nisq','star','ft','hybrid'}


def test_restore_uses_artifact_directory_and_preserves_existing_bytes(tmp_path,monkeypatch):
    import package_results
    data=b'module { // bounded restoration fixture\n}\n'
    packed=gzip.compress(data,mtime=0)
    (tmp_path/'case.mlir.gz').write_bytes(packed)
    digest=lambda b:hashlib.sha256(b).hexdigest()
    record={'files':{'case.mlir':{'stored_as':'case.mlir.gz',
            'sha256':digest(data),'stored_sha256':digest(packed)}}}
    (tmp_path/'storage.json').write_text(json.dumps(record))
    monkeypatch.setattr(package_results,'ROOT',tmp_path)
    package_results.restore()
    package_results.restore()
    assert (tmp_path/'case.mlir').read_bytes()==data
