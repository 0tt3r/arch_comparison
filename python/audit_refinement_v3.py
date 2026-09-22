"""Replay compact v3 evidence and freeze or verify its provenance manifest.

Run in the science environment. Creating the manifest queries the local CUDA-Q
environment for package versions but does not compile or execute quantum jobs.
"""

from project_paths import ROOT, project_file
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.metadata
import json
from math import isclose
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

from hubbard_circuits import Operation
from model import workload
from nisq_observables import parity_acceptance
from run_nisq_refinement import describe

OUT=ROOT/'results/refinement_v3'
SOURCE_NAMES=(
    'parallel_nisq.py','run_parallel_refinement.py','small_dynamics.py',
    'noisy_small_observable.py','cultivation_model.py','star_parallel.py',
    'cudaq_star_gadget.py','cudaq_schedule_reproducer.py','cudaq_minimal_timing.py',
    'cudaq_encoded_reference.py','cudaq_compile.py','compiled_model.py','model.py',
    'hubbard_circuits.py','nisq_compilation.py','nisq_observables.py','run_nisq_refinement.py',
    'plot_refinement_v3.py','audit_refinement_v3.py','test_refinement_v3.py',
    'test_model.py','test_hubbard_circuits.py','test_compiled_model.py','test_nisq_refinement.py',
    'assumptions.json','cultivation_assumptions.json','requirements-cudaq.lock',
    'refinement_v3.tex','refinement_v3.pdf','refinement_v3_references.bib',
    'references.bib','refinement_references.bib','REFINEMENT_V3.md','README.md','WORKING_DECISIONS.md')


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(actual,expected):
    if isinstance(expected,dict):
        assert actual.keys()==expected.keys()
        for key in expected:close(actual[key],expected[key])
    elif isinstance(expected,(list,tuple)):
        assert len(actual)==len(expected)
        for a,b in zip(actual,expected):close(a,b)
    elif isinstance(expected,float):
        assert isclose(actual,expected,rel_tol=2e-12,abs_tol=2e-13),(actual,expected)
    else:assert actual==expected,(actual,expected)


def replay():
    paths=sorted(OUT.glob('parallel*/result.json'));assert len(paths)==9
    sources=OUT/'sources'
    for path in paths:
        record=read(path);c=record['configuration']
        with gzip.open(path.parent/'native_gates.json.gz','rt') as f:
            gates=[Operation(**{**g,'qubits':tuple(g['qubits'])}) for g in json.load(f)]
        mapping={int(k):v for k,v in record['mode_to_wire'].items()}
        w=workload(record['L'],record['tau'],record['steps'],record['algorithm'],c['u_over_j'])
        repeated,events=describe(gates,mapping,record['observable_modes'],c,w)
        for key in ('schedule','observable_exposure','full_circuit_chi','costs_by_idle_ratio'):
            close(repeated[key],record[key])
        close(parity_acceptance(events,repeated['schedule'],mapping,c),record['postselection'])
        with gzip.open(path.parent/'native_p0.mlir.gz','rt') as f:mlir=f.read()
        actual=Counter('qlx_standard_'+name for name in re.findall(r'qlx\.apply #qlx\.action<([^>]+)>',mlir))
        expected=Counter('qlx_standard_'+('pauli_rotation' if g.name in ('rz','rx','rzz') else g.name) for g in gates)
        assert actual==expected==record['cudaq_logical']['actions']
        assert mlir.count('qlx.prepare "zero"')==record['schedule']['width']
        for name,digest in record['sources'].items():
            assert sha(sources/(digest+'__'+name))==digest,(path,name)

    dynamics=read(OUT/'small_dynamics.json')
    assert len(dynamics['trotter'])==112 and len(dynamics['curves'])==164
    assert len(dynamics['truncation'])==4
    for row in dynamics['trotter']:
        close(row['absolute_trotter_error'],abs(row['doublon_trotter']-row['doublon_exact']))
        assert row['norm_error']<1e-10
    for row in dynamics['curves']:
        assert row['norm_error']<1e-10 and -1e-12<=row['doublon']<=1+1e-12
    for row in dynamics['truncation']:
        close(row['finite_reference_difference'],abs(row['reference_D']-row['truncated_D']))
    noisy=read(OUT/'noisy_small.json');assert len(noisy)==3
    for row in noisy:
        suffix='' if row['r']==4 else f"_r{row['r']}"
        native=read(OUT/f'parallel_L2_observable{suffix}'/'result.json')
        close(row['spin_parity']['acceptance'],native['postselection']['acceptance'])
        assert row['trace_error']<2e-11 and row['hermiticity_error']<1e-12
        assert row['min_eigenvalue']>=-1e-12
        for selection in ('unselected','spin_parity','spin_particle_number'):
            values=row[selection]
            close(values['bias_from_noiseless_circuit'],values['doublon']-row['ideal_circuit_D'])
            close(values['bias_from_exact_dynamics'],values['doublon']-row['exact_dynamics_D'])

    from compiled_model import compiled_star,encoded
    from cultivation_model import cultivated_encoded
    from star_parallel import parallel_star
    c=read(project_file('assumptions.json'));source=read(project_file('cultivation_assumptions.json'))
    cultivated=read(OUT/'cultivation.json');assert len(cultivated)==88
    for row in cultivated:
        record=read(ROOT/f"results/cudaq_v1/plaquette_L{row['L']}_tau1_r4/result.json")
        scenario={**c,'relative_error':row['absolute_error']/c['doublon_reference']}
        if row['source']=='STAR reference':expected=compiled_star(record,scenario)
        elif row['source']=='conventional distillation':expected=encoded(record,scenario,row['arm'],injection_error=46*c['p']/15)
        else:
            expected=cultivated_encoded(record,scenario,row['arm'],source,row['handoff_error'],
                                       fixed_slot=row['source']=='cultivation fixed slot')
        close({key:row[key] for key in expected},expected)
    parallel=read(OUT/'star_parallel.json');assert len(parallel)==64
    for row in parallel:
        scenario={**c,'relative_error':row['absolute_error']/c['doublon_reference']}
        expected=parallel_star(row['L'],scenario,row['lanes'],row['logical_bias'])
        close({key:row[key] for key in expected},expected)
    gadgets=read(OUT/'star_gadget/results.json')
    assert [g['cap'] for g in gadgets]==[1,2,4,8,24]
    for row in gadgets:
        close(row['ideal_exhaustion_probability'],2.**(-row['cap']))
        with gzip.open(OUT/f"star_gadget/angle_doubling_cap{row['cap']}_p0.mlir.gz",'rt') as f:
            text=f.read()
        assert 'qlx.program' in text and len(text)>1000
    backend=read(OUT/'schedule_reproducer/summary.json')
    assert [r['T_count'] for r in backend]==[4,8,16,32]
    rejected=[r for r in backend if not r['timing_sanity_passed']]
    assert len(rejected)==1 and rejected[0]['T_count']==16
    with gzip.open(OUT/'schedule_reproducer/repeated_4/failure.txt.gz','rt') as f:
        assert 'nested schedule entry must lie within its parent interval' in f.read()
    controls=read(OUT/'minimal_timing/summary.json');assert len(controls)==4
    assert all(r['passes'] for r in controls)
    suites=ET.parse(OUT/'tests.xml').getroot().findall('testsuite')
    tests=sum(int(s.attrib['tests']) for s in suites)
    assert tests>=73 and all(int(s.attrib['failures'])+int(s.attrib['errors'])==0 for s in suites)
    return {'native_compiled_cases':len(paths),'trotter_comparisons':112,'exact_curve_points':164,
            'finite_spatial_comparisons':4,'noisy_density_cases':3,'cultivation_scenarios':88,
            'star_parallel_scenarios':64,'logical_retry_caps':5,'scheduler_diagnostics':4,
            'scheduler_positive_controls':4,'tests_passed':tests,
            'status':'PASS: saved-model consistency and provenance; physical encoded schedules unvalidated'}


def write_manifest():
    snapshots=OUT/'sources';snapshots.mkdir(exist_ok=True)
    for name in SOURCE_NAMES:
        path=project_file(name)
        (snapshots/(sha(path)+'__'+name)).write_bytes(path.read_bytes())
    audit=replay()
    (OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    packages={name:importlib.metadata.version(name) for name in ('numpy','scipy','matplotlib','pytest')}
    query="import json,sys,importlib.metadata as m; print(json.dumps({'python':sys.version,'packages':{n:m.version(n) for n in ['cudaq-logical','cuda-quantum-cu12','numpy','scipy']}}))"
    compiler=json.loads(subprocess.check_output([str(ROOT/'.venv-cudaq/bin/python'),'-c',query],text=True))
    inputs=[ROOT/f'results/cudaq_v1/plaquette_L{L}_tau1_r4/result.json' for L in (2,4,6,8)]
    inputs += [ROOT/f'results/nisq_v2/plaquette_L{L}_r4_fused_cone/result.json' for L in (2,4,6,8)]
    result={'created_utc':datetime.now(timezone.utc).isoformat(),
            'base_git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'science_environment':{'python':sys.version,'packages':packages},
            'compiler_environment':compiler,
            'sources':{name:sha(project_file(name)) for name in SOURCE_NAMES},
            'prior_result_inputs':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
            'artifacts':{str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob('*'))
                         if p.is_file() and p.name!='manifest.json'},
            'scope':audit['status']}
    (OUT/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')


def verify_manifest():
    manifest=read(OUT/'manifest.json')
    historical_sources=0
    for key,base in [('sources',ROOT),('prior_result_inputs',ROOT),('artifacts',OUT)]:
        for name,digest in manifest[key].items():
            path=project_file(name) if key=='sources' else base/name
            if key=='sources' and (not path.exists() or sha(path)!=digest):
                # The v3 evidence remains frozen across repository reorganization.
                # Verify its original source, then replay with the current code.
                path=OUT/'sources'/(digest+'__'+Path(name).name)
                historical_sources+=1
            assert sha(path)==digest,(key,name)
    current={str(p.relative_to(OUT)) for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'}
    assert current==manifest['artifacts'].keys(),'Unmanifested or missing artifacts'
    close(replay(),read(OUT/'audit.json'))
    print(json.dumps(read(OUT/'audit.json'),indent=2))
    print('Verified',len(manifest['artifacts']),'artifact hashes and',len(manifest['sources']),
          'source hashes;',historical_sources,'use the archived original; current code replay passed')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-manifest',action='store_true')
    args=parser.parse_args()
    if args.write_manifest:write_manifest()
    verify_manifest()
