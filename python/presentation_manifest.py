"""Verify the presentation against frozen scientific inputs and current sources."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import sys

from plot_capacity import costs,maximum_tested_side,nisq_sampling_audit
from project_paths import ROOT

OUT=ROOT/'results/presentation_v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_content():
    saved=json.loads((OUT/'capacity_thresholds.json').read_text())
    assert saved==costs() and len(saved)==32
    assert json.loads((OUT/'nisq_sampling_audit.json').read_text())==nisq_sampling_audit()
    frontier=json.loads((OUT/'capacity_frontiers.json').read_text())
    for row in frontier:
        points=[r for r in saved if r['arm']==row['arm'] and r['absolute_error']==row['absolute_error']]
        expected=maximum_tested_side(points,row['Q'])
        assert row['L_max_tested']==expected or (row['L_max_tested'] is None and expected!=expected)
    manuscript=(ROOT/'tex/main.tex').read_text()
    assert '\\input{' not in manuscript,'Keep the complete main narrative in main.tex'
    assert '\\begin{thebibliography}' not in manuscript,'Use the dedicated bib file'
    labels=re.findall(r'\\label\{([^}]+)\}',manuscript)
    assert len(labels)==len(set(labels)),'Duplicate TeX labels'
    for label in re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',manuscript):assert label in labels,label
    entries=set(re.findall(r'@\w+\{([^,]+),',(ROOT/'tex/references.bib').read_text()))
    for group in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}',manuscript):
        for key in group.split(','):assert key.strip() in entries,key
    figures=re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',manuscript)
    for name in figures:assert (ROOT/'figures'/name).exists() or (ROOT/name).exists(),name
    for name in ('README.md','docs/REPRODUCING.md','tex/archive/README.md','results/README.md'):
        path=ROOT/name
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if '://' in target or target.startswith('#'):continue
            assert (path.parent/target.split('#')[0]).exists(),(name,target)
    assert (ROOT/'tex/main.pdf').read_bytes().startswith(b'%PDF-')
    discussion=json.loads((ROOT/'results/discussion_v3/validation.json').read_text())
    assert discussion['status'].startswith('PASS:')
    assert discussion['scenario_points']==len(json.loads((ROOT/'results/discussion_v3/costs.json').read_text()))
    delivery=json.loads((ROOT/'results/delivery_v1/validation.json').read_text())
    assert delivery['status'].startswith('PASS:')
    assert delivery['rows']==len(json.loads((ROOT/'results/delivery_v1/costs.json').read_text()))
    large=json.loads((ROOT/'results/dynamics_4x4/validation.json').read_text())
    assert large['status'].startswith('PASS:')
    followup=json.loads((ROOT/'results/dynamics_followup_v1/validation.json').read_text())
    assert followup['status'].startswith('PASS:')
    physics=json.loads((ROOT/'results/physics_sweep_v1/validation.json').read_text())
    assert physics['status'].startswith('PASS:')
    assert physics['sweep_rows']==len(json.loads((ROOT/'results/physics_sweep_v1/costs.json').read_text()))
    implementation=json.loads((ROOT/'results/implementation_v1/validation.json').read_text())
    assert implementation['status'].startswith('PASS:')
    assert implementation['growth_shots']==10000000
    frontier_check=json.loads((ROOT/'results/frontier_v1/validation.json').read_text())
    assert frontier_check['status'].startswith('PASS:') and frontier_check['growth_shots']==24000000
    return {'frontier_algorithm_costs':frontier_check['algorithm_costs'],
            'frontier_growth_shots':frontier_check['growth_shots'],'threshold_rows':len(saved),'frontier_points':len(frontier),
            'discussion_scenario_points':discussion['scenario_points'],
            'discussion_calibration_cases':discussion['calibration_cases'],
            'delivery_scenario_points':delivery['rows'],'large_patch_measurements':large['measured_rows'],
            'direct_allocation_checks':followup['allocation_checks'],'patch_sensitivity_rows':followup['cost_rows'],
            'physical_error_sweep_rows':physics['sweep_rows'],'algorithm_alternative_rows':physics['algorithm_rows'],
            'held_out_estimator_checks':implementation['held_out_estimates'],
            'decoded_growth_shots':implementation['growth_shots'],'explicit_ft_stages':implementation['ft_stages'],
            'manuscript_figures':len(figures),'status':'PASS: discrete inversion, manuscript references, and artifact paths'}


def inventory():
    sources=[ROOT/'Makefile',ROOT/'pytest.ini',ROOT/'README.md']
    for folder,pattern in [('python','*.py'),('tests','*.py'),('configs','*'),('tex','*.tex'),
                           ('tex','*.bib'),('docs','*.md'),('cpp','*.cpp'),('pbs','*.pbs')]:
        sources.extend(p for p in (ROOT/folder).glob(pattern) if p.is_file())
    sources.extend([ROOT/'results/README.md',ROOT/'tex/archive/README.md'])
    inputs=[ROOT/'results/refinement_v3'/name for name in (
        'manifest.json','cultivation.json','star_parallel.json','small_dynamics.json',
        'noisy_small.json','accuracy_calibration.pdf','postselection_bias.pdf')]
    inputs.extend(ROOT/f'results/refinement_v3/parallel_L{L}_observable/result.json' for L in (2,4,6,8))
    artifacts=[ROOT/'tex/main.pdf']+sorted((ROOT/'figures').glob('*'))
    for folder in sorted((ROOT/'results').glob('discussion_v*'))+[ROOT/'results/delivery_v1',ROOT/'results/dynamics_4x4',ROOT/'results/dynamics_followup_v1',ROOT/'results/physics_sweep_v1',ROOT/'results/workshop_reconciliation_v1',ROOT/'results/implementation_v1',ROOT/'results/frontier_v1']:
        artifacts.extend(p for p in sorted(folder.rglob('*')) if p.is_file())
    artifacts.extend(p for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='manifest.json')
    return {key:{str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))}
            for key,paths in [('sources',sources),('scientific_inputs',inputs),('artifacts',artifacts)]}


def main(write=False):
    report=check_content();path=OUT/'manifest.json'
    current=inventory()
    if write:
        record={**current,'checks':report,'python':sys.version,
                'packages':{n:importlib.metadata.version(n) for n in ('numpy','scipy','matplotlib','pytest')},
                'scope':'Current discussion calculations, presentation, and repository organization; prior versioned scientific artifacts remain unchanged'}
        path.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    manifest=json.loads(path.read_text())
    for category,values in current.items():
        assert values==manifest[category],f'{category} changed; regenerate the paper and presentation manifest'
    assert report==manifest['checks']
    print(json.dumps(report,indent=2))
    print('Verified',sum(len(v) for v in current.values()),'source, input, and artifact hashes')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args();main(args.write)
