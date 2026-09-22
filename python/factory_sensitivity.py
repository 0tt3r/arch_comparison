"""Separate physical error, injected-state error, precision, and replica rounding."""

from project_paths import ROOT, project_file
import csv
import json
from math import log, log10
from pathlib import Path

from compiled_model import compiled_star, encoded
from model import accuracy

OUT=ROOT/'results/nisq_v2'


def continuous_area(row,c):
    """Relax replica integrality for the chosen configuration, retaining >=1 copy."""
    _,shots=accuracy(c)
    return log10(row['qubits_per_copy'])+max(0,log10(shots)+row['log10_gamma2']
                       +log10(row['seconds_per_attempt']/c['week_seconds']))


def main():
    base=json.loads((project_file('assumptions.json')).read_text());rows=[]
    for L in (2,4,6,8):
        record=json.loads((ROOT/f'results/cudaq_v1/plaquette_L{L}_tau1_r4/result.json').read_text())
        for eps in (.001,.01):
            c={**base,'relative_error':eps/base['doublon_reference']}
            # Keeping the compiled word fixed isolates error allocation and
            # factory switching. It is conservative at the looser precision.
            star=compiled_star(record,c)
            rows.append({**star,'scenario':'STAR','absolute_error':eps,
                         'log10_Q_continuous_selected':continuous_area(star,c)})
            for factor in (1.,34/15,46/15,2/15):
                for arm in ('ft','hybrid'):
                    row=encoded(record,c,arm,injection_error=factor*c['p'])
                    rows.append({**row,'scenario':arm,'absolute_error':eps,'injection_factor':factor,
                        'injection_source':('Akahoshi Eq.6 prepared-state sensitivity; supply uncalibrated'
                                            if factor==2/15 else 'Akahoshi conventional injection or v1 baseline'),
                        'log10_Q_continuous_selected':continuous_area(row,c),
                        'gain_over_STAR':10**(star['log10_Q_min']-row['log10_Q_min']),
                        'compiled_precision':'same words compiled for epsilon=0.001',
                        'magic_budget':eps*c['magic_fraction']})
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'factory_sensitivity.json').write_text(json.dumps(rows,indent=2)+'\n')
    fields=sorted(set().union(*(row.keys() for row in rows)))
    with (OUT/'factory_sensitivity.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    for row in rows:
        if row['L']==4 and row.get('injection_factor') in (1.,46/15):
            print(row['scenario'],'eps',row['absolute_error'],'p_in/p',round(row['injection_factor'],4),
                  'levels',row['factory_level'],'Q',round(10**row['log10_Q_min']),
                  'continuous Q',round(10**row['log10_Q_continuous_selected']))


if __name__=='__main__':main()
