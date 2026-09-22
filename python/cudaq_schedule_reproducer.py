"""Small fixed-angle circuits that exercise the rejected encoded timing path."""

from project_paths import ROOT, project_file
from contextlib import redirect_stdout
from dataclasses import asdict
import gzip
from io import StringIO
import json
from math import pi
from pathlib import Path

from cudaq_encoded_reference import run
from hubbard_circuits import Operation


def main():
    root=(ROOT / 'results/refinement_v3/schedule_reproducer');root.mkdir(exist_ok=True,parents=True)
    pattern=[Operation('rotate',(0,1),pi/8,'ZZ','fixed_T'),
             Operation('rotate',(1,),pi/8,'X','fixed_T'),
             Operation('rotate',(1,2),pi/8,'ZZ','fixed_T'),
             Operation('rotate',(2,),pi/8,'X','fixed_T')]
    rows=[]
    for repeat in (1,2,4,8):
        out=root/f'repeated_{repeat}'
        out.mkdir(exist_ok=True,parents=True)
        (out/'input_operations.json').write_text(json.dumps([asdict(op) for op in pattern*repeat],indent=2)+'\n')
        captured=StringIO()
        try:
            with redirect_stdout(captured):
                record=run(operations=pattern*repeat,output_directory=out)
            entries=json.loads((out/'schedule.json').read_text())
            top=[r for r in entries if r['parent'] is None]
            makespan=max(r['start_ns']+r['duration_ns'] for r in top)
            summary={'repetitions':repeat,'T_count':record['T_count'],'makespan_ns':makespan,
                     'loose_serial_upper_bound_ns':record['timing_sanity_bound_ns'],
                     'timing_sanity_passed':record['timing_sanity_passed'],
                     'status':record['status'],'largest_call_duration_ns':max(r['duration_ns'] for r in entries)}
        except Exception as exc:
            detail=str(exc)+'\nCAUSE:\n'+str(exc.__cause__)
            with gzip.open(out/'failure.txt.gz','wt') as f:f.write(detail)
            summary={'repetitions':repeat,'T_count':4*repeat,'timing_sanity_passed':False,
                     'status':'REJECTED by scheduler','exception':type(exc).__name__,
                     'message':str(exc).splitlines()[0],
                     'diagnostic':'nested schedule entry must lie within its parent interval'
                                  if 'nested schedule entry' in detail else str(exc).splitlines()[0]}
        summary['scope']='same surface macro and factory as v1; four fixed-angle operations repeated'
        rows.append(summary)
        # Keep compact exact IR; scheduling records stay readable JSON.
        raw=out/'physical_p3.mlir'
        with gzip.open(out/'physical_p3.mlir.gz','wb') as f:f.write(raw.read_bytes())
        raw.unlink()  # this script's just-created redundant raw copy
        print(json.dumps(summary),flush=True)
    (root/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':main()
