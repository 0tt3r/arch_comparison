"""Compile an explicit angle-doubling STAR injection protocol at logical P0.

The output includes an exhaustion flag; exhausted attempts must be discarded.
This is a logical gadget and does not replace surface-code injection/growth.
"""

from project_paths import ROOT, project_file
import gzip
import json
from pathlib import Path

import cudaq.logical as cql


def inject(data,theta,remaining):
    ancilla=cql.allocate(1,state=cql.types.zero)[0]
    ancilla=cql.h(ancilla)
    ancilla=cql.rz(ancilla,2*theta)
    data,ancilla=cql.cx(data,ancilla)
    failed=cql.measure_z(ancilla)
    if remaining>1:
        data,failed=cql.ops.cond(failed,
            then=lambda live,flag:inject(live,2*theta,remaining-1),
            else_=lambda live,flag:(live,flag),carries=(data,failed))
    return data,failed


def program(theta,cap):
    @cql.program
    def rotation() -> tuple[bool,bool]:
        data=cql.allocate(1,state=cql.types.zero)[0]
        data=cql.h(data)
        data,exhausted=inject(data,theta,cap)
        data=cql.h(data)
        return cql.measure_z(data),exhausted
    return rotation


def run():
    out=(ROOT / 'results/refinement_v3/star_gadget');out.mkdir(exist_ok=True,parents=True)
    rows=[]
    for cap in (1,2,4,8,24):
        build=cql.compile(program(.137,cap))
        counts=cql.estimate(build,tier=cql.estimate.Tier.LOGICAL).to_dict()
        with gzip.open(out/f'angle_doubling_cap{cap}_p0.mlir.gz','wt') as f:f.write(build.to_mlir())
        rows.append({'cap':cap,'theta':.137,'ideal_exhaustion_probability':2.**(-cap),
                     'worst_case_prepared_states':cap,'expected_attempts':2*(1-2.**(-cap)),
                     'cudaq_logical':counts,'scope':'logical injection and conditional angle doubling; physical preparation and QEC unspecified'})
        print('Compiled STAR RUS cap',cap,flush=True)
    (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':run()
