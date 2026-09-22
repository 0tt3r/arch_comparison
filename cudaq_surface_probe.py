"""Exercise the installed surface recipe without treating it as our NN target."""
import json
from pathlib import Path

import cudaq.logical as cql
from cudaq.logical.targets import recipes

from cudaq_compile import program_from_operations
from hubbard_circuits import Operation


def run():
    architecture = recipes.surface_architecture(3)
    builder = cql.devices.DeviceBuilder("ReferenceSurfaceProbe")
    region = builder.logical.add_compute(capacity=1)
    encoded = builder.qec.bind(region, architecture=architecture)
    physical = builder.physical.add_qubits(
        architecture.code.block.size,
        native_actions=cql.architecture.physical_actions.clifford_set(),
        native_instruments=(cql.architecture.physical_instruments.MZ,
                            cql.architecture.physical_instruments.MPP))
    builder.physical.bind(encoded,to=physical)
    builder.physical.set_operating_point(timing={"cycle_ns":240.})
    device = builder.build()
    records = []
    for name, ops in [("zero_memory", []), ("clifford_h", [Operation("h",(0,))]),
                      ("T_gate", [Operation("t",(0,))]),
                      ("analog_rz", [Operation("rz",(0,),.1234)])]:
        try:
            logical = cql.compile(program_from_operations(ops,1,name))
            placed = cql.compiler.place(logical,device=device)
            encoded = cql.compile(placed,pipeline=cql.compiler.pipelines.qec(),device=device)
            physical = cql.compile(encoded,pipeline=cql.compiler.pipelines.physical(),device=device)
            scheduled = cql.compiler.schedule(physical)
            resources = cql.estimate(scheduled,tier=cql.estimate.Tier.SCHEDULE,
                                     p_phys=.001,failure_budget=.1,cycle_time=240e-9)
            record = {"case":name,"status":"compiled","estimate":resources.to_dict()}
        except Exception as exc:
            record = {"case":name,"status":"unsupported_by_reference_recipe",
                      "exception":type(exc).__name__,"message":str(exc)}
        records.append(record)
        print(json.dumps(record),flush=True)
    out=Path("results/cudaq_v1/surface_probe.json")
    out.write_text(json.dumps(records,indent=2)+"\n")


if __name__ == "__main__":
    run()
