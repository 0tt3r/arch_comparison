"""Two-resource scheduler reproducer with no QEC, factory, or Hubbard circuit."""
from dataclasses import asdict
import json
from pathlib import Path

import cudaq.logical as cql
import cudaq.mlir.ir as ir
from cudaq.logical.programs import DefinitionHandle


def graph(reuses):
    text='''module attributes {qlx.profiles = ["p3"]} {
  phys.machine @arch {
    phys.resource_class @q {count = 2 : i64, kind = "qubit", native_actions = []}
  }
  phys.resource @q0 {index = 0 : i64, kind = "qubit", resource_class = @q}
  phys.resource @q1 {index = 1 : i64, kind = "qubit", resource_class = @q}
  phys.graph @g on @arch : () -> () {
    %a = phys.acquire [@q0] {event_id = "acquire0"} : !phys.state<@q0>
    %b = phys.acquire [@q1] {event_id = "acquire1"} : !phys.state<@q1>
    %busy = phys.delay %a {duration_ns = 1.000000e+03 : f64, event_id = "unrelated_delay"}
      : (!phys.state<@q0>) -> !phys.state<@q0>
    %c0 = "phys.call"(%b) <{callee = @work, event_id = "call0", instance = "root.work.call0"}> ({
    ^bb0(%current: !phys.state<@q1>):
      %next = phys.delay %current {duration_ns = 3.000000e+00 : f64, event_id = "short_delay"}
        : (!phys.state<@q1>) -> !phys.state<@q1>
      phys.yield %next : !phys.state<@q1>
    }) : (!phys.state<@q1>) -> !phys.state<@q1>
'''
    for i in range(1,reuses+1):
        text+=f'''    %c{i} = phys.call_template %c{i-1} {{callee = @work, event_id = "call{i}",
      instance = "root.work.call{i}", template_event = "call0"}}
      : (!phys.state<@q1>) -> !phys.state<@q1>
'''
    return text+f'''    phys.release %busy {{event_id = "release0"}} : !phys.state<@q0>
    phys.release %c{reuses} {{event_id = "release1"}} : !phys.state<@q1>
    phys.return
  }}
}}
'''


def run():
    out=Path('results/refinement_v3/minimal_timing');out.mkdir(exist_ok=True,parents=True);rows=[]
    for count in (0,1,2,4):
        text=graph(count);(out/f'reuses{count}.mlir').write_text(text)
        context=ir.Context();module=ir.Module.parse(text,context);assert module.operation.verify()
        build=cql.compiler.Build(context=context,module=module,root=DefinitionHandle('g','phys.graph','p3'),profile='p3',pipeline=None)
        scheduled=cql.compiler.schedule(build)
        record={'reuses':count,'expected_makespan_ns':max(1000.,3.*(count+1)),
                'reported_makespan_ns':scheduled.makespan_ns,
                'expected_call_duration_ns':3.,
                'reported_call_durations_ns':[e.duration_ns for e in scheduled.entries if e.kind in ('call','call_template')]}
        record['passes']=record['reported_makespan_ns']==record['expected_makespan_ns'] and all(t==3 for t in record['reported_call_durations_ns'])
        rows.append(record)
        (out/f'reuses{count}_schedule.json').write_text(json.dumps([asdict(e) for e in scheduled.entries],indent=2)+'\n')
        print(json.dumps(record),flush=True)
    (out/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':run()
