module attributes {qlx.profiles = ["p3"]} {
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
    phys.release %busy {event_id = "release0"} : !phys.state<@q0>
    phys.release %c0 {event_id = "release1"} : !phys.state<@q1>
    phys.return
  }
}
