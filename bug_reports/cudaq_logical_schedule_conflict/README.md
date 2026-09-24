# cudaq-logical: `phys-schedule` rejects some repetition counts of a small program

## Summary

A three-qubit program of four π/8 Pauli rotations (Z₀Z₁, X₁, Z₁Z₂, X₂), repeated n times, is synthesized
to Clifford+T, placed on a distance-3 surface-code device with one T-state factory, and compiled through the
`qec()` and `physical()` pipelines. `cql.compiler.schedule(p3)` succeeds for n = 1, 2, 5, 6, 8 and fails for
n = 3 and 4 with

```
ScheduleConflict: native greedy-ASAP scheduling failed for phys.graph @prog_placed_qec_physical:
Pass pipeline failed: phys-schedule{graph=prog_placed_qec_physical result=prog_placed_qec_physical_schedule}
loc(unknown): 'phys.schedule' op nested schedule entry must lie within its parent interval
```

| n | T count | Result | Makespan (ns) |
|---:|---:|---|---:|
| 1 | 4 | scheduled | 515,557 |
| 2 | 8 | scheduled | 1,029,914 |
| 3 | 12 | **ScheduleConflict** | — |
| 4 | 16 | **ScheduleConflict** | — |
| 5 | 20 | scheduled | 2,572,984 |
| 6 | 24 | scheduled | 3,087,341 |
| 8 | 32 | scheduled | 4,116,055 |

Passing schedules are factory-bound, as expected: about 128.5 µs per T, which is the declared factory
output interval (178.5 × d cycles × 240 ns). The failure is not monotone in program size, which suggests an
interval-bookkeeping problem in the scheduler, perhaps where factory service overlaps nested patch calls,
rather than a resource limit.

## Reproduce

```bash
python repro.py            # n = 1 2 3 4 5 6 8, a few minutes on one CPU core
python repro.py 3          # a single failing case
```

`repro.py` is self-contained: it needs only `cudaq` with `cudaq.logical`. Outputs go to `output/`: per case
`result.json`, the lowered `physical_p3.mlir` (written before scheduling, so failing cases have it too), and
`traceback.txt` with the full exception and the rejected `phys.schedule` operation.

## Environment

See `environment.txt`: cudaq-logical 0.1.1.post1, CUDA-Q 0.16.0 (62fce8b), Python 3.12.0, RHEL 8 x86-64, CPU
only.

## How the device is built (possibly relevant)

- `cql.architectures.pinnacle.for_code(surface.definitions(3).code, cycle_rounds=3)`, with its
  `pauli_product_measurement` link root **replaced** by a custom `QECLowering` (`macro_mpp_lowering` in
  `repro.py`) that emits `cql.ops.mpp` over separate surface patches. We did this because the supplied
  lowering accepted only a single packed block, and our programs need Pauli products across patches.
- One logical factory producing `T_STATE` through a trivial `@cql.protocol`, bound to a `factory_lane`
  resource with `FactoryModel(startup_cycles = output_interval_cycles = 178.5 × 3)`.
- Physical patch resources with `NATIVE_PAULI_PRODUCT_ROTATION` and the MPP instrument; operating point
  `cycle_ns = 240`, `rpp_ns = 9·d·240`, `mpp_ns = d·240`, `condition_ns = 200`.
- `sys.setrecursionlimit(50000)`, because the ownership verifier recurses along SSA chains and larger
  programs exceed Python's default depth.

## Questions for the developers

1. Is the n = 3, 4 `ScheduleConflict` a known scheduler issue? If it comes from our device construction
   instead, what is wrong with it?
2. Is there, or will there be, a supported surface-code lattice-surgery lowering for multi-patch
   Pauli-product measurements and rotations? Our macro makes the schedule replay declared durations
   instead of compiling surgery. How should a multi-qubit Pauli rotation with T-state consumption be expressed
   so that the physical backend schedules real lattice surgery?
3. Can factory placement, T-state transport, and routing constraints enter the physical schedule, so that
   congestion and factory-to-data transport appear in the makespan?

## Context

This came from a resource comparison of fault-tolerant and partially fault-tolerant architectures for
Fermi–Hubbard dynamics. CUDA-Q Logical's compiled Clifford+T counts agree with our analytical ledger to
within 9% through an 8×8 lattice. The encoded physical timing is the part we could not use, because of this
failure and because the multi-patch surgery is a macro rather than a compiled gadget.
