# CUDA-Q Logical: first compiled comparison

22 September 2026. Read the [updated resource comparison](results/cudaq_v1/compiled_comparison.png)
alongside the [gate-count and duration comparison](results/cudaq_v1/compilation_diagnostics.png).
The former combines actual circuit counts with the declared error and encoded
timing models. The latter shows what the physical compilation itself establishes.

The accepted precision is unchanged: **1% relative error at the provisional
signal scale D* = 0.1**, hence absolute error 0.001. The statistical allocation
is 0.0005, giving 7,377,759 ideal Hoeffding shots at 95% confidence. All points
here have tau = 1 and r = 4. Neither the signal scale nor the Trotter error has
been calibrated by a dynamics simulation.

## What now runs

An isolated `.venv-cudaq` contains `cudaq-logical==0.1.1.post1` and
`cuda-quantum-cu12==0.16.0.post1`; the complete package lock is
[requirements-cudaq.lock](requirements-cudaq.lock). These compiler runs use
the CPU. They do not require building the local CUDA-Q checkout or a GPU.

[hubbard_circuits.py](hubbard_circuits.py) implements the agreed open-boundary,
uncontrolled Hubbard workload, with spin-interleaved modes on a site snake.
Both termwise Strang and the Campbell-style plaquette construction are compiled
at L = 2, 4, 6, 8 through [cudaq_compile.py](cudaq_compile.py):

- **Physical:** explicit nearest-neighbor CNOT and one-qubit circuits, compiled
  to CUDA-Q Logical P0. Two implementations are retained: Jordan–Wigner parity
  ladders on a square-grid Hamiltonian path, and greedy adjacent fermionic-SWAP
  routing with a persistent mode permutation. A project-side ASAP scheduler
  counts idle qubit-time and verifies every CNOT is a grid edge. The latter
  route tracks the changed readout wires. Neither route is optimized.
- **STAR:** the analog Pauli-rotation program is compiled without synthesis,
  retaining signed angles, Pauli support, fixed-angle gates, and preparation
  and measurement operations. Injection, growth, retries, and protected
  Clifford timing still use the analytical reference, not compiled STAR gadgets.
- **FT:** CUDA-Q GridSynth actually synthesizes each signed angle class at
  the allocated operator precision. Whole-program Clifford+T synthesis at
  L = 2 and 4 reproduces the class-weighted T counts exactly. These full
  programs also pass CUDA-Q's Pauli-based-computation transformation. L = 6
  and 8 use the class-weighted synthesis counts, before inter-gate optimization.

This division follows CUDA-Q's [synthesis contract](https://nvidia.github.io/cuda-quantum/latest/preview/logical/use-cases/logical-synthesis.html):
Clifford+T synthesis and PBC normalization are device-independent P0 passes;
placement, encoded gadgets, and physical scheduling are additional stages.
The adapter explicitly converts our exp(-i theta P) convention to CUDA-Q's
exp(-i phi P/2), with phi = 2 theta.

## What changes the comparison

For L = 4, tau = 1, r = 4:

| Quantity | Termwise Strang | Plaquette |
|---|---:|---:|
| Variable Pauli rotations | 889 | 288 |
| Fixed pi/8 Pauli rotations | 0 | 576 |
| Dynamic total angle | 64 | 44 |
| Actual synthesized T count | 61,178 | 19,584 |
| CNOTs, JW parity route | 7,922 | 11,872 |
| CNOTs, greedy fSWAP route | 18,374 | 11,848 |
| Attempt duration, JW route | 101.32 us | 77.77 us |
| Attempt duration, fSWAP route | 51.16 us | 24.73 us |

Durations include 100 ns reset and 100 ns final readout. All native gates take
10 ns; disjoint gates may run concurrently. The square copy uses 36 physical
qubits for 32 modes. Gate counts include preparation. Final measurement reads
all modes; the local-doublon SPAM model charges its two readout bits.

**The synthesis fit was reasonably close; the physical gate proxy was not.**
For the plaquette point, GridSynth uses 66 T gates per variable rotation,
giving 19,584 total rather than the fit's 20,736. Its FT reference attempt
falls from 1.392 s to 1.315 s. Required area remains 13,791,250 qubits because
both cases round to 17 replicas of 811,250 qubits at d = 25.

The FT reference is actually factory-limited at this point: with 16 serial
two-level factories, supply costs about 11.16 clocks per T, exceeding the
9-clock consumption assumption. Changing consumption alone to one clock does
not improve this configuration. The CSV includes that sensitivity. Factory
parallelism and internal pipelining need to be varied together with consumption.

By contrast, the original physical proxy gave only 1,152 two-qubit gates for
this plaquette point. Explicit JW strings give 11,872. fSWAP routing gives
nearly the same CNOT count but substantially less idle time. For termwise
Strang this greedy fSWAP order increases the gate count while reducing depth.
Ordering, native hopping synthesis, and use of the full two-dimensional array
are therefore essential refinements before interpreting a physical crossover.

The new physical default is **p per 100 ns of idle time**; the old p-per-10-ns-layer
convention remains a named sensitivity. For one-qubit Pauli depolarization,
the composed idle probability is

    p_idle(dt) = (3/4) [1 - (1 - 4p/3)^(dt / 100 ns)].

At L = 4 plaquette, the log variance multiplier chi = log(Gamma^2) separates as:

| Exposure | JW, p/100 ns | fSWAP, p/100 ns |
|---|---:|---:|
| Two-qubit gates | 47.49 | 47.40 |
| One-qubit gates | 2.35 | 4.37 |
| Idle | 87.50 | 17.56 |

The old per-layer convention gives JW idle exposure 874.72. Even after changing
that convention, whole-circuit PEC on these unoptimized implementations gives
enormous sampling overhead. The plotted costs characterize those circuits and
that mitigation strategy; they do not establish a NISQ limit. Observable-aware
mitigation, locality, postselection, and better compilation remain open.
The physical panels truncate log10 Q at 100 and mark omitted points; the CSV
retains their full logarithmic values. The surface logical-error fit remains
provisional and has not been recalibrated to this continuous-time idle model.

## The first hybrid test

[compiled_model.py](compiled_model.py) adds distilled fixed pi/8 rotations plus
analog dynamic rotations, charging the factory, residual magic-state error,
analog injection area, logical faults, and replicas. This removes the fixed
gates' original-injection PEC contribution but keeps the dynamic angle and floor
exposure. It uses the same reference timing and error models as the other curves.

At the L = 4 plaquette point, all-angle STAR needs 166,698 qubits, whereas the
hybrid needs 66,150, about a **2.52-fold area improvement**. The hybrid uses one
single-level factory and d = 21; its attempt duration is 52.44 ms. At L = 2
the factory overhead outweighs the benefit. At L = 8 the hybrid switches to
two-level distillation under the specified error allocation. These are
conditional results of the current supply model, not optimized architecture
claims. For termwise Strang, there are no fixed angles, so the hybrid and STAR
curves coincide, as they should.

## C9: when logical faults are amplified

The review identifies a valid improvement under an additional condition.
If a branch-independent fault channel F is outside a complete PEC cancellation
block, or commutes through the inverted channel E, then E^-1 F E = F. If F is
a stochastic fault of probability q_L, a [0,1] observable has bias at most q_L.
Pauli fault channels and Pauli noise channels commute locally, so this is a
useful case for the proposed noise model.

Being a fixed linear map alone does not imply that equality. A small explicit
counterexample uses Bloch matrices E = diag(s,s,1) and
F = (1-q) I + q H, where H is the Hadamard conjugation channel. On the initial
state |0>, measuring (I+X)/2 after E^-1 F E gives bias q/(2s), exceeding q for
s < 1/2. The fault is stochastic and independent of the PEC branch. The issue
is its position and noncommutation, not the linearity of the estimator.

The reference retains Gamma q_L as a sufficient worst-case bound.
`STAR_commuting_faults` and `hybrid_commuting_faults` in the CSV remove this
factor under the stated condition. A compiled STAR gadget must establish where
residual faults occur relative to cancellation and whether retry branches share
the same residual channel before choosing the less conservative bound. This
condition has not yet been established for injection/growth/RUS schedules.

## Encoded backend status

[cudaq_surface_probe.py](cudaq_surface_probe.py) exercises the installed minimal
surface recipe. Preparation and readout compile, but its registered implementations
do not cover H, T, or arbitrary Rz. This does not mean CUDA-Q cannot implement
surface-code computation; this particular recipe is insufficient.

[cudaq_encoded_reference.py](cudaq_encoded_reference.py) adapts the supplied
Gidney–Ekera standalone example to the L = 2 plaquette circuit. It compiles
P0 -> P1 -> P2 -> P3 using surface patches, a typed T-state source, and explicitly
declared native encoded MPP/RPP operations. A small custom MPP adapter is needed
because the supplied Pinnacle WSC lowering only accepts one packed code block,
whereas the final PBC measurements span separate surface patches. The selected
d = 3 is an integration test, not a distance meeting our accuracy budget. Native
MPP/RPP, patch footprints, and factory service in this test are macro assumptions,
not a nearest-neighbor syndrome-extraction schedule.

**The Hubbard backend timing is rejected.** At r = 4, its 2,224-T program reports
2.19e18 ns, exceeding even a deliberately loose serial macro bound of 7.35e8 ns.
At r = 1, the estimator also rejects its own peak-concurrency accounting. In
contrast, a 20-rotation factory test returns 2.574 ms, consistent with the
128.53 us mean service interval. The root cause in the adapter/backend interaction
is not isolated. Raw physical IR, schedule entries, and explicit rejection
statuses are retained in `results/cudaq_v1/encoded_reference*/` and the positive
control in `factory_probe_20/`. None of these Hubbard backend timings enters a
resource curve. The validated P0 counts remain independent of this problem.

## Validation and next refinements

Forty tests pass in `qb_gsee`. New checks include full four-mode Fock-space
identities for the fermionic Hadamard and plaquette block, equality after NN
parity/fSWAP decomposition, common analytical count/angle identities, grid-edge
and idle accounting, and the logical-fault counterexample. Independently
reconstructing every emitted one-qubit Clifford+T word verifies its projective
operator error is within the allocated tolerance. These are circuit checks,
not exact Hubbard dynamics or Trotter-error calibration.

The next substantive replacements are:

1. Improve physical mode ordering, use both grid directions, and fuse hopping
   with fermionic swaps; compare the resulting gate/idle exposure with an
   observable-aware mitigation model.
2. Resolve the encoded scheduling discrepancy with a reduced reproducer before
   importing backend timing. Compile an explicit nearest-neighbor surface-code
   schedule and characterize its noise model.
3. Implement and schedule STAR injection/growth/retry gadgets, including parallel
   supply. Use their residual channels to settle the applicable C9 bound.
4. Optimize the hybrid and FT factory supply/consumption together.

Trotter calibration and qubitization retain their agreed later position.

The original analytical note, dedicated `references.bib`, and
`results/analytical_v1/` have been left unchanged. Per-case JSON records carry
the configuration and circuit-source hashes; the new manifest covers the
compiler package versions, source files, raw IR, diagnostics, and plots.
The living `WORKING_DECISIONS.md` has been updated; its hash in the older
manifest therefore predates this compilation pass. The original scientific
inputs and analytical code still match that manifest.
