# Plan: NISQ vs partial-FT (STAR) vs full-FT cost of Fermi-Hubbard dynamics on a 2D sea of qubits

Written 2026-09-22 by Claude (Fable 5.1) after reading the three refs and surveying
`~/soft/cuda-quantum` and `~/soft/pyLIQTR`.

## 0. Verdict

Doable in a hackathon. The core deliverable is a ~300-line analytic cost model plus
a phase diagram; the tools (CUDA-Q Logical, pyLIQTR, Azure QRE) calibrate its
constants rather than replace it. Main risk is building CUDA-Q Logical; Azure QRE
(already installed) takes it off the critical path.

## 1. What the references actually say (and how it changes the plan)

**`refs/2603.13093v1.pdf` (Chung et al., Mar 2026)** uses **TopQAD** (1QBit), not
"topQAC", and only for the full-FT side. It compares STAR vs full FT on p-benzyne and
2D Hubbard (L = 4..10, U/t = 4, 8; ε = 0.005 L² and ε = 0.01) via QCELS phase
estimation. **There is no NISQ / error-mitigated track**, so our NISQ arm is new.
Findings to reuse:
- PEC dominates STAR runtime. Hubbard at ε = 0.01 on their "target" hardware
  (2q error 5e-4, Λ = 9.3) is feasible only to L = 4 (5.2 days, 7.1 h without PEC);
  L = 6 is 16 years.
- "Goldilocks zone" ≈ 1e5–1e6 small-angle rotations where STAR beats distillation
  on factory footprint; below ~1e5 magic-state cultivation wins; above ~1e6 PEC is
  prohibitive.
- Full-FT Hubbard (Trotter+HW phasing, or qubitization) at L = 8: ~1e6 physical
  qubits, seconds-to-minutes runtime, T count roughly L-independent when ε ∝ L².
- Hardware table (Table I) and the surface-code law μ d² Λ^{-(d+1)/2} are a
  ready-made second parameter set for sensitivity checks.

**`refs/PRXQuantum.5.010337.pdf` (Akahoshi et al. 2024, STAR)**:
- Rotation via injected |m_θ⟩ = R_Z(θ)|+⟩ into a [[4,1,1,2]] → d=2 → d surface
  code, consumed by lattice-surgery Pauli-product measurement, repeat-until-success
  (mean 2 tries).
- Leading-order error **angle-independent**: P_Z = 2p/15 per attempt, ≈ 4p/15 per
  rotation incl. RUS; X error O(p²).
- Rotation error is a pure phase flip → exact PEC; overhead for N rotations
  **γ^{2N} ≈ exp(16 p N / 15)**. Budget N_rot ≈ 15/(4p) (e⁴ ≈ 55 overhead).
- Surface-code memory fit under their circuit-level noise: P_L = C (p/p_th)^{(d+1)/2}
  with C ≈ 0.07, **p_th ≈ 0.0039** (so p = 1e-3 is only 4× below threshold).
- Time unit "1 clock = d rounds"; one analog Pauli rotation ≈ 18 clocks average.
- 2D Hubbard benchmark: JW snake, **11N rotations per Trotter step** (N sites).
  Their device: 1e4 qubits at p = 1e-4 → 64 logical qubits at d = 7, 3.75e4
  rotations, 1.7e7 Cliffords.

**Later STAR variant** (Toshio et al. PRX 15, 021057 (2025); used by 2603.13093):
transversal injection gives per-rotation error **∝ angle**, p_analog ≈ 0.4 k |θ| p
(k = ⌈d/2⌉), plus a code-growth floor 2·I_growth(d_i, d_f) that is ~1e-5–1e-8
depending on d_i. PEC overhead becomes exp(4 α_RUS p θ_tot) with θ_tot ≈ T‖H‖₁/2.
**Under this model refining the Trotter step is free at leading order.**

**`refs/CUDA-Q Logical...pdf`**: MLIR compiler (P0 logical → P1 placement → P2 QEC →
P3 physical schedule) with tiered estimators (`LogicalProfile`, `FabricCounts`,
`FabricEstimate`, `ScheduleEstimate`). Inputs are CUDA-Q kernels (constant angles)
or its Python DSL; **not OpenQASM**. Non-Cliffords only via Clifford+T RUS synthesis;
**no partial FT, no analog injection, no PEC, no NISQ mode**. Has a worked
plaquette-Trotter Hubbard PE example (L = 8/16/32 on Pinnacle qLDPC; L = 16: 513
logical qubits, 7.4e6 T states, ~300 s). Physical knobs: p_phys, scaling law
(prefactor, threshold), cycle time, reaction/decode time, failure budget. Gate vs
measurement time are not separately exposed; only a cycle time.

## 2. Tool status on Improv (checked today)

| Tool | Status |
|---|---|
| CUDA-Q (`~/soft/cuda-quantum`) | upstream main `v0.15.2_base-247`, clean, **not built**; no `cudaq` installed in any env |
| CUDA-Q Logical | in `preview/logical/` (`import cudaq.logical`); needs py ≥ 3.11 + `pip install cuda-quantum-cu12` + CMake build. Ships `cudaq/logical/algorithms/fermi_hubbard.py` (2D JW first-order Trotter, rotation classes with multiplicities), `examples/03_fermi_hubbard.py`, `examples/02_surface_code_resource_estimate.py` (`surface_physical_target(distance, p_phys, failure_budget, cycle_time)` → physical qubits, makespan) |
| Base CUDA-Q extras | `cudaq.estimate_resources`, `compiler-bench-ftqc-clifford-t` target (gridsynth T counts), `cudaq.synth.gridsynth`, Kraus noise models (no gate times) |
| pyLIQTR 1.4.2 | works in `~/soft/mambaforge/envs/qb_gsee` (py 3.12, qualtran 0.4.0, cirq 1.4.0). `FermiHubbard(shape, J=-t, U)`, `FermiHubbardSquare` encoding (Babbush 2018), `qsvt_dynamics`, `estimate_resources`, `schedule_circuit`. **No Trotter builder, no surface-code/physical layer.** Rotation synthesis law: T ≈ 3.02 log2(1/ε) + 0.77 |
| Azure QRE (`qsharp`) | also in `qb_gsee`; takes logical counts + qubit params → physical qubits, runtime. Fastest full-FT path |
| stim / pymatching / sinter | `~/soft/mambaforge/envs/stim`, `~/soft/venvs/stim314_ns` |
| GPU | none on login node; everything above is CPU |

## 3. Problem definition (fix on day 1)

- Lattice L×L, N = L² sites, 2N qubits (JW, snake or fswap network), open BCs.
- H = −t Σ (c†c + h.c.) + U Σ n↑n↓, U/t = 4. In Pauli form: hopping (XX+YY)/2 with
  Z strings, onsite ZZ and Z. Per Trotter step: 11N rotations (STAR count),
  ~10N two-qubit gates in NISQ.
- Dynamics: quench from checkerboard (CUDA-Q module's default), evolve to time T
  (units 1/t), measure a **local** observable (site density or double occupancy)
  to additive precision ε = 0.01.
- Second-order Trotter, r steps. For a local observable the error is L-independent
  → r ∝ T^{3/2} ε^{-1/2}. Calibrate the prefactor on 2×2, 2×3 exactly.
- **Decouple device size M×M (physical) from problem size L×L.** "Same sea"
  means: for a given M, what (L, T, ε) is reachable by each architecture in a week.
- Hardware: t_g = 10 ns, t_meas = 100 ns, p = 1e-3 (also sweep 1e-4), 1-week
  cutoff (6.05e5 s). Surface-code cycle ≈ t_meas + 4 t_g + reset ≈ 150–200 ns.

## 4. Cost model per arm (the laws to test)

Let S = shots achievable in a week, ε target, α = effective exponent per error
location under the mitigation strategy (α ≈ 1 for post-selection of detectable
errors, α ≈ 4 for PEC γ² = e^{4pG}).

**NISQ**: qubits = 2N; shot time ≈ r · c · L · t_g + t_meas (JW strings / fswap
depth ∝ L); samples = e^{α p G}/ε² with G = total error locations ≈ 10N·r (+ idles).
Wall: **G ≲ ln(S ε²)/(α p) ≈ 1e4 at p = 1e-3**, logarithmic in the time budget.
Reachable region: L² T^{3/2} ε^{-1/2} ≲ const. Device size enters only as
L ≤ M/√2.

**STAR, 2024 model**: qubits ≈ (2N..4N patches + routing)·2d²; d from the Clifford
volume condition N_patch · N_clock · P_L(d) ≤ 0.01. Rotation budget
**N_R = 11 N r ≲ 15 ln(S ε²)/(16 p)** — only ~4× more Trotter steps than NISQ at
the same p, for 2d² more qubits per site. Time: 11N r / (parallelism) × 18 clocks ×
d × t_cycle per shot, × PEC overhead.

**STAR, small-angle (2025) model**: PEC overhead exp(4 α_RUS p θ_tot + 8 N_R I_growth),
θ_tot ≈ T‖H‖₁/2, ‖H‖₁ ≈ 4N t + N U/4. Wall: **T ‖H‖₁ ≲ 2/(α_RUS p)** (a hyperbola
T L² ≲ const in the (L,T) plane), independent of ε to leading order, plus the
growth-floor term which does depend on r. Needs large d_i/d_f → big qubit cost.

**Full FT**: N_T = N_R × (1.15 log2(1/ε_synth) + 9.2) [or pyLIQTR's 3.02 log2 + 0.77];
ε_synth ≈ ε/N_R. d from N_ops · P_L(d) ≤ ε_budget with P_L = 0.1 (p/0.01)^{(d+1)/2}
(also try the STAR paper's C = 0.07, p_th = 0.0039 — at p = 1e-3 that moves d from
~15 to ~25). Qubits = (2N + ancilla) · 2d² · routing + factories; time = max(T-limited,
Clifford-limited). No exponential wall; large fixed cost.

**Hypotheses stated as laws**
1. NISQ ceiling in total gates ∝ (1/p) ln(budget): weekly vs daily budget changes it
   by ~10 %. Exponential sampling wall is real and nearly time-budget-independent.
2. STAR window at p = 1e-3 is thin (rotation error 4p/15 vs 2q gate p, threshold
   0.4 %); it opens at p = 1e-4. The "angle accuracy" exponential is really
   exp(p × Σ|θ|) or exp(p × N_R) depending on protocol — locate it by sweep.
3. Full FT crosses over STAR at N_R ~ 1e5–1e6 rotations (consistent with 2603.13093)
   and NISQ at G ~ 1e4; the qubit cost of that crossover is 1e5–1e6 physical qubits
   at p = 1e-3.

## 5. Work packages

1. **Cost model** (`model.py`): three functions `nisq(L,T,eps,p,M,budget)`,
   `star(...)`, `ftqc(...)` returning physical qubits, wall time, sampling
   overhead, feasibility. Shared Trotter step count and Hubbard term counts.
   Do this first; everything else supplies constants.
2. **NISQ calibration**: stim or CUDA-Q density-matrix on 2×2 / 2×3 under
   depolarizing noise; measure bias and variance vs pG under particle-number
   post-selection and PEC/ZNE → pin α. Same runs check the Trotter constant.
3. **STAR module**: rotation/Clifford counts from the term list (or CUDA-Q's
   `fermi_hubbard_stats`), both error models, d from the Clifford-volume
   condition, clocks → seconds with t_cycle.
4. **Full-FT arm**: (a) Azure QRE on logical counts (T, rotations, logical qubits,
   depth) with our qubit params — hours, no build. (b) In parallel, build
   CUDA-Q Logical (`module load python/3.11.4`, venv, `pip install cuda-quantum-cu12
   cudaq-devel nanobind cmake ninja`, `cmake -S preview/logical -B build/...`), run
   `fermi_hubbard_jordan_wigner(L, trotter_steps=r, evolution_time=T)` through
   `surface_physical_target` as the schedule-level cross-check. (c) pyLIQTR
   `qsvt_dynamics` on `FermiHubbardSquare` for an FT-only qubitization curve
   (queries ≈ 1.4 α T + ln(1/ε), α = 2N|t| + N U/2); qubitization angles are
   large and are not STAR-compatible, so this stays FT-only.
5. **Sweeps and figures**:
   - cost vs L at fixed (T, ε): physical qubits and wall time, three arms.
   - **phase diagram in (L, T) at ε = 0.01** for M ∈ {1e2, 1e3, 1e4, 1e5}: reachable
     regions and crossovers, 1-week cutoff.
   - same at p = 1e-4 (STAR window should open).
   - Trotter vs qubitization crossover inside the FT arm.

Keep second-order Trotter as the common algorithm across all three arms so the
comparison is apples-to-apples.

## 6. Risks / open choices

- CUDA-Q Logical build (half a day; QRE fallback).
- Which STAR error model to headline: model both, report both.
- Surface-code threshold assumption (1 % vs 0.4 %) moves d by ~10 at p = 1e-3 and
  is the biggest single knob in the FT and STAR qubit counts.
- Trotter error bound for local observables: use empirical small-L calibration plus
  the light-cone argument, not the global commutator bound (which is ∝ N and would
  wrongly make r grow with L).
- JW layout on a 2D grid: snake strings (depth ∝ L per step) vs fermionic swap
  network vs local encodings (Derby–Klassen, 1.5× qubits). Affects NISQ depth and
  STAR Clifford count; pick snake + strings for simplicity and note the alternative.
