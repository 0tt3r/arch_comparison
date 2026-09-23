# Working plan: protection, mitigation, and algorithm choice for Fermi-Hubbard dynamics on a common qubit array

2026-09-22. Merged working plan. Sources: PLAN_CODEX.md, PLAN_CLAUDE.md, PLAN_CODEX2.md,
PLAN_CLAUDE2.md, and the two rounds of cross-critique. Framing and organizing questions
follow PLAN_CODEX2.md; conventions table, per-arm models, and constants follow
PLAN_CLAUDE2.md with the corrections from the second Codex critique applied. No resource
curves or crossover locations have been established. Every number is a reference
constant with provenance or arithmetic on one.

---

## 1. Purpose

Combine the perspectives of application, mitigation, partial-FT, and full-FT experts to
identify what governs the resource cost of useful Fermi-Hubbard dynamics on a common
physical-qubit array, and to find better ways of spending that array than the three
reference approaches:

1. Physical qubits with error suppression, post-selection, and mitigation.
2. STAR: surface-code Cliffords, injected analog rotations, PEC.
3. Full surface-code FT with Clifford+T synthesis and magic-state supply.

Their ordering, and whether STAR has a useful window, are open. Intermediate strategies
(different protection per operation class, hybrid analog/distilled rotations,
observable-aware allocation of correction and mitigation) are in scope when a mechanism
motivates them and their implementation cost is included. Scope is set by scientific
value, not by software effort or elapsed time.

## 1a. Stipulated hardware and variables

- Square 2D nearest-neighbor array of Q = M² physical qubits (Q is the swept device size).
- t_g = 10 ns (1q and 2q gates), t_meas = 100 ns, reset 100 ns (declared, C6), p = 1e-3
  (nominal, swept to 1e-4), B = 604800 s answer-time cutoff including sampling.
- N = L² sites, 2N modes before overhead, τ = Jt/ħ, ε absolute observable error,
  δ estimator failure probability, T_answer elapsed time to the specified answer.

## 2. Organizing questions

| Question | What resolves it |
|---|---|
| How much of the system must remain accurate for the requested answer? | Observable support, exact circuit causality, physical propagation with a budgeted truncation error, the region affected by relevant noise |
| Where does error correction buy the most useful computation? | Reduction in estimator cost vs. added area, execution time, lost replica throughput |
| What limits STAR for dynamics? | Total analog angle, angle-independent floors, preparation acceptance, resource-state delivery rate, Clifford volume |
| What makes full FT favorable? | Protected volume, synthesis and factory demand, longer accessible evolution, or a different algorithm |
| Which apparent laws survive implementation changes? | Controlled changes of mapping, ordering, boundary conditions, observable, error model |
| Can a tailored intermediate strategy beat all three? | Optimize protection and rotation choices at fixed task and error target |

The expert contributions meet at interfaces: the application fixes which errors matter,
mitigation fixes the attainable bias–variance trade-off, QEC fixes the cost of
suppressing particular errors, compilation fixes whether the assumed operations and
parallelism are realizable. The consequential assumptions live at these interfaces.

## 3. The candidate unifying picture (hypothesis, with its failure modes as results)

In each regime the dominant cost is set by a noise exposure χ = p × (weighted volume),
where only operations that influence the observable count. NISQ pays per gate, STAR
pays per radian, FT pays per logarithm:

| Regime | Exposure | Enters as | Weight |
|---|---|---|---|
| NISQ + PEC | χ_phys = Σ_j c_j p_j over included locations | Γ² = e^{χ_phys} (c_j ≈ 4 for PEC on 2q depolarizing) | per gate |
| STAR small-angle | χ_angle = 4 p Σ_j α_j \|θ_j\| | Γ² = e^{χ_angle} | per radian |
| STAR floor | χ_floor = 4 Σ_j q_floor,j | Γ² = e^{χ_angle + χ_floor} | per rotation |
| STAR original | χ = (16/15) p N_rot | Γ² = e^{χ} | per rotation |
| Full FT | protected volume V_L, synthesis N_rot log(N_rot/ε_synth) | d ∝ log(V_L/q_L)/log(p_th/p); area ∝ d², time ∝ d | logarithmic |

Positive-time first- and second-order product formulas preserve Σ\|θ_j\| = τ Σ_l \|h_l\| as
the step count changes, so refining the Trotter step multiplies χ_phys by r, leaves
χ_angle fixed, raises χ_floor, Clifford work, preparation demand, and execution time,
and only raises FT T-count and time. Higher-order formulas with signed stages do not
preserve the total angle and must be recounted (relevant because Chung et al. suggest
higher-order Trotter to shrink angles). A fixed α and fixed q_floor is a slice through
an optimization over protocol and distances, not an intrinsic property of STAR.

Whether these variables organize the results is the test in §7. Do not assume collapse.

## 4. Conventions (settle in writing, one owner each, before any constant is fitted)

Proposed defaults in bold; alternatives are sensitivity parameters.

| # | Convention | Proposed default | Why it matters |
|---|---|---|---|
| C1 | Hamiltonian, geometry | H = −J Σ(c†c + h.c.) + U Σ n↑n↓, **U/J = 4, even L×L, open BC, zero flux** | Open vs. periodic changes term counts and JW strings; the experiment and Toshio use a torus |
| C2 | Initial state | **Half-filled Néel product state**; dimer/holon-doublon state as extension | Both are constant-depth preparations whose reset, gates, and errors are counted; the dimer state is fermionic-magic (3 Rzz per dimer) |
| C3 | Observable | **Local doublon density D_i at a central site** | Local charge density is pinned by particle-hole and spin symmetry for this quench |
| C4 | Accuracy | **Absolute error ≤ 0.01 at 95 % CI on total error**: statistical + truncation + Trotter + synthesis + mitigation bias + logical-fault bias + SPAM; δ = estimator failure probability named explicitly; split allocated and optimized per arm | A statistical interval alone is not the criterion; FT still pays sampling |
| C5 | Angle convention | **U_P(θ) = e^{−iθP}; Θ = Σ\|θ_j\| = τ Σ_l \|h_l\|** before Clifford extraction. Converting the same gate to e^{−iφP/2} gives φ = 2θ and halves the angle-error coefficient. Toshio's e^{iθZ} has the same magnitude as ours; α_RUS ≈ 1.5 is quoted in it | A factor of two in α is a factor of two in the STAR wall |
| C6 | Noise model | Circuit-level stochastic Pauli, one scalar p: **2q depolarizing p, 1q depolarizing p/10, init and readout flip p, idle depolarizing p per idle *layer* (Stim circuit-level convention), reset 100 ns, feedforward 200 ns**. Sensitivity: time-proportional idle at rate p per 100 ns, under which a data qubit waiting through a 200 ns measure+reset interval accumulates ≈ 2e-3, i.e. a cycle-level exposure above the nominal threshold. Coherent/correlated components added where they could change the inference | Idle noise during measurement is the most regime-asymmetric assumption and the one that can put the operating point above the effective threshold; prep/readout must not silently inherit the gate convention |
| C7 | Syndrome cycle | **Explicit schedule: 4 CZ layers + measurement + reset ≈ 240 ns**, plus decode/feedforward latency where it gates a logical operation; one clock = d cycles | Cannot be 10 ns; Toshio uses 1 μs, TopQAD 25 ns gates with 100 ns readout |
| C8 | Suppression law | **One law for both encoded arms**, p_L(d) = A (p/p_th)^{(d+1)/2} per patch per cycle. **A = 0.1, p_th = 0.01 is a reference curve only** until §7.3 calibrates (A, p_th) under the C6/C7 schedule; sensitivity A = 0.07, p_th = 0.0039 (Akahoshi circuit-level fit). State the interval convention. The calibration must allow the outcome that at this operating point increasing d does not help | C6–C8 are not yet a consistent calibrated pair; this is the most consequential remaining assumption. The fits differ by ~10 in d at p = 1e-3; different fits per arm manufacture an advantage |
| C9 | Logical-fault criterion | **Per-shot logical fault probability q_L budgeted inside C4 as bias ≤ q_L for the ordinary estimator and ≤ γ q_L under PEC** (quasi-probability norm amplifies worst-case bias). Not zero faults over the campaign | Zero-fault-campaign inflates d by 2–4 |
| C10 | Fermion mapping | **JW snake with Z strings** baseline; **fermionic swap network** (merged hop+FSWAP) as locality-conscious comparison; compact local encodings noted with their extra qubits and checks | A mapping artifact must not be read as an architecture law |
| C11 | Locality | Two distinct constructions, both budgeted: (a) **exact circuit cone** of the observable, which depends on r, ordering, mapping, routing and is not shrunk by small angles; used for exact locality-aware PEC. (b) **Physical truncation** to a region outside which influence is weak, with an explicit truncation error ε_trunc inside C4. Available to all three arms | The two give different resource estimates; conflating them was an error in PLAN_CLAUDE2 |
| C12 | Trotter order and step law | **Second order**; r(τ, ε) fitted for the doublon observable on 2×2, 2×3, 3×3, 4×3 exact dynamics with a stated extrapolation, supported by analytical bounds and convergence evidence. Not the global commutator bound | The τ exponent is unknown a priori because support grows with τ |
| C13 | Replicas and area | **m = ⌊Q / Q_copy⌋ subject to actual placement, every arm**; encoded arms trade area among data, routing, factories, injection patches, replicas | Device size enters the NISQ budget as ln(m S₁/N₀) |
| C14 | Size/time families | (1) **fixed τ ∈ {1, 2, 4}** (local physics; optimized cost may saturate in L); (2) **τ = cL** (system-scale dynamics; relevant region grows); (3) an extended observable or profile if needed to test a mechanism, with its own normalization and coverage criterion | Three scientifically different questions |
| C15 | Term counts | Open L×L, E = 2L(L−1) bonds: uncombined first-order 4E + 3N = 11N − 8L terms; **drop the 2N single-Z terms** (sector-dependent phase) → 9N − 8L; **λ = Σ\|h_l\| = 5N − 4L** at U/J = 4. Second-order rotation count per step is **not** the first-order count: record the executed multiplicities after ordering and merging | Quoted 11N and λ = 5N are periodic |

## 5. Accounting

Cost the production of one observable estimate. Principal outputs:

    Q_min(N, τ, ε, δ; T_answer ≤ B)      and      τ_max(N, Q, ε, δ; T_answer ≤ B)

with B = 6.05e5 s, retaining the area/runtime Pareto frontier and its bottleneck
(sample throughput, circuit latency, resource-state supply). Where acceptance and
variance factorize,

    T_answer ≈ N₀(ε_stat, δ) · Γ² · t_attempt / (a · m)

N₀ = ideal shots (normal-approximation reference: 1.96² · (1/4) / 0.005² ≈ 3.8e4 for
worst-case Bernoulli variance at ε_stat = 0.005 and 95 %; a smaller specified variance or
a different statistical allocation lowers it; not a finite-sample guarantee), Γ² =
variance multiplier, a = whole-circuit acceptance, m = replicas. When post-selection and
mitigation interact, use the actual conditional estimator. Factory rejection and RUS
retries live in the supply model, not as circuit rejection. Include setup, calibration,
training, and online classical latency where material, and state the boundary.

The mitigation wall under an exponential model: χ ≲ ln(m S₁/N₀), with **S₁ computed from
each arm's own attempt time** (a STAR attempt is 1e3–1e5 times longer than a NISQ
attempt, so its S₁ and its allowance are correspondingly smaller and depend on N_rot
itself). This gives a logarithmic gain from budget and replicas once mitigation
dominates; it is not a universal gate ceiling.

## 6. Per-arm models

Common: N = L²; N_cone(τ) from C11(a) or C11(b); r(τ, ε) from C12; λ_cone from C15
restricted to the region; N_rot from executed multiplicities.

**Physical qubits.** Width 2N (+ swap-network ancilla). t_attempt = r layers × depth ×
t_g + t_meas, order 10 μs. Strategy menu, each with its own (a, Γ², bias): raw; symmetry
post-selection on (N↑, N↓) (acceptance a ≈ e^{−f p G} with f the detectable fraction,
plus residual bias from undetected faults; hurt at 160 shots in the experiment, may help
at 1e4); whole-circuit PEC (Γ² ≈ e^{4pG}, upper envelope); cone-restricted PEC
(Tran–Sharma–Temme); training-based (TFLO) + smoothing as an **explicitly optimistic
comparison curve**: it enters the claimed Q_min frontier only once its residual bias is
shown to meet C4 with stated evidence. Executed calibration point: Phasecraft/Quantinuum
7×4 torus, 56 qubits, 2415 two-qubit gates at 2q error 1e-3 (pG ≈ 2.4, XEB 1–5 %),
160 shots/point, pseudo-twirling + spin-reflection averaging + TFLO + GPR, no
post-selection, 19 % of shots at correct Hamming weight, ~5 s/shot. Their accuracy
statement: doublon number and weight-≤2 Z/ZZ observables at U = 0 against the exact
NearFLO reference, mitigated-minus-exact residuals ≈ 0.05 read from their Figs. 21–23,
with uncertainties per their Eq. D1 (twirl-instance total variance); no U = 4 exact
reference exists. Recompile for the NN array (their SWAPs are free relabelings) before
comparing.

**STAR.** Patches ≈ 1.5(2N + 6) sequential (Toshio floor plan) or Litinski compact block
(Akahoshi), plus injection patches (2–3 per resource state per clock at p = 1e-3, Toshio
Fig. 10), plus parallel-rotation layouts (Akahoshi et al. arXiv:2408.14929). Rotation
error, RUS folded in: original 4p/15 angle-independent; Toshio P̃_L = α_RUS|θ|p,
α_RUS ≈ 1.5 with protocol switching, no injection floor at full distance, coherent
miscalibration Σφ²/2 to be budgeted; TopQAD 0.4k|θ|p + 2 I_growth(d_i, d_f) with
I_growth = μ d_i² Λ^{−(d_i+1)/2}. **The operating point (m, k, d_i, d_f, growth,
acceptance, delivery rate) is an optimization output** for our workload, budget, and
throughput, not imported from Toshio's QPE table. Prediction to check: full-distance
injection at k = d_f/m has ideal acceptance p_ideal = sin^{2k}θ + cos^{2k}θ with physical
angle θ = θ*^{1/k}; this tends to 1 as θ* → 0 at fixed k, but over our workload's angle
distribution (θ* ≈ JΔτ and UΔτ/4, i.e. 0.01–0.5 rad for r ≈ 4–40 steps per unit τ) it is
small for k ≳ 5 (e.g. θ* = 0.01, k = 12 gives θ ≈ 0.68 rad and p_ideal ≈ 3e-3). Where
growth becomes favorable is therefore decided by the optimization over the angle
distribution, required d_f, and preparation capacity; the expectation is that at
p = 1e-3 the floor question is live. Clifford budget N_patch · N_clock · C · p_L(d) ≤ q_L (C9) sets d_f. Time per
rotation: C_av ≈ 4.6 clocks sequential (Toshio) or 18 (compact + RUS), × d × t_cycle.
Floor dominates the exponent when N_rot q_floor ~ α p Θ; N_rot ~ 1/q_floor marks
order-one accumulated floor error, a different point. Qubitization on STAR is expensive
(Toffolis → π/4 analog rotations), not excluded.

**Full FT.** Logical qubits 2N + O(1) (Trotter) or 2N + O(log N) (qubitization).
Synthesis default: **deterministic gridsynth (Ross–Selinger), ≈ 3 log₂(1/ε_synth) + O(1)
T per rotation**, which is what the inspected CUDA-Q `clifford-t-synthesis` pass and the
pyLIQTR fit (3.02 log₂ + 0.77) implement; RUS synthesis (≈ 1.15 log₂ + 9.2, CUDA-Q
Logical paper) is a separate implementation with its own scheduling cost. ε_synth =
ε_syn,tot/N_rot. d from V_L p_L(d) ≤ q_L. Factories: 15-to-1 at one or two levels at
p = 1e-3, or cultivation (≲ 1e5 T at 1e-9 per T, Chung Fig. 16); rejection and RUS in
throughput. Time = max(T-supply, Clifford, reaction) × N₀/m. Qubits = 1.5 × 2N × 2d² +
factories. Coherent estimation is a separately costed extension if precision becomes
the driver.

**Marginal criterion for any hybrid or added protection** (fixed N₀, factorized model):

    Δ log t_attempt + Δ log Γ² − Δ log a − Δ log m < 0

including conversion, routing, and synchronization costs.

## 7. Tests of the picture (primary result)

Plotting costs computed from Γ² = e^{χ} against χ illustrates the assumed law. The test
uses **independently obtained** quantities:

1. Estimator variance and bias from noisy density-matrix or trajectory simulation of
   the actual estimator (post-selection, PEC, cone-restricted PEC) on 2×2 … 4×3
   workloads under C6, varying r, ordering, mapping, τ. Plot measured Γ² and bias
   against χ_phys; the model curve is overlaid, not fitted.
2. STAR rotation and growth channels under C6/C7: the growth and memory parts are
   stabilizer and run in Stim directly; the angle-dependent injection channel is
   reduced to Clifford simulation by syndrome-subspace sampling (Toshio App. C.2,
   Eqs. C3–C5; Chung App. A.3), which ordinary Stim alone does not provide. Measured
   P̃_L(θ, d_i, d_f) and acceptance against the α\|θ\|p + floor form.
3. Logical memory and lattice-surgery p_L(d) from Stim under C6/C7 against C8.
4. Hold-out workloads outside the calibration set (different observable, initial state,
   boundary condition, mapping) to test whether the same χ organizes them.
5. The executed 7×4 point as an external check on the NISQ favorable envelope.

Then the resource outputs:

- Q_min vs N under the week cutoff for the C14 families; τ_max at fixed Q.
- Reachable (L, τ) regions at Q ∈ {1e2 … 1e6}, three arms overlaid, implementation and
  algorithm labeled.
- Decomposition of cost near each crossover: sampling, rejection, routing, protection,
  synthesis, supply.
- Sensitivity strips: p ∈ {1e-4, 1e-3}, C8 law, floor model, mapping, PEC coefficient,
  replicas, cone construction (C11a vs C11b), signed-stage Trotter.
- At least one concrete protocol or algorithm change that moves a boundary, if supported.

Testable mechanism questions: does STAR fail first from total angle, floor, or delivery
rate? Does FT gain most from longer evolution, precision, or algorithm change? Does
physical-qubit computation stay competitive through observable structure after a
full-circuit estimate predicts failure? Does the three-regime picture fail for a local
task and return for an extended observable?

## 8. Workloads

| Purpose | Workload |
|---|---|
| Main study | C1–C3 baseline, L ∈ {4, 6, 8, 10, 12, 16}, C14 families |
| Executed reference | Phasecraft 7×4 torus, π flux, dimer state, U/J = 4, τ ≤ 2, recompiled for the NN array; undertaken when it resolves a question in §7 |
| Initial-state sensitivity | Square-lattice dimer/holon-doublon quench |
| Algorithm sensitivity | Qubitized dynamics (pyLIQTR `FermiHubbardSquare`; α = 4N\|J\| + N\|U\| in site count, periodic; queries ≈ 1.4 α τ + ln(1/ε)); matched and optimized comparisons reported separately |

## 9. Tools (calibrate constants; none is the model)

| Tool | Status (Improv, 2026-09-22) | Role | Adaptation required |
|---|---|---|---|
| Cost model `model.py` | to write | The deliverable; consumes workload summaries (width, depth, angle histograms, Θ, Clifford volume, placements) and returns feasible configurations and bottlenecks | — |
| Exact small-L dynamics | `qutip` env, numpy; openfermion in `qb_gsee` | C12, C11, §7.1 | Stim cannot run generic Hubbard rotations |
| Stim + PyMatching (`stim` env) | works | §7.2–7.3 | Stabilizer parts only |
| CUDA-Q Logical (`~/soft/cuda-quantum/preview/logical`) | source only, not built; py ≥ 3.11 + `cuda-quantum-cu12` + CMake | FT logical counts, gridsynth T counts, physical tiers | `fermi_hubbard_jordan_wigner` is first-order, controlled, periodic; `surface_physical_target` is a compact reference realization with native multi-Pauli measurements, not our NN array |
| Base CUDA-Q | not installed | `cudaq.synth.gridsynth`, `compiler-bench-ftqc-clifford-t` | pip install; CPU |
| pyLIQTR 1.4.2 (`qb_gsee`) | works | Qubitization counts, scheduler | No Trotter builder, no physical layer |
| Q# resource estimator (`qsharp` 0.28 in `qb_gsee`) | installed, not validated | Independent FT estimate | Compare budgets, connectivity, synthesis, factory assumptions before reading agreement |

Agreement between packages under the same assumptions is a consistency check, not
validation. A build need not hold up analytical work.

## 10. Sequence (logical order, not a schedule; concurrent where assumptions allow)

1. Agree the physical question and C1–C15 in writing, with owners.
2. Populate reference models: both STAR generations, the physical-qubit strategy menu,
   FT reference; identify the parameters the conclusion is most sensitive to.
3. Construct common workloads and alternative mappings; verify conventions, angle
   totals, multiplicities; produce fixed-circuit and optimized-step comparisons.
4. Resolve consequential uncertainties by the cheapest adequate method: mitigation
   bias–variance, routed schedules, syndrome timing, STAR preparation/growth, FT supply.
5. Run §7; optimize and challenge the apparent regimes; test a hybrid if a specific
   opportunity appears; revisit any law that fails to organize the data.
6. Extract the general result: controlling variables, assumptions, supported
   improvements, and the boundaries still sensitive to unresolved choices.

## 11. Reference constants (inputs with provenance; re-derive under C5–C8)

| Constant | Value | Source |
|---|---|---|
| Original STAR rotation error | 2p/15 per attempt, ≈ 4p/15 per RUS rotation | Akahoshi 2024 Eqs. 6, 12–15 |
| Original STAR PEC | e^{16pN/15} | Akahoshi 2024 Sec. IV.B |
| Compact-block timing | ≈ 18 clocks per rotation | Akahoshi 2024 p. 16 |
| Akahoshi surface fit | C_Z = 0.0679, p_th,Z = 0.00385; C_X = 0.0819, p_th,X = 0.00416 | Akahoshi 2024 Table I |
| Small-angle error | P̃_L = α_RUS\|θ\|p, α_RUS ≈ 1.5 (0.4k fixed-weight before switching), e^{iθZ} convention | Toshio 2025 Eq. 25, Sec. IV.D |
| Small-angle PEC | e^{4α_RUS Θ p}; Trotter τ ≲ 1/(α_RUS λ p) | Toshio 2025 Eqs. 29–30, 37 |
| P_ud | (k/15)p native ZZ; (2k/15)p CNOT–Rz–CNOT | Toshio 2025 p. 10 |
| Sequential clocks per rotation | C_av ≈ 4.6 | Toshio 2025 Sec. VII.C |
| Toshio Hubbard 8×8 QPE, ε = 0.01, PBC | p = 1e-3: d = 25, 2.5e5 qubits, 3.3e4 h; p = 1e-4: d = 13, 6.8e4 qubits, 254 h | Toshio 2025 Table II (their workload, not ours) |
| Growth floor | 2 I_growth, I_growth = μ d_i² Λ^{−(d_i+1)/2}; 1.7e-5 (d_i = 9) … 4e-8 (d_i = 15) at 2q 5e-4 | Chung 2026 Eq. 5, Fig. 10b |
| TopQAD favorable region | 1e5–1e6 small-angle rotations, 100 logical qubits, θ ≤ 1e-4 | Chung 2026 Sec. IV.D (reference result, not our crossover) |
| Full-FT Hubbard 8×8, ε = 0.01, qubitization | ~1e6 qubits, 4.7 min time-optimal | Chung 2026 Table V |
| Fowler–Gidney fit | 0.1(100p)^{(d+1)/2} | Toshio p. 27; CUDA-Q Logical |
| T per rotation | 3.02 log₂(1/ε) + 0.77 (gridsynth fit); 1.15 log₂(1/ε) + 9.2 (RUS) | pyLIQTR `rotation_gates.py`; CUDA-Q Logical Eq. 2 |
| Executed NISQ point | 7×4, 2415 2q gates, p₂ = 1e-3, 160 shots, ≈ 0.05 accuracy, 4 second-order steps to τ = 2 | 2510.26300 Tables I–III, App. C.3, Figs. 21–23 |

## 12. References

- R1 `refs/PRXQuantum.5.010337.pdf` Akahoshi et al. 2024 — original STAR.
- R2 `refs/PhysRevX.15.021057.pdf` Toshio et al. 2025 — small-angle STAR, total-angle bound.
- R3 `refs/2603.13093v1.pdf` Chung et al. 2026 — STAR with growth vs. TopQAD full FT; ground-state task.
- R4 `refs/CUDA-Q Logical ... .pdf` — compilation-linked FT estimation.
- R5 `refs/2510.26300v2.pdf` Phasecraft/Quantinuum 2025 — executed 7×4 quench.
- arXiv:2303.06496 Tran, Sharma, Temme — locality-aware mitigation.
- arXiv:1912.08854 Childs et al. — Trotter error theory.
- Takagi et al., npj QI 8, 114 (2022) — mitigation limits under stated assumptions.
- arXiv:2408.14929 Akahoshi et al. — parallel STAR rotations.
