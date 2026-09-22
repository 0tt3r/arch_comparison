# Synthesis plan: what controls the cost of Fermi-Hubbard dynamics on a common qubit array, across NISQ, STAR, and full FT

Written 2026-09-22 by Claude (Fable 5.1). Synthesizes PLAN_CLAUDE.md, PLAN_CODEX.md, the
Codex critique of PLAN_CLAUDE.md, and all five references in `refs/`. Supersedes
PLAN_CLAUDE.md. Audience: the authors of the literature being synthesized. No numerical
results have been produced; every number below is either a reference constant (flagged
as such) or arithmetic on one.

---

## 1. The question, stated once

For a square nearest-neighbor array of Q physical qubits with t_g = 10 ns, t_meas = 100 ns,
p = 1e-3, and one week of wall time including sampling: what Fermi-Hubbard dynamics
(L×L sites, evolution time τ = Jt, local observable to absolute error 0.01) can be
computed by (i) physical qubits with error suppression, post-selection, and mitigation,
(ii) STAR (surface-code Cliffords + injected analog rotations + PEC), (iii) full
surface-code FT with Clifford+T synthesis and magic-state supply?

The deliverable is not the three cost curves. It is the identification of the
**controlling quantity in each regime**, a test of whether cost collapses onto it, and
the location of the crossovers as a consequence. Everyone in the room already knows the
per-regime story; the new content is the joints between the stories, where each
community's default convention is invisible to itself and wrong for the others.

---

## 2. The unifying picture (hypothesis to be tested, not assumed)

In every regime the dominant cost is set by a **noise exposure**: p times a weighted
volume of the circuit, where only the part inside the observable's backward light cone
counts. The regimes differ in the weight per operation and in whether the exposure
sits in an exponent or a logarithm.

| Regime | Exposure X | Where X enters | Weight per operation |
|---|---|---|---|
| NISQ + PEC | X = p · Σ_{gates in light cone} w_g, w_g ≈ 1 per 2q gate, ~0.1 per 1q | Γ² = e^{αX}, α ≈ 4 for PEC, smaller for post-selection/training envelopes | every gate, once |
| STAR (small-angle, Toshio 2025) | X = p · Σ_{rotations in light cone} α_RUS|θ_i| = α_RUS p Θ_eff, α_RUS ≈ 1.5 | Γ² = e^{4X} (e^{2X} with symmetric Clifford twirling) | rotation by |θ| |
| STAR + growth floor (TopQAD 2026) | X = α_RUS p Θ_eff + N_rot,eff · q_floor, q_floor = 2 I_growth(d_i, d_f) | Γ² = e^{4X} | rotation by |θ| plus a constant |
| STAR (original, Akahoshi 2024) | X = (4/15) p N_rot,eff | Γ² = e^{4X} = e^{16pN/15} | rotation, angle-independent |
| Full FT | protected volume V_L = N_patch · N_cycles; synthesis demand N_rot log(N_rot/ε) | d ≈ 2 log(V_L/ε_L)/log(p_th/p), qubits ∝ d², time ∝ d | logarithmic |

The physics in one sentence: **NISQ pays per gate, STAR pays per radian, FT pays per
logarithm.** Refining the Trotter step multiplies the NISQ exposure by r, leaves the
STAR angle exposure fixed at λτ (Toshio Eq. 37: τ ≲ 1/(α_RUS λ p) independent of step
count), and multiplies only the FT T-count and wall time. The growth floor is what
reintroduces a per-rotation charge into STAR, and whether the floor or the angle
controls the STAR wall at p = 1e-3 is the sharpest open question below.

Under this picture the three walls have different shapes in the (L, τ) plane:

- NISQ: X ∝ p c₂ N_lc(τ) r(τ, ε), where N_lc is the number of sites in the light cone
  (saturates at N = L² when the cone fills the lattice). Log-in-time-budget wall:
  X ≲ ln(m S₁ ε²)/α with m = replicas, S₁ = shots per replica per week.
- STAR: X ∝ 1.5 p λ_lc τ with λ_lc = 5 N_lc at U/J = 4 (after dropping number-conserving
  terms). Same log-in-budget wall, but r-independent.
- FT: no exponential wall; polynomial cost with a large fixed offset (d ≈ 25 at p = 1e-3
  under the Fowler-Gidney fit, i.e. ~1.9e3 physical qubits per patch after routing).

Arithmetic on reference constants, for orientation only (whole lattice, no light cone,
PEC α = 4, c₂ ≈ 20 two-qubit gates per site per second-order step from the Phasecraft
compilation): the ratio of NISQ to STAR exponents is (c₂ r)/(1.5 · 5 τ) ≈ 2.7 r/τ. At
r/τ ≈ 2–4 steps per unit time STAR gains a factor 5–10 in exponent for roughly 2d² · 1.5
≈ 2e3 physical qubits per logical qubit at d = 25. That is the trade the study
quantifies; the collapse test decides whether the exponent ratio is really set by
c₂ r versus 7.5 τ or by something else (floor, Clifford volume, acceptance).

---

## 3. Conventions to settle before any constant is fitted

These are the joints. Each is a decision one person at the table can make in minutes
and each moves a crossover by an order of magnitude. Proposed defaults in bold;
alternatives become sensitivity parameters.

| # | Convention | Proposed default | Why it matters / who owns it |
|---|---|---|---|
| C1 | Hamiltonian, geometry | H = −J Σ(c†c + h.c.) + U Σ n↑n↓, **U/J = 4, even L×L, open BC, zero flux** | Open vs periodic changes term counts (11N periodic vs 11N − 8L open, uncombined) and JW string structure. Periodic torus is what the experiment and Toshio use. |
| C2 | Initial state | **Half-filled Néel product state**; dimer/holon-doublon state (Phasecraft) as extension | Product state has zero prep cost; dimer state is fermionic-magic and has O(1) depth. |
| C3 | Observable | **Local doublon density D_i = ⟨n_i↑ n_i↓⟩ at a central site** | Local charge density is pinned by particle-hole symmetry for the Néel quench. Doublon ∈ {0,1} so variance ≤ 1/4. |
| C4 | Accuracy criterion | **Absolute error ≤ 0.01 at 95 % CI on the total error**: statistical + Trotter + synthesis + mitigation bias + logical-fault bias + SPAM. Allocate and optimize the split per arm. | A narrow statistical interval alone is not the criterion. FT still pays sampling. |
| C5 | Angle convention | **Rotations written e^{−iθP}; Θ = Σ|θ_i| = τ Σ_l |h_l| = λτ** for a positive-time product formula. λ = ‖H‖₁ in the Pauli basis after dropping conserved terms; λ = 5N at U/J = 4. | Toshio's e^{iθZ} convention is the same magnitude as this one; α_RUS ≈ 1.5 is quoted in it. The "usual" e^{−iθP/2} convention halves θ and doubles α. |
| C6 | Noise model | Circuit-level stochastic Pauli, single parameter p: **2q depolarizing p, 1q depolarizing p/10, init/meas flip p, idle error rate p per (2q gate time) applied to elapsed time**; reset time and feedforward latency declared (propose 100 ns reset, 200 ns feedforward). | Idle noise during the 100 ns measurement is a large fraction of the syndrome-cycle budget and is the single most regime-asymmetric assumption. STAR/Toshio use 1q depolarizing p, not p/10; pick one. |
| C7 | Syndrome cycle | **Explicit schedule: 4 CZ layers + measurement + reset ≈ 4·10 + 100 + 100 = 240 ns**, plus decode/feedforward latency where it gates a logical operation. One clock = d cycles. | Cannot be 10 ns. Toshio assumes 1 μs; TopQAD assumes 25 ns gates with 100 ns measurement. |
| C8 | Suppression law | **One law for both encoded arms**: p_L(d) = A (p/p_th)^{(d+1)/2} per patch per cycle, **A = 0.1, p_th = 0.01** (Fowler-Gidney, used by Toshio and CUDA-Q Logical). Sensitivity: STAR-paper fit A = 0.07, p_th = 0.0039 (their own circuit-level noise). | At p = 1e-3 the two fits differ by ~10 in d. Assigning different fits to STAR and FT manufactures an advantage. |
| C9 | Logical-fault criterion | **Per-shot logical fault probability q_L contributes bias ≤ q_L to a bounded observable; budget it inside C4** (e.g. q_L ≤ 1e-3). Not "zero faults over the campaign". | Zero-fault-campaign is a stronger, different requirement and inflates d by ~2–4. Toshio's Eq. 42 uses 1 % per deepest circuit. |
| C10 | Fermion mapping | **JW snake with Z strings** as baseline; **fermionic swap network** (Phasecraft-style, merged hop+FSWAP) as the locality-conscious comparison; compact local encodings (Derby-Klassen, 1.5× qubits) noted. | Square-lattice fermions are not automatically geometrically local. A mapping artifact must not be read as an architecture law. Affects NISQ depth, STAR Clifford volume, FT routing. |
| C11 | Light cone | **V_eff and Θ_eff count only operations in the observable's backward light cone**, defined by exact circuit causality (not Lieb-Robinson bounds), and available to all three arms. Locality-aware PEC per Tran-Sharma-Temme. | This is what makes fixed-τ cost saturate in L. It must be applied uniformly or NISQ is penalized by fiat. |
| C12 | Trotter order and step law | **Second order**, r(τ, ε) fitted on 2×2, 2×3, 3×3, 4×3 exact dynamics for the doublon observable, with an explicit extrapolation assumption. Do not use the global commutator bound (∝ N). | Phasecraft: local Trotter error roughly size-independent, 4 steps to τ = 2 comparable to hardware noise; grows for τ > 1.5. The τ exponent is not known a priori because the light-cone region grows with τ. |
| C13 | Replicas and area | **m = ⌊Q / q_arch(N)⌋ concurrent replicas in every arm**; encoded arms may trade area between data, routing, factories, injection patches, replicas. | Device size enters the NISQ budget as ln(m S₁ ε²). |
| C14 | Time families | **Fixed τ ∈ {1, 2, 4}** and **τ ∝ L** (cone crosses the lattice; v_LR ≈ 2J as a placeholder to be fitted) | Fixed-τ slices test L-saturation; τ ∝ L is where the walls actually differ. |
| C15 | Term counts | Open BC, uncombined: 8N − 8L hopping Pauli terms + N ZZ (+ 2N single-Z, **dropped** as number-conserving) = **9N − 8L** rotations per first-order pass; second-order merging halves the boundary layers. Cross-check against CUDA-Q `fermi_hubbard_stats` (periodic default). | The STAR paper's 11N is periodic and keeps the Z terms. |

---

## 4. Cost model per arm

Common: N = L², N_lc(τ) sites in the light cone, r(τ, ε) from C12, λ_lc = 5 N_lc,
N_rot,eff = (9N_lc − 8L_lc) r, Θ_eff = λ_lc τ, N₀ = (1.96)² · Var/ε_stat² ≈ 1e4 for
ε_stat = 0.005, budget B = 6.05e5 s.

Runtime factorization (used only where the factors are separable):

    T_answer ≈ N₀ · Γ² · t_attempt / (a · m)

a = whole-circuit acceptance, m = replicas, Γ² = mitigation variance multiplier.
Post-selection and mitigation interact; estimate acceptance, conditional bias and
estimator variance together on the small-L noisy calibrations rather than multiplying.

### 4.1 Physical qubits (NISQ)

- Width 2N (+ ancilla for swap network if any); q_arch = 2N.
- Attempt time: r layers of depth ≈ c_d · L (JW strings) or O(1) per layer (swap network)
  at t_g, plus t_meas. Order 10 μs.
- Strategy menu, each evaluated, with its own (a, Γ², bias):
  1. Raw, where the bias is inside budget (only for X ≲ 0.01).
  2. Symmetry post-selection on (N↑, N↓): acceptance a ≈ e^{−f p G}, f = detectable
     fraction; residual bias from undetected faults. Phasecraft found this *hurts* at
     160 shots; at our N₀ ≈ 1e4 it helps. Compute, don't assume.
  3. PEC, whole-circuit: Γ² ≈ e^{4 p G}. Reference model, upper envelope of cost.
  4. PEC, light-cone: Γ² ≈ e^{4 p G_lc}. Tran-Sharma-Temme.
  5. Training-based (TFLO on U = 0) + smoothing: bias not rigorously bounded; include
     as the favorable envelope with the Phasecraft calibration point (below).
- Wall: X_NISQ = p G_lc ≲ ln(m S₁ ε²)/α. With S₁ ≈ B/t_attempt ≈ 6e10 and ε = 0.01,
  ln(S₁ε²) ≈ 15.6 and each factor of 10 in m or B adds 2.3. So the NISQ ceiling is
  G_lc ≲ (4–8)e3 two-qubit gates at p = 1e-3, α = 2–4, nearly independent of budget.
- Calibration point (Phasecraft/Quantinuum, 2510.26300): 7×4 torus, 56 qubits, 2415
  two-qubit gates at 2q error 1e-3 (p G ≈ 2.4, XEB 1–5 %), 160 shots per point,
  pseudo-twirling + spin-reflection averaging + TFLO + GPR, doublon-type observables to
  ≈ 0.05 (read off their Figs. 21–23), no post-selection, 19 % of shots at correct
  Hamming weight, ~5 s per shot on H2. This is an executed point on the NISQ favorable
  envelope at 5× our error target with 60× fewer shots than N₀. Recompile it for the
  nearest-neighbor array (their SWAPs are free relabelings) before comparing.

### 4.2 STAR

- Patches: N_patch ≈ 1.5 (2N + 6) (Toshio floor plan, sequential rotations) or the
  Litinski compact block (Akahoshi); q_arch = N_patch · 2d² plus injection patches
  (2–3 per resource state per clock at p = 1e-3, Toshio Fig. 10).
- Rotation error models (all three from §2), with the RUS chain folded in:
  - Original: P = 4p/15 per rotation, angle-independent.
  - Toshio: P̃_L = α_RUS |θ| p, α_RUS ≈ 1.5 (protocol switching above ~0.1 rad), no
    angle-independent injection floor when the injection distance equals the data
    distance. Residual coherent miscalibration Σφ²/2 not budgeted in the paper; budget it.
  - TopQAD: injection at d_i < d_f then growth; per-rotation error 0.4 k |θ| p (k = ⌈d_i/2⌉)
    + 2 I_growth(d_i, d_f), with I_growth = μ d_i² Λ^{−(d_i+1)/2}. At d ≈ 25 the
    full-distance transversal protocol needs k ≈ 12 physical rotations of angle θ^{1/k}
    ≈ 1 rad and p_ideal collapses, so growth is *forced* at p = 1e-3. The floor is
    therefore not optional at our p; the question is its size.
- Clifford budget: N_patch · N_clock · C · p_L(d) ≤ q_L (C9), C = cycles per clock
  operation; d follows. Toshio Table II at p = 1e-3, PBC 8×8, QPE: d = 25, 2.5e5 qubits.
- Time per rotation: sequential lattice surgery C_av ≈ 4.6 clocks (Toshio) or 18 clocks
  including RUS (Akahoshi), × d × t_cycle (C7) ≈ 25 × 240 ns × 4.6 ≈ 28 μs per rotation
  at d = 25. Parallel rotations (Akahoshi et al. arXiv:2408.14929) as the alternative.
- Wall: X_STAR = 1.5 p Θ_eff + N_rot,eff q_floor ≲ ln(m S₁ ε²)/4. At p = 1e-3 the
  angle term alone gives Θ_eff ≲ 2.6e3 rad, i.e. λ_lc τ ≲ 2.6e3, or N_lc τ ≲ 520.
- Qubitization on STAR: possible (Toffolis become π/4 analog rotations with
  angle-proportional error), just expensive. Evaluate rather than exclude.

### 4.3 Full FT

- Logical qubits 2N + O(1) (Trotter) or 2N + O(log N) (qubitization).
- Synthesis: T per rotation ≈ 1.15 log₂(1/ε_synth) + 9.2 (RUS, used by CUDA-Q Logical)
  or 3.02 log₂(1/ε_synth) + 0.77 (pyLIQTR, Ross-Selinger fit); ε_synth = ε_syn,tot/N_rot.
  Both arms of the FT number must use one law.
- d from V_L p_L(d) ≤ q_L (C8, C9). Factories: 15-to-1 at the same d (one or two levels
  at p = 1e-3), or cultivation (viable ≲ 1e5 T at 1e-9 per T per TopQAD Fig. 16).
  Include factory rejection and RUS delay in throughput, not as circuit rejection.
- Time = max(T-supply-limited, Clifford-limited, reaction-limited) × N₀/m.
- Qubits = 1.5 × (2N) × 2d² + factories. Reference: TopQAD full-FT Hubbard 8×8 at
  ε = 0.01, qubitization: ~1e6 qubits, minutes (target hardware, 2q 5e-4). CUDA-Q
  Logical L = 16 plaquette-Trotter QPE on Pinnacle qLDPC: 5.9e4 qubits, ~300 s at p = 1e-3
  (different code; not directly comparable).

---

## 5. The collapse test (primary deliverable)

For each arm, compute cost over a grid in (L, τ, r, mapping, strategy) and plot:

1. **Γ² (or total shots) against X** for NISQ and STAR, and **d against log(V_L/q_L)**
   for FT, all points from the grid on one axis per arm. Collapse = the hypothesis holds.
   Where it fails, the failure mode is the result:
   - NISQ: acceptance vs. conditional bias under post-selection; light-cone truncation
     error competing with the exponent; idle noise during measurement.
   - STAR: floor term overtaking the angle term (predicted at N_rot,eff ≳ 1/q_floor,
     TopQAD's 1e5–1e6 Goldilocks); Clifford volume forcing d up faster than the angle
     budget shrinks; injection supply rate at p = 1e-3.
   - FT: factory-limited vs. Clifford-limited crossover; reaction-time limit.
2. **Q_min(N, τ, ε = 0.01; T ≤ 1 week)** vs N for the τ families in C14.
3. **Reachable (L, τ) regions at fixed Q ∈ {1e2, 1e3, 1e4, 1e5, 1e6}**, three arms overlaid.
4. **Sensitivity strips** on: p ∈ {1e-4, 1e-3}; suppression law (C8); floor model; mapping
   (C10); PEC α; replicas on/off; light cone on/off.

Retain the qubit/runtime trade-off rather than collapsing to a space-time volume.

---

## 6. Questions only this group can settle

Framed as decisions with proposed defaults, because the literature disagrees on them by
assumption rather than by calculation.

1. **Does angle-proportional STAR error survive at p = 1e-3 under the shared noise
   model (C6)?** Toshio's numerics are at 1e-4 with a 1 μs cycle; TopQAD's growth floors
   are at 2q error 5e-4. Default: adopt the TopQAD floor model with our (A, p_th) and
   measure its size by Stim on the growth step alone (that part *is* Clifford).
2. **Is light-cone PEC a fair envelope for NISQ, or an optimistic one?** It is exact
   for the observable's causal cone but the cone of a second-order Trotter circuit at
   τ = 2 already covers ~(4vτ)² sites. Default: compute both, report both.
3. **One suppression fit for both encoded arms (C8): is that honest?** STAR's lattice
   surgery and FT's are the same operations; the Akahoshi fit differs because of their
   noise model, not their architecture. Default: yes, one fit.
4. **What is the effective α in NISQ?** PEC gives 4; symmetry post-selection with
   residual bias gives an effective 1–2 for the detectable fraction; training-based
   methods are not characterized by an exponent at all. Default: report the menu.
5. **Where does the Trotter step law come from?** Empirical small-L fit with stated
   extrapolation (C12), not commutator bounds. Someone owns the 4×3 exact runs.
6. **Do we charge STAR for coherent miscalibration?** Toshio does not. Default: budget
   Σφ²/2 with |φ| ≈ √p per physical rotation, as a sensitivity.

---

## 7. Workloads

| Purpose | Workload |
|---|---|
| Main scaling study | C1–C3 baseline, L ∈ {4, 6, 8, 10, 12, 16}, τ families C14 |
| Executed reference point | Phasecraft 7×4 torus, π flux, dimer state, U/J = 4, τ ≤ 2, recompiled for the NN array |
| Initial-state sensitivity | 2–3 square-lattice cases with dimer/holon-doublon preparation |
| Algorithm sensitivity (FT, and STAR if cheap enough) | qubitized dynamics via pyLIQTR `FermiHubbardSquare` (queries ≈ 1.4 α τ + ln(1/ε), α = 4N|J| + NU in site count) |

Keep matched-algorithm and optimized-algorithm results separately identifiable: STAR's
small-angle advantage and FT's freedom to change algorithm are both real and both
belong in the paper.

---

## 8. Tools: roles, status, and the adaptations they need

All tool output calibrates prefactors in the inspectable cost model; none of it is the
model. Status checked 2026-09-22 on Improv (CPU login node, no GPU).

| Tool | Status | Role | Needed adaptation |
|---|---|---|---|
| Cost model (`model.py`, to be written) | — | The deliverable | Consumes workload summaries (width, depth, gate/rotation histograms, Θ, Clifford volume) from any source |
| Exact small-L dynamics | `qutip` env or numpy; `qb_gsee` has openfermion | C12 Trotter law, light-cone check, mitigation α via density-matrix noise | **Stim cannot run generic Hubbard rotations**; use a density-matrix/trajectory simulator for ≤ 18 qubits |
| Stim + PyMatching (`stim` env) | works | STAR growth-floor and memory p_L(d) under C6/C7; nothing else | Growth step and Clifford part only |
| CUDA-Q Logical (`~/soft/cuda-quantum/preview/logical`) | source only, **not built**; needs py ≥ 3.11 + `cuda-quantum-cu12` + CMake | FT arm: logical counts, Clifford+T synthesis (`clifford_t_target`), analytical/schedule physical estimates | `fermi_hubbard_jordan_wigner` is **first-order, controlled, periodic**; rewrite for ordinary second-order dynamics. `surface_physical_target` is a compact reference realization with native multi-Pauli measurements, **not our NN array**; supply a device realization or use only its logical tier |
| Base CUDA-Q | not installed | `cudaq.synth.gridsynth` for exact T counts per rotation class; `compiler-bench-ftqc-clifford-t` | pip install; CPU targets suffice |
| pyLIQTR 1.4.2 (`qb_gsee` env) | works | Qubitization workload counts; scheduler for T-depth | No Trotter builder, no physical layer; α convention above |
| Azure QRE (`qsharp` in `qb_gsee`) | installed, not validated | Independent FT physical estimate from logical counts with our (t_g, t_meas, p) | Validate on a known case first |

---

## 9. Sequence

Not time-boxed. Each step produces something the next one needs; the order is logical,
not a schedule.

1. **Conventions memo.** Settle C1–C15 and §6 in writing, one line each, with owners.
   Nothing numerical before this.
2. **Cost model with exposures.** Implement §4 with X_NISQ, X_STAR, V_L as explicit
   intermediate outputs, all reference constants in one table with provenance.
   Produce the first collapse plots from the model alone; label them as model predictions.
3. **Small-L exact calibration.** Trotter law r(τ, ε) for the doublon observable;
   light-cone size vs τ; noisy density-matrix runs under C6 for α, acceptance, bias.
4. **Architecture calibration.** Stim: p_L(d) and growth floor under C6/C7. CUDA-Q
   Logical or QRE: T counts, factory throughput, d, makespan for the FT arm on the
   matched workload. STAR: clocks per rotation and Clifford volume for both floor plans.
5. **Sweeps and collapse test.** §5 figures 1–3.
6. **Sensitivity.** §5 figure 4; the experimental reference point; qubitization.
7. **Write-up** around the collapse test and the crossovers it implies, with the
   conventions table as the methods section.

---

## 10. Reference constants (inputs, with provenance; all to be re-derived under C6–C8)

| Constant | Value | Source |
|---|---|---|
| STAR original rotation error | 2p/15 per attempt, 4p/15 per RUS rotation | Akahoshi 2024 Eqs. 6, 12–15 |
| STAR original PEC overhead | e^{16pN/15} | Akahoshi 2024 Sec. IV.B |
| Small-angle rotation error | P̃_L = α_RUS|θ|p, α_RUS ≈ 1.5 (e^{iθZ} convention) | Toshio 2025 Eq. 25, Sec. IV.D |
| Small-angle PEC overhead | e^{4 α_RUS Θ p}, Θ = λτ for Trotter | Toshio 2025 Eqs. 29–30, 37 |
| Undetectable-error constant | P_ud = (k/15)p native ZZ, (2k/15)p with CNOT–Rz–CNOT | Toshio 2025 p. 10 |
| Growth floor | I_growth = μ d_i² Λ^{−(d_i+1)/2}; 1.7e-5 (d_i = 9) to 4e-8 (d_i = 15) at 2q 5e-4 | Chung 2026 Eq. 5, Fig. 10b |
| Goldilocks window | 1e5–1e6 small-angle rotations | Chung 2026 Sec. IV.D |
| Surface-code fits | 0.1(100p)^{(d+1)/2} (Fowler-Gidney); 0.07(p/0.0039)^{(d+1)/2} (Akahoshi noise) | Toshio p. 27; Akahoshi Table I |
| Clocks per rotation | 4.6 (sequential, Toshio); 18 (compact block + RUS, Akahoshi) | Toshio Sec. VII.C; Akahoshi p. 16 |
| STAR Hubbard 8×8, p = 1e-3, QPE ε = 0.01 | d = 25, 2.5e5 qubits, 3.3e4 h sequential | Toshio Table II |
| Hubbard λ | 5N at U/J = 4 (Z terms dropped) | Toshio App. G |
| T per rotation | 1.15 log₂(1/ε) + 9.2 (RUS); 3.02 log₂(1/ε) + 0.77 (Ross-Selinger fit) | CUDA-Q Logical Eq. 2; pyLIQTR `rotation_gates.py` |
| Full-FT Hubbard 8×8 (qubitization, ε = 0.01) | ~1e6 qubits, 4.7 min time-optimal | Chung 2026 Table V |
| Executed NISQ point | 7×4, 2415 2q gates, p₂ = 1e-3, 160 shots, ≈ 0.05 accuracy | 2510.26300 Tables I–III, Figs. 21–23 |
| Second-order step count, executed | 4 steps to τ = 2, local Trotter error size-independent | 2510.26300 App. C.3 |

---

## 11. References and roles

- `refs/PRXQuantum.5.010337.pdf` Akahoshi et al. 2024 — original STAR; angle-independent error; Hubbard 11N-per-step count (PBC).
- `refs/PhysRevX.15.021057.pdf` Toshio et al. 2025 — small-angle STAR; total-angle bound; Hubbard QPE tables at p = 1e-3, 1e-4.
- `refs/2603.13093v1.pdf` Chung et al. 2026 — STAR with growth vs TopQAD full FT; floors; Goldilocks window. Ground-state task, not dynamics.
- `refs/CUDA-Q Logical ... .pdf` — FT compilation and tiered estimation; Hubbard plaquette-Trotter QPE example.
- `refs/2510.26300v2.pdf` Phasecraft/Quantinuum 2025 — executed 7×4 Hubbard quench; NISQ compilation and training-based mitigation.
- arXiv:2303.06496 Tran, Sharma, Temme — locality-aware mitigation cost (C11).
- Takagi et al., npj QI 8, 114 (2022) — fundamental mitigation limits; use for the exponential-envelope claim's scope.
- arXiv:1912.08854 Childs et al. — Trotter error theory; cite for why the global bound is the wrong tool here.
- arXiv:2408.14929 Akahoshi et al. — parallel STAR rotations (alternative to sequential floor plan).
