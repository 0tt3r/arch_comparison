# Scientific review of the integrated comparison (Opus 5.5, 2026-09-23)

Scope: whether `results/integrated_v1` and the manuscript are *right*, *defensible* and
*best-possible* (within ~2×), not only internally consistent; and whether our methodology,
`fh_iqus_idears/` (IQUS, current model `1206aa66c52e`) and `crossover_martin/` (identical to
`nisq_ftqc_crossover_v4.zip`) are the same or different. Five independent reviews were run
(NISQ, STAR/hybrid, FT, shared task/metric, cross-group methodology). No repository files were
modified; the re-runs quoted below are scratch calculations (monkeypatched copies of our code
and read-only copies of the partner code), not saved results. They are indicative and should be
reproduced in the repository before being quoted.

Baseline state at review time: `make test` 223 passed; `integrated_comparison.py --check` and
`workshop_reconciliation.py --check` pass. The reconciliation check passes only because it
replays `source_snapshot.zip`: most IQUS `fhcost/*.py` hashes no longer match
`results/workshop_reconciliation_v1/comparison.json`.

---

## 0. Corrections after cross-review (read first)

A second reviewer (Codex) checked this document against the code, the saved results, the
Richardson arithmetic, the raw STAR channels and the free-fermion calculation. The corrections
below were verified against `results/integrated_v1/costs.json` and
`integrated_comparison.py`, and they supersede the corresponding statements in Sections 1–8.
Sections 1–8 are left as originally written so the record is intact.

**What stands**
- The fixed-τ = 1 capacity in L should not be read as more useful physics once the local
  observable has converged. The headline should be the time-growing workload or τ_max(Q), with
  fixed-time capacity kept as a benchmark of full-lattice implementations. That remains a valid
  comparison of specified implementations.
- Several controlling choices are fixed rather than optimized: STAR lanes, FT consumption ports
  and factories together, the HWP batch cap, the `d ≥ source_end_distance` clamp at p = 10⁻⁴, and
  STAR lane tiles charged to FT.
- **PCEC (§4.1) is the strongest new result; Codex reproduced it.** Without PCEC the exposure
  falls from 0.685 to 0.452 at 4×4, with a coherent-bias bound of 1.18×10⁻⁵. At 16×16 it falls
  from 12.60 to 8.02, with a bound of 2.38×10⁻⁴. The cancellation depends on the assumed
  success/failure channels and retry probabilities. The fact that we measure an expectation
  value does not justify it on its own. Raw RUS belongs in the optimizer as an option, with its
  coherent residue and retry tails budgeted, and it should be compared against the published
  PCEC construction (Toshio et al., PRX 15, 021057). *Update:* [docs/RAW_RUS_NOTE.md](docs/RAW_RUS_NOTE.md)
  derives the first-order cancellation and its assumptions. It notes that the paper's own
  branch-averaged Appendix E finds raw RUS *slightly worse* than PCEC, so the instrument and the
  paper disagree. Raw RUS remains a budgeted model option pending reconciliation with the authors.
- The zero-failure growth bound should not be presented as a physical rate. Replacing it with
  zero or an extrapolated rate is a different assumption, not a measurement.
- The three groups share one skeleton and answer different questions. The cross-group
  attribution in §6 and the partner issues in §6.2 stand.

**What was wrong or overstated**
1. **FT baseline description (§5.1) is wrong.** At ε = 10⁻³ the optimizer already selects the
   fast one-clock block, not the nine-clock compact block. At L = 8 and 16 the binding
   constraint is supply, not consumption (p = 10⁻³, τ = 1, component scenario):

   | L | ε | Layout | Factories | d | Supply (s) | Consumption (s) | Bottleneck |
   |---|---|---|---|---|---|---|---|
   | 4 | 10⁻² | compact | 1 | 21 | 2.02 | 2.68 | consumption |
   | 16 | 10⁻² | compact | 1 | 25 | 18.7 | 28.1 | consumption |
   | 8 | 10⁻³ | fast | 64 | 27 | 3.68 | 2.51 | supply |
   | 16 | 10⁻³ | fast | 64 | 29 | 15.55 | 11.02 | supply |
   | 32 | 10⁻³ | fast | 128 | 31 | 33.0 | 48.5 | consumption |

   Adding consumption ports alone cannot remove a supply bottleneck.
2. **Additional factories are not cheap (§5.1, and the "~1k qubits" in the chat summary).** With
   one cleanup level, each factory is charged the 1024-qubit source plus 13 protected tiles
   (`tiles = data_tiles + factories*(1 + 12*cleanup)`): 1024 + 13·2d². At the strict 16×16 point
   (d = 29) that is **22,890 qubits per factory**; going from 64 to 128 factories adds about
   **1.46 million qubits**. The copy there is 3.41M qubits (1989 protected tiles = 1157 data +
   64×13). Joint port/factory optimization is the right test, but more ports carry a large area
   cost. **The 10–100× FT figure in §5.1, §7 and §8 is not established.**
3. **Amplitude estimation savings are proposals.** The query-scaling advantage is established
   (e.g. IQAE, arXiv:1912.05559). The coherent A/A† calls, reflections, maximum circuit depth
   and the extra protection they need have not been costed, so the 30–100× figure is not a
   result.
4. **Richardson gain (§2.2, §3.2) is not automatic.** The residuals reproduce: about 1.25×10⁻⁵
   for the 8/16-step pair and 7.84×10⁻⁷ for 16/32. But ‖c‖₁² = 25/9 is not the full cost. With
   equal branch variances, runtime proportional to steps, and optimal shot allocation, the cost
   is (1/3 + (4/3)√2)² = 4.92 times the base-r runtime. So the 8/16 pair costs about 39
   step-equivalents, against 37 for the current calculation. Any gain has to come from
   reallocating the error budget. For NISQ, each branch also pays its own exponential
   mitigation cost. The "2–3× for encoded arms" claim does not follow.
5. **Interacting convergence (§2.1) is unsupported.** The U = 0 finite-size numbers reproduce.
   But a uniform Lieb–Robinson bound does not make the free-fermion truncation error an upper
   bound on U = 4 dynamics. The U = 0 numbers motivate testing 6×6–8×8 patches at U = 4; they do
   not establish L_conv. Interacting bulk convergence stays explicitly unresolved.
6. **"Our NISQ arm is pessimistic, not optimistic" (§3.6) overreaches.** Attenuation and
   post-selection acceptance can motivate lower-overhead scenarios. They do not determine the
   variance and accuracy of the complete estimator. The matched Willow result contains partial
   cancellation between circuit error and experimental discrepancy, so it cannot calibrate a
   universal κ, particularly at strict precision. The data justify a wider κ band, not a
   direction. The inverse-requirement framing remains the sound choice.
7. **"Every arm is pessimistic" and "within 2× once corrected" are not established.** The review
   shows which knobs are not optimized. It does not show the size of the improvement or that the
   corrected comparison is within a factor of two.
8. **STAR delivery wording (§6.3, §8.4).** STAR delivery is a *major unresolved assumption that
   can flip* the STAR/FT ordering. It is not *the one* assumption: injection channels, growth,
   precision and FT synthesis also move that boundary.

**Agreed next pass (supersedes §8)**

This is a bounded optimization round. Its deliverable is *which choices change the winner, and
why*; improvement factors are measured afterwards, not assumed.
1. Optimize raw RUS against PCEC, charging coherent residue and retry tails.
2. Jointly optimize STAR lanes and FT ports, factories and Hamming-weight batches, each with its
   own area charge.
3. Allow separate source and data distances, charging the transfer between them (this replaces
   the `d ≥ source_end_distance` clamp).
4. Add Richardson with branch-specific runtime, variance, mitigation and systematic-error
   accounting.
5. Refresh all arms together, and lead with the time-growing workload. Keep fixed-time capacity
   as an implementation benchmark.

**Status (2026-09-23):** done. See `results/optimization_round_v1/summary.md` and
`docs/OPTIMIZATION_ROUND.md`, which include the leader-change attribution.

Carried over from §8, outside this round: the Stim-calibrated p_L(d) (§5.4; *done*:
[docs/SURFACE_LAW.md](docs/SURFACE_LAW.md) gives Λ = 7.41 [7.20, 7.62] at p = 10⁻³, i.e. +4–6 in d,
not yet propagated through the sweep), refreshing the IQUS
reconciliation and sending the partner issues (§6.2), and restructuring the manuscript once the
refreshed crossovers exist.

---

## 1. Bottom line

- **Arithmetic: correct.** Every operating point checked by hand reproduces (NISQ 4×4 PEC,
  STAR 4×4 60,648 qubits / 0.0636 s, FT 4×4 75,994 qubits / 2.7 s / d = 21). Hamiltonian,
  angle and term bookkeeping have no factor-2 errors. Ours is the most complete of the three
  models.
- **Headline metric: not defensible as physics.** Capacity in L at fixed τ = 1 measures
  lattice area beyond the point where the local observable has converged (L ≈ 6–8).
- **Best-possible: no, in every arm, and the biases point in different directions.** NISQ,
  STAR and FT are each pessimistic by more than 2× in their controlling quantity, so the
  crossover locations are not yet within 2×.
- **Qualitative picture survives:** NISQ saturates logarithmically in Q; STAR saturates at a
  wall set by its angle exponent; FT has a large fixed cost and then scales; at p = 10⁻³ FT
  eventually beats STAR. The *locations* (NISQ→encoded crossover, STAR→FT crossover, the
  p = 10⁻⁴ STAR window) are not yet defensible.
- **Three groups:** the same skeleton answering different questions. With the task aligned,
  NISQ-PEC and FT agree within 1.3–2.3× in sites; STAR spans 36–237 sites, set by the
  parallel-injection assumption.

---

## 2. Shared task, algorithm and figure of merit

### 2.1 L at τ = 1 is past convergence (most serious)
- What we do: `figures/integrated_capacity.png`, `results/integrated_v1/summary.md` report the
  largest L at τ = 1 (up to 48) for the central-site doublon after a Néel quench; the circuit is
  full-lattice (`physics_sweep.py:54-58`, `integrated_comparison.py:198-212`).
- Check: exact U = 0 Néel quench of the central doublon (factorizes at U = 0) against L = 60.
  Interactions do not raise the Lieb–Robinson velocity above the hopping bound, so this is a
  reasonable proxy.

  | τ | Finite-size error |
  |---|---|
  | 1 | L=4: 3.3e-2; L=6: 2.4e-3; L=8: 2.4e-4; L=10: 1e-5 |
  | 2 | L=6: 1.4e-3; L=8: 2.2e-4 |
  | 4 | L=10: 1.5e-2; L=12: 2e-4; L=16: 1.4e-3; L=20: 8e-5 |

  Converged at L ≈ 6 (ε = 0.01) and L ≈ 8 (ε = 0.001) for τ = 1; roughly L_conv ≈ 3τ + 3 at
  ε = 0.01.
- The repository does not state this. `results/frontier_v1/locality_certificate.json` gives a
  loose bound (radius 78, 157×157 patch at τ = 1) that the docs correctly say is not necessary;
  `results/implementation_v1/locality_comparisons.csv` compares nested patches against a finite
  reference one site from the boundary (0.015 for 4×3 vs 4×4), not against the bulk.
- Impact: the FT lead ("L = 32 at 10⁷") and the STAR plateau (L ≈ 10–12) both lie in the region
  that buys no physics at τ = 1. Restricted to L ≈ 6–8, all encoded arms arrive at ~1–3×10⁵
  qubits and are genuinely tied. At p = 10⁻³, ε = 0.01 the NISQ effective line reaches L = 6
  from Q ≈ 10², so at τ = 1 NISQ-vs-encoded reduces entirely to the residual-bias assumption.
  τ = 1 on 8×8 is also likely classically tractable (light-cone / tensor-network); we show no
  classical line (IQUS does). IQUS OPEN_ITEMS O7 raises the same point.
- Fix: headline τ_max(Q, ε) on a patch L ≥ L_conv(τ, ε) + margin, with light-cone truncation
  error charged to the budget (PLAN.md C11(b), C14 family 2). The existing τ ≈ L/4 family is
  close for L ≳ 12. Add a classical reference.

### 2.2 Trotter-error extrapolation missing (largest efficiency lever)
- What we do: second-order formula, r from C/r² within a share of ε; planning C(τ=1) = 0.29
  (plaquette), 0.21 (termwise) from the 3×3/4×3/4×4 envelope; r = 10–13 at ε = 0.01,
  28–41 at ε = 0.001.
- Our own 4×4 data (`results/dynamics_followup_v1/measurements.json`) are cleanly 1/r²:
  plaquette τ = 1 errors 4.138e-3 (r=8), 1.025e-3 (r=16), 2.557e-4 (r=32). Richardson on
  {8,16} leaves 1.3e-5; on {16,32} 8e-7. r→2r ratios 4.04 and 4.05 at τ = 2, 4.
- Impact: at ε = 0.001, r ≤ 16 instead of ~37 at variance cost ‖c‖₁² ≈ 2.8. Encoded arms
  ~2–3× cheaper in time; NISQ and per-rotation STAR exponents drop >2× (exponential in cost).
  IQUS independently finds multiproduct formulas are their largest lever.
- Higher-order Suzuki is the wrong option for STAR (signed stages raise Θ ≈ 2.3×); MPF /
  Richardson preserves Θ per circuit. Refs: Low–Kliuchnikov–Wiebe 2019; Carrera Vazquez et al.
  2023; Endo et al. 2019; Watson 2024.

### 2.3 Transfer of the Trotter coefficient to larger L
- Size-independence of C for a local observable at fixed τ is correct in principle (Childs et
  al., PRX 2021). Patch spread at τ = 1 is ~1.4× in C (≤1.2× in r): fine for encoded arms,
  1.2× in the NISQ exponent.
- C grows ~τ^2.2–2.5 (0.29, 1.25, 5.8 at τ = 1, 2, 4). A τ = cL headline needs larger-patch or
  free-fermion-anchored checks; IQUS O9 sees ~3× spread at τ = 2 across sizes.

### 2.4 Shot count and observable
- N₀ = 3.84 × 0.09 / 0.004² ≈ 21,600 (Gaussian, D* = 0.1). The exact 4×4 D(τ=1) = 0.145
  (variance 0.124), so shots are ~1.4× low for every arm; ordering unaffected.
- The PEC reference line uses Hoeffding while other curves use Gaussian planning
  (`discussion_model.py:26-27,103`): ~2–5× more shots, small in χ, but an inconsistency to
  state or remove.
- Single-site readout: averaging over bulk sites of a converged patch (as Phasecraft does) cuts
  shots by roughly the number of weakly correlated interior sites for every arm.
- Error budget: linear addition of biases, 95% on the statistical part only — consistent across
  arms and conservative.

### 2.5 Hamiltonian bookkeeping — correct
Hopping (J/2)(XX+YY) per spin per bond → 2J per bond, 4N−4L total; U n↑n↓ → (U/4)ZZ (1 per site
at U = 4); single-Z terms legitimately dropped (number conservation). λ = 5N−4L; 11N−8L terms
(9N−8L without single Z). Θ = τλ = 64 for 4×4 at τ = 1 in the e^{−iθP} convention, matching the
reconciliation ledger.

### 2.6 Replicas and one-week budget — fair
m = ⌊Q/Q_copy⌋ for every arm under the same deadline; NISQ's gain from ~10⁷ copies is
correctly logarithmic.

---

## 3. NISQ arm

Provenance: `plot_integrated.py:42-45` replays `results/physics_sweep_v1/costs.json`
(`discussion_model.py:93-116`). G = 4372·(r/3)·(L/6)³, depth = 294·(r/3)·(L/6), anchored to
Willow 6×6, 3 steps (`configs/discussion_model.json:12-18`). Dotted PEC line:
χ = G·2 ln(1+2p/s) ≈ 4pG, full circuit, two-qubit noise, Hoeffding. Dashed effective line:
χ = ln(100)(p/10⁻³)G/4372, i.e. κ ≡ χ/(pG) = 1.05.

1. **Effective-mitigation κ is a chosen constant; evidence points to 3–6× smaller.**
   The overhead 100 is chosen (config line 18); the κ = 0.53/1.58 scenarios (lines 35-37) are
   not drawn. Phasecraft/Quantinuum measured observable attenuation Λ = 0.20 at pG = 2.4 (IQUS
   O5) → rescaling κ ≈ 2Λ/(pG) ≈ 0.17. Willow 6×6 strict number-selection acceptance 0.02 at
   pG ≈ 13 (hardware p ≈ 3×10⁻³) → κ_post ≈ 0.3. Damping acts only on the support of the
   Heisenberg-evolved operator (IQUS weight test ZZ/Z = 1.91), so full-circuit fidelity
   overcharges. Re-run at Q = 10⁶, p = 10⁻³:

   | κ | ε = 0.01 (r = 8) | ε = 0.001 (r = 21) |
   |---:|---:|---:|
   | 1.05 (ours) | L = 6 | L = 4 |
   | 0.4 | L = 10 | L = 6 |
   | 0.17 | L = 12 | L = 8 |

   At κ = 0.17–0.4 the NISQ plateau meets the STAR/hybrid plateau at p = 10⁻³; the encoded
   crossover moves from ~10⁵ to ~10⁶–10⁷ qubits. Caveat: no evidence controls rescaling/TFLO
   bias at 0.01 for U = 4; our matched Willow reference comes close (4×4, 3 steps, RMS 0.0103,
   71% within 0.01; `results/willow_reference_v1/summary.md`). Verdict: defensible as a
   scenario, not a best estimate; the largest single uncertainty in NISQ-vs-encoded.
2. **No step extrapolation.** r is optimized (8 at ε = 0.01, 19–21 at 0.001), consistent with
   Willow's measured error. MPF-4 (‖c‖₁² ≈ 2.8) takes r ≈ 8 → 3 at equal κ: exponent ~2.7×
   smaller, L 6 → ~8–10. Compounds with item 1. IQUS O6: the model-vs-experiment gap is
   entirely the step count.
3. **PEC arithmetic correct; the dotted line is incomplete relative to our noise model.**
   2 ln γ = 4.0003p per depolarizing 2q gate at p = 10⁻³ (IQUS 4.000268p; matches van den Berg
   et al. 2023 first order). The dotted line omits 1q, idle and SPAM that C6 specifies; the
   earlier 4×4 decomposition was χ = 10.6 (2q) + 4.1 (1q) + 5.2 (idle) + 0.14 (SPAM), so a
   consistent compiled PEC is ~2× the plotted value; the idle convention (p per 10 ns layer) is
   in turn pessimistic for dynamically decoupled hardware. Exact circuit-cone PEC gains little
   with a swap network (modes traverse the array each step); IQUS's geometric cone fraction is
   optimistic for this compilation. Net ≤ 2× in χ, ≤ one even-L step
   (`model.py:86-93`, `discussion_model.py:102-103`).
4. **Compilation within 2× at plateau sizes.** Willow compiler: 400r+20 (4×4), 1440r+52 (6×6)
   CZ, ratio consistent with L³. Compact local encodings (Derby–Klassen; Cade et al. 2020) give
   O(N) gates per step at ~1.5N qubits; swap networks O(N^{3/2}). At L = 6–8 the gap is ≤1.33×.
   L_max ∝ (κr)^{−1/3} makes the plateau insensitive to constants; with G ∝ L² the plateau
   moves 6 → 8 (κ = 1.05) and 12 → 20 (κ = 0.17). Show compact encoding as the best-known
   alternative (`integrated_comparison.py:319-324`, `discussion_model.py:97-98`).
5. **Hand check (4×4, PEC, ε = 0.01, r = 8):** G = 3455, χ = 13.82, γ² = 1.0×10⁶; depth 523
   layers → 5.43 µs; 9.4×10¹⁰ shots (N₀ ≈ 9×10⁴ Hoeffding) → 5.1×10⁵ s on one copy → Q_min = 36,
   matching log₁₀ Q_min = 1.556. Correct.
6. **Experiments.** Willow 4×4 (1220 CZ, p ≈ 3× ours, 2×10⁵ shots/time) is within ~0.01 of
   continuous dynamics; our effective model gives χ = 1.3 at p = 10⁻³ for it. Willow 6×6 U = 4
   misses 0.01 by ~5× against an uncertain TDVP reference. Phasecraft reached ~0.05 residuals at
   U = 0 with 160 shots. Neither demonstrates our ε contract at scale; neither contradicts
   κ ≈ 0.2–0.4 for variance. Our NISQ arm is pessimistic, not optimistic.
7. **Other groups.** Martin: F = exp(−ΣpG), variance 1/F² (κ = 2, full circuit, no bias), wrong
   quantity for a local observable and omits bias. IQUS: PEC κ = 4 on cone gates by default,
   though their own O5 measurement shows ~13× overcharge.

NISQ verdict: arithmetic right; cost enters through ln(budget), so a factor 2 in cost is
irrelevant, but the defensible κ range spans ~25× (0.17–4) and MPF adds ~2.5×. The plateau at
Q = 10⁶, p = 10⁻³ could be L = 4–12 (central 6) at ε = 0.01 and 2–8 at ε = 0.001. "NISQ walls
below the encoded arms at p = 10⁻³" is not robust; "FT eventually wins (L ≥ 16 at Q ≥ 10⁷)" is
unaffected. Best practice: anchor κ to the attenuation measured in data we already hold, add
MPF, draw NISQ as a κ band with its bias requirement (the requirement framing in
`nisq_requirements.png` is already right).

---

## 4. STAR and hybrid arms

1. **Noisy PCEC raises χ ~1.55× for no benefit to an expectation value.**
   `mode='noisy_pcec'` is hard-coded (`integrated_comparison.py:284`); ideal PCEC composes the
   channel with its conjugate, doubling the dephasing (`star_channel_study.py:77-92`). Toshio et
   al. (PRX 15, 021057, Sec. IV C) use PCEC to reduce worst-case (diamond) error per trial. For
   an expectation value the RUS-averaged channel already cancels first-order coherent error for
   any over-rotation sequence δ_m (success probability 1/2, failure conjugates):
   Σ_m δ_m[2^{−m−1} − Σ_{n>m}2^{−n−1}] = 0. Our raw-RUS channel confirms: θ = 0.01, d = 9 gives
   Im c = 3.5e-9 vs dephasing 1e-5; the workload coherent-bias bound is 1.2e-5 (4×4) and 2.4e-4
   (16×16), ≤2.4% of ε = 0.01. Without PCEC χ(4×4) 0.685 → 0.452, χ(16×16) 12.6 → 8.0,
   effective α 2.7 → 1.7 (Toshio 1.5). At p = 10⁻³, ε = 0.01: log₁₀ Q_min L=12 6.91 → 6.07,
   L=16 10.08 → 8.04, L=20 13.81 → 10.60. Same for the hybrid's analog rotations.
   **This first-order cancellation argument should be confirmed by the STAR authors.**
2. **Central growth error is a shot-count bound, not a rate.** The component scenario uses
   `pauli_total_upper` from `results/frontier_v1/growth_summary.json` = 1.2347534e-05 for every
   zero-failure campaign (9→13, 9→25 at both p; 5→25 at 10⁻⁴). At p = 10⁻³ the physical
   estimate A(p/p_th)^{(d+1)/2}×rounds is ~1e-5, so it is plausible; at p = 10⁻⁴ it is ~1e-10.
   At p = 10⁻³ growth is about half of α (1.77 with vs ~0.95 injection-only). At p = 10⁻⁴ it is
   ~6× the injection error, and pushes the optimizer to seed 25 at acceptance 5e-3
   (38–157 s/attempt). The central p = 10⁻⁴ STAR curve is an artifact; seed-9 zero-growth is the
   realistic line:

   | p = 10⁻⁴, ε = 0.01 | Q = 10⁶ | 10⁷ | 10⁸ |
   |---|---|---|---|
   | STAR central (L) | 16 | 20 | 20 |
   | STAR zero-growth (L) | 24 | 32 | 40 |
   | FT (L) | 24 | 48 | 48 |

   (`integrated_comparison.py:235-243`, `configs/integrated_comparison.json` scenarios.) Use the
   bound as the pessimistic band edge and a physical estimate as the centre.
3. **Four analog lanes fixed, not optimized.** `parallel_lanes = 4`, 6 support tiles per lane
   (`configs/discussion_model.json:22`); attempt time ∝ N·r/4. Spatially separate rotations can
   run concurrently (Toshio Table I; Akahoshi et al. arXiv:2408.14929). With item 1 at 64 lanes
   (p = 10⁻³, ε = 0.01): L = 16 at 10^6.86 instead of 10^8.04; L = 20 at 10^9.38 instead of
   10^10.6. Items 1+3 move STAR capacity 12 → 16 at 10⁷ and 12 → ~20 at 10⁹. Lanes should be an
   optimization variable with area charged per concurrent injection (2–3 patches).
4. **Hand check (4×4, p = 10⁻³, ε = 0.01): correct.** 2221 rotations at r = 10;
   Θ = 1920×0.025 + 282×0.05 + 19×0.1 = 64. Trial = 4/acceptance + (d_data+1) + 9d_data rounds at
   240 ns ≈ 36 µs; ~2 RUS trials → 72 µs/rotation; 2221 × 72 µs × (19/13)/4 + 5 ms handoff =
   63.6 ms (stored 0.0636 s). Tiles 1.5(32+6)+3 = 60 + 24 lane = 84; 84 × 2 × 19² = 60,648.
   Clifford budget 84 × 2.65e5 × 0.1 × 0.1^10 = 2.2e-4 ≤ ε × 0.05 × e^{−χ/2} = 3.5e-4. ~4.5 clocks
   per rotation vs Toshio sequential ~9.2 (C_av = 4.6 per trial).
5. **Rotation-error law: correct and the most complete of the three.** Exact enumeration of the
   transversal multirotation construction (undetectable generator faults p/15 per pair factor,
   p/30 final factor, P_ud ≈ kp/15, Toshio Eq. 22), workload angles, acceptance, RUS, wrapping mod
   π/4, sampling multiplier Π(1−2q)^{−2}; logical-fault bias amplified by e^{χ/2} (conservative,
   C9). α is angle- and seed-dependent (~0.8–2 at d = 9, ~2–2.4 at d = 25, injection only).
   Limitation: other syndrome/preparation faults enter only through acceptance (as in Toshio);
   not validated end-to-end.
6. **The STAR wall.** χ ≈ 4αpτλ ≈ 0.049L² at α ≈ 2.7 (L=12: 6.96; L=16: 12.6); d rises through
   e^{χ/2}. L_wall ≈ √[ln(B·m/(N₀·t_att))/(20αpτ)]. Best-estimate α ≈ 1.7 (no PCEC) or ~1.0
   (negligible growth) puts the wall at 1.5–2.7× more sites (1.3–1.7× in L). The size grid
   (…12, 16, 20, 24; 1.56–1.78× in N per step) is itself near 2× resolution.
7. **Hybrid: well defined and fairly priced** (plaquette circuit; fixed T rotations from the FT
   cultivation source at 2e-9 with e^{χ/2} amplification; dynamic rotations analog). Inherits
   items 1–3.
8. **Other groups.** Martin χ = 0.4pR (per rotation, angle-independent, 2d cycles/rotation): χ
   scales with r not Θ; 4×4 gives 0.89 (ε = 0.01, close by coincidence) and ~2.5 at ε = 0.001 vs
   our 1.12. IQUS χ = 1.066pR is the original Akahoshi 2024 4p/15 protocol, superseded by STAR
   ver.3 small-angle injection; 2.37 at 4×4, 3.5× ours and ~9× the best-estimate (~0.25).

STAR verdict: mechanics correct; systematically pessimistic (PCEC ~1.5×; growth bound up to
~1.8× at 10⁻³ and a controlling artifact at 10⁻⁴; fixed lanes). At p = 10⁻³ the wall is within
2× in L, not in N (≈12 now vs ≈16–20 best estimate). STAR still saturates and FT still overtakes
it at p = 10⁻³, about half a decade later (~10⁶ → ~3×10⁶ qubits). At p = 10⁻⁴ STAR's leading
region widens (but see FT item 2).

---

## 5. Full FT arm

Provenance: `integrated_comparison.optimize_ft` (line 185) → `physics_sweep.solve` →
`budget_optimization.optimize` → `discussion_model.estimate` (lines 147–188). At ε = 0.01 every
FT point needs one copy (footprint-bound: 2·tiles·d² + factories; runtime enters only via d).
At ε = 10⁻³, L ≥ 6 needs 4–118 copies (time-bound; runtime converts directly to qubits).

1. **ε = 10⁻³: 10–100× pessimistic in Q for L ≥ 8.** Serial T consumption through one port
   (`T_consumption_ports: 1`, `configs/integrated_comparison.json:19-20`), compact block at
   9 clocks/T, and standard sampling N₀ = z²σ²/ε_stat² ≈ 10⁶ shots
   (`discussion_model.py:25-31,160,165-169`); per-shot times 0.5–105 s (L = 16: 3.4M qubits/copy,
   15 copies). Rotations in a Trotter layer commute and are disjoint, so T consumption
   parallelizes with a matching factory bank (Litinski 2019; Beverland et al. 2022). The doublon
   is a projector, suitable for iterative/robust amplitude estimation (~c/ε calls, c ≈ 3–10).
   IQUS and Martin both parallelize (13–56 ms vs our 590 ms per shot at the reconciliation 4×4
   point; 6 ms vs 2.7 s under IQUS's ARCH_COMPARISON preset). Parallel ports ≥10× in runtime for
   a few ~1k-qubit factories; amplitude estimation ~10⁴ coherent calls instead of 1.2×10⁶ shots
   (30–100×) for +4–6 in d (~1.5× footprint). Combined: ε = 10⁻³ points at L ≥ 8 collapse to
   1–2 copies; L = 16 ~5×10⁷ → ~5×10⁶; L = 32 5.4×10⁸ → ~2×10⁷. The STAR–FT crossover at
   p = 10⁻³, ε = 10⁻³ moves from ~5×10⁶ to ≲10⁶. Amplitude estimation is incompatible with PEC,
   so it is a genuine FT-only lever. The NISQ-requirements comparator uses these FT runtimes and
   is correspondingly too lenient to NISQ.
2. **p = 10⁻⁴: distance clamped to d ≥ 15 (~2× against FT).**
   `d = max(d, source_end_distance)` with 15 (`discussion_model.py:172`;
   `integrated_comparison.py:159` sets `minimum_d = max(target, source_end_distance)`). The
   cultivation output distance does not set the data distance. Re-run unclamped:

   | p = 10⁻⁴, ε = 0.01 | d now → unclamped | Q/copy now → unclamped |
   |---|---|---|
   | L = 4 | 15 → 9 | 39.3k → 17.8k |
   | L = 8 | 15 → 11 | 104k → 56k |
   | L = 12 | 15 → 11 | 212k → 116k |

   STAR at L = 8 needs 77k (itself clamped at d = 13 by its growth target). The summary's
   "STAR leads at p = 10⁻⁴, Q = 10⁵, ε = 0.01" is very likely an artifact. Never binds at 10⁻³.
3. **Hamming-weight phasing capped at batch 16.** `hwp_batches: [1,4,8,16]`; Hubbard angles are
   uniform across sites/plaquettes, so whole-group phasing is available (Gidney 2018; Kivlichan
   et al. 2020; Campbell 2022 Thm 2, used by IQUS). Our `synthesis_plan`: L = 32, ε = 0.01:
   1.51M T (16) → 0.72M (64) → 0.36M (1024); ε = 10⁻³: 6.36M → 1.30M (1024). Workspace ~N logical
   (+25–50% tiles). Negligible where footprint-bound (ε = 0.01), compounds with item 1 at 10⁻³.
   Gridsynth T = 3.02 log₂(1/η) + 0.77 (69–102 T/rotation, `resource_extensions.py:76-77`);
   mixed-fallback/probabilistic synthesis (Kliuchnikov et al. 2023; Hastings/Campbell) or a
   phase-gradient register + adder gives ~2–3× fewer T per rotation.
4. **p_L law uncalibrated but acceptable (±2 in d).** A = 0.1, p_th = 0.01 (Λ = 10 at 10⁻³;
   `configs/assumptions.json:20-21`), shared with STAR so it cancels in STAR–FT. C8 Stim
   calibration was never done. Quick Stim+PyMatching memory at p = 10⁻³: SD6-style d = 3/5/7 →
   2.6e-4, 2.9e-5, 2.7e-6 per round (Λ ≈ 9–11, A ≈ 0.03); our saved 240 ns schedule circuits
   (d = 3, 5 no-growth controls) ~1.5× worse (Λ ≈ 7–9, A ≈ 0.04–0.05). Extrapolated to
   d ≈ 21–27: correct to ~5× optimistic in p_L → 0–2 in d, ≤20% in qubits for both encoded arms.
   Cheap to replace with the measured fit.
5. **Hand check (4×4, p = 10⁻³, ε = 0.01): correct, with small overheads.** 832 synthesized
   rotations + 1728 fixed T (plaquette) at 69 T/rotation → 59,136 T. Tiles 1.5·38+3 = 60 + 24
   lane + 1 factory = 85. 59,136 × 9 × 21 × 240 ns = 2.68 s (1.1×10⁷ cycles);
   0.1·10^{−(d+1)/2} ≤ 10⁻³/(85×1.1×10⁷) → d = 21; 2·85·21² + 1024 = 75,994. HWP (22.8k T)
   correctly rejected (workspace > time saved). The 24 STAR lane tiles are charged to FT
   (`discussion_model.py:76`: `+ extra_tiles_per_lane*parallel_lanes` in the shared tile
   function; +28% at L = 4, negligible at L ≥ 16), and the serial runtime pushes d from ~17–19 to
   21 (~1.2–1.5×).
6. **Magic-state supply reasonable.** Cultivation at 2e-9, 1024 qubits, ~100 cycles/T (~10⁵
   qubit-cycles), optional 15-to-1 cleanup — within ~1.5× of Gidney–Shutty–Jones 2024 volume.
   Source error held at the 10⁻³ value at p = 10⁻⁴ (pessimistic, not binding). Factories never
   bind at ε = 0.01, so the source menu does not move the headline curves.
7. **Other groups.** Per-copy footprints at the shared 4×4 point: Martin 62k, ours 79–89k, IQUS
   126k — within 2×. Time per shot: ours 10–50× slower (serial T). IQUS "FT = 0 at ε = 0.01" is a
   task difference (ε relative to a decaying signal, t_max = √m/v, commutator-bound r by
   default, 1 µs cycle; ε_abs = 3.7×10⁻³ at m = 16 and shrinking; ~10⁷–10⁸ shots). Under their
   ARCH_COMPARISON preset their FT 4×4 is d = 19, 78k qubits, 6 ms/shot — footprint matches ours.

FT verdict: ε = 0.01, p = 10⁻³ correct and within ~1.2–1.6× of best (biased against FT); the
p = 10⁻³ STAR–FT crossover near 10⁶ is defensible on FT grounds. p = 10⁻⁴ ~2× against FT
(clamp); the STAR window there is probably not real. ε = 10⁻³ is 1–2 orders of magnitude in Q
from best-possible for L ≥ 8. The systematic bias is against FT throughout.

Note the tension between STAR item 2 (zero-growth STAR widens its p = 10⁻⁴ window) and FT
item 2 (unclamped FT closes it). Both corrections must be applied together before any
p = 10⁻⁴ ordering is stated.

---

## 6. Methodology comparison of the three efforts

### 6.1 Controlling choices

| Axis | Ours (integrated_v1) | IQUS (current) | Martin |
|---|---|---|---|
| Observable, initial state | Central doublon, Néel | NN Czz, triplet-dimer | ⟨Z⟩ from Néel (generic ±1) |
| Time | Fixed τ = 1 (and 2) | t = √m/v_B, v_B = 2, plus 20 paid time points (`budget.py:64`, `ftqc.py:330`) | Fixed tJ = 4 |
| Accuracy | Absolute ε = 0.01/0.001; variance D(1−D) = 0.09 | 10% of smallest signal on the grid, floor 0.02 → ε_abs ≈ 0.0037 (`hubbard.py:441`); variance 1 | Absolute ε = 0.1; variance 1; z = 2 |
| Trotter | Measured 3×3/4×3/4×4 envelope, size-independent; C(1) ≈ 0.21–0.29 | Measured W_eff(τ), size-independent; 0.19 at τ = 1 (`hubbard.py:39`) | Guessed extensive 0.1·N·t³/r² (`benchmarks/hubbard.py:37`) |
| NISQ | PEC reference + inverse requirement (variance and bias); no boundary claimed | PEC exp(4pG) over light cone; ZNE and experiment-calibrated variants | 1/F², F = exp(−pG) global depolarizing (`architectures.py:92`) |
| STAR error | Full angle-dependent instrument: acceptance, RUS, measured growth, coherent bound | Count law 1.066pR (default) or 4αpΘ + 4fR | 0.2p per rotation (χ = 0.4pR) |
| STAR timing | 18 clocks per rotation, 4 lanes (~26 µs/rotation) | 20·d rounds per Trotter step regardless of rotation count; no injection area (`ftqc.py:445`) | 2d cycles per rotation layer (17 layers/step); no injection area |
| FT | 3 log₂ gridsynth or HWP; cultivation + cleanup; one consumption port | HWP over the whole group; factories sized for throughput; 2× bias per failure | 3 log₂ gridsynth; Litinski table; parallel rotation layers |
| p_L | 0.1(p/0.01)^((d+1)/2) | same | same |
| Syndrome cycle | 240 ns | 1 µs (240 ns option) | 1 µs; NISQ 40 ns gates, 1 µs readout |
| Replicas | ⌊Q/area⌋ | ⌊n/area⌋ | ⌊n/area⌋ |
| Size axis | Side L of fully evolved lattice | Sites m (lattice grows with time) | Sites L² |

### 6.2 Classification of differences

**Structural (different questions)**
- Time convention dominates. At fixed τ a local observable converges in L once L exceeds
  ~2v_Bτ plus a buffer (~6–8 at τ = 1). Our figure and Martin's measure how large a lattice
  each arm can fully evolve, not the largest converged answer. IQUS's cone switch shows the
  pathology directly: fixed τ with a light-cone discount lets NISQ "reach" 333,333 sites
  (`presets.py` bisection). Tying τ to L (IQUS) or scanning τ at converged L avoids it.
- Accuracy contract: IQUS requires 20 time points at 10% of a signal decaying to ~0.04 (~3×
  stricter than ε = 0.01, ×20 in sampling). Martin's ε = 0.1 is 10× looser than ours.
- NISQ framing: we solve an inverse requirement; IQUS and Martin assert a model.

**Constants and assumptions**
- Trotter: ours and IQUS agree within ~1.5× at τ = 1 (different observables; both find size
  independence).
- STAR delivery is the biggest open constant. At the harmonized IQUS point STAR takes 0.8 ms per
  shot for ~39k rotations (m = 237); ours is ~0.4–0.6 s per shot at L = 10–12 — ~10³× apart.
  The allowed exponent grows as ln(Q·B/(q·t_shot·N₀)), so 10³× in t_shot shifts it by ~6.9,
  more than our whole χ = 4.8 at L = 10. The truth lies between our 4 serial lanes and free
  parallelism (rotations in a layer commute, but each concurrent injection costs ~2–3 patches).
- FT delivery: our serial port gives 9.6 s/shot at L = 12 vs 17.6 ms in IQUS. FT capacity at
  10⁶ is area-bound in both, so capacity barely moves, but our NISQ-requirements comparator
  (encoded runtime) is made too easy.
- Clocks: IQUS and Martin default to off-specification clocks.

**Errors**
- Martin:
  - Extensive Trotter law r ∝ √N for a local observable is wrong in kind; ~6× step overcharge
    at L = 10.
  - FT logical-error charge covers logic cycles only and omits factory-wait exposure
    (`architectures.py:173` vs `:193`).
  - STAR 0.2p is 2.7× below the canonical 16/15·pR exponent.
  - 1/F² from global depolarizing is not PEC (2pG vs 4pG).
- IQUS:
  - `ARCH_COMPARISON` preset labels our headline as "1% relative, eps_abs = 1e-3"; ours is 0.01
    absolute. Its s_sig = 0.1 is ignored because signal_regime stays "curve", giving
    eps_abs = 7.2e-4 (`hubbard.py:354`); the "alt_model" row does not reproduce our model.
  - Doublon costed with variance 1 instead of 0.09 (~11× too many shots).
  - STAR time per step independent of rotation count; no injection area.
  - STAR logical-error criterion lacks the ×2 bias factor applied to FT (`ftqc.py:446` vs
    `:387`).
  - Docstrings say 1.15 log₂ and exponent d/2; code uses 3.0 and (d+1)/2 (docs only).
  - STAR default is the superseded 2024 count law.
- Ours: no arithmetic error found in the compared parts; issues are the fixed-τ size axis and
  conservative delivery (STAR lanes, FT port), plus the items in Sections 3–5.

### 6.3 Attribution of the headline gaps

IQUS at n = 10⁶, cumulative changes with their code (sites; NISQ / STAR / FT):

| Step | NISQ | STAR | FT | Main effect |
|---|---|---|---|---|
| IQUS default | 13 | 22 | 12.5 | — |
| 20 time points → 1 | 14 | 29 | 32 | STAR ×1.3, FT ×2.6 |
| Doublon, fixed signal, 0.01 absolute | 20 | 42 | 44 | ×1.5, ×1.4 |
| τ fixed at 1 | 333k (cone artifact) | 432 | 186 | ×10, ×4.2 |
| STAR α = 1.5 angle law with floor | — | 237 | — | ×0.55 |
| Routing G ∝ L³, 1q/idle/SPAM noise | 10 | 237 | 186 | — |

Against ours (STAR 100, FT 144, NISQ-PEC 16 sites): FT within 1.3×; STAR's remaining 2.4× is
free parallel injection and variance 1 against our full instrument and serial lanes; the
240 ns round has no effect after the τ change (STAR sampling-limited). The disagreement is mostly
τ ∝ √m (~4–10×), then 20 time points (1.3–2.6×) and the accuracy definition (~1.5×).

Martin at 10⁶, cumulative (sites):

| Step | NISQ | STAR | FT |
|---|---|---|---|
| Default | 4 | 16 | 324 |
| ε = 0.01 absolute | 0 | 4 | 25 |
| τ = 1 | 9 | 25 | 81 |
| Measured Trotter coefficient | 9 | 81 | 256 |
| Specification clocks | 25 | 100 | 324 |
| Canonical STAR law | 25 | 36 | 324 |
| PEC one-norm | 16 | 36 | 324 |

Harmonized NISQ (16) matches ours; FT 324 is 2.25× ours (parallel rotation layers, no HWP
workspace or layout overhead); STAR 36 under the canonical law, 100 under 0.2p (= ours).
Martin's "encoded beats NISQ at 3.5×10⁴" at ε = 0.1 moves to ~1.3–1.8×10⁴ harmonized, loosely
consistent with our ~10⁵ given our stricter NISQ exponent.

Net: within 2× once aligned for NISQ-PEC and FT (1.3–2.3× in sites). Not for STAR: 36 → 100 →
237 sites (Martin → ours → IQUS), and the STAR/FT ordering at 10⁶ can flip on the STAR delivery
assumption (IQUS STAR 237 > FT 186; ours and Martin FT ahead).

### 6.4 Soundest choice per axis

| Axis | Soundest | Adopt |
|---|---|---|
| Figure of merit | IQUS in spirit | Couple L to τ, or τ_max(Q) at L ≥ L_conv(τ, ε) with measured finite-size error |
| Accuracy | Ours | Absolute ε with the observable's actual variance; one or an explicitly declared number of time points |
| Trotter | Ours and IQUS (they agree) | Measure along (m, τ = √m/v) as IQUS O9 proposes; drop the extensive law |
| NISQ | Ours | Inverse-requirement framing, exact PEC reference, IQUS support-damping measurement as empirical evidence; 1/F² only as an optimistic bound |
| STAR | Ours for the channel | Keep the angle-dependent instrument; optimize concurrent injections against 2–3 patches of area each; serial and area-charged parallel as band edges; free parallelism is an upper bound only |
| FT | Ours for synthesis/workspace, IQUS for throughput | Keep synthesis selection and workspace; size ports and factories for throughput; include factory waiting in the logical budget; ×2 bias factor consistently |
| Clocks | Specification | 10 ns gates, 100 ns measurement, 240 ns cycle for all |

---

## 7. What survives, and what does not

Robust: NISQ saturates logarithmically in Q; STAR saturates at a wall set by its angle
exponent; FT pays a large fixed overhead and then scales; FT eventually overtakes STAR at
p = 10⁻³; ε = 0.001 is qualitatively harder for sampling-based arms than for FT.

Not yet defensible: capacities beyond L ≈ 8 at τ = 1 as physics; the NISQ→encoded crossover
(~10⁵ vs 10⁶–10⁷ depending on κ and MPF); the STAR→FT crossover (~10⁶ vs ~3×10⁶, and ≲10⁶ at
ε = 10⁻³ once FT is best-practice); any p = 10⁻⁴ STAR window; the NISQ-requirements comparator
(uses non-optimal FT runtimes).

---

## 8. Recommended work, in order

1. **Headline metric.** τ_max(Q, ε) on a converged patch, truncation error budgeted, classical
   reference line; state L_conv(τ, ε) explicitly (the U = 0 free-fermion check is a cheap
   anchor; confirm with interacting ED on the largest feasible patches).
2. **Best-practice options per arm, all through the same optimizer:**
   - Shared: Richardson/MPF Trotter extrapolation; bulk-site averaging; consistent Gaussian
     shot planning for the PEC line; actual doublon variance.
   - NISQ: κ band anchored to Willow/Phasecraft attenuation data; MPF; compact-encoding
     alternative; complete PEC (1q/idle/SPAM) reference.
   - STAR/hybrid: drop PCEC for expectation values (after expert confirmation); physical growth
     rate at the centre, shot bound as band edge; optimize lanes with area per injection.
   - FT: parallel T ports with matching factories; amplitude estimation; whole-group HWP;
     cheaper synthesis option; remove the d ≥ source_end_distance clamp; stop charging STAR lane
     tiles to FT.
3. **Stim-calibrated p_L(d)** under the 240 ns schedule for both encoded arms (C8).
4. **STAR parallel-injection area cost** agreed with the STAR authors: a major unresolved quantity that can
   decide STAR vs FT across the three groups.
5. **Refresh the workshop reconciliation** against current IQUS code and send partners their
   specific issues (Section 6.2).
6. Then rewrite the abstract, figure titles and a conclusions section around the surviving laws
   and the corrected crossover locations; restructure the manuscript from its current
   chronological layout into one narrative.
