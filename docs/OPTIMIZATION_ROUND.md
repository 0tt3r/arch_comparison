# Optimization round: methods

This round implements the five changes agreed in [OPUS55_REVIEW.md](../OPUS55_REVIEW.md) §0.
Its deliverable is **which choices change the capacity leader, and why**; improvement factors
are read off afterwards rather than assumed. Results: [summary](../results/optimization_round_v1/summary.md).

## Round v2: calibrated surface-code law

`--round v2` (`configs/optimization_round_v2.json`, `results/optimization_round_v2/`) repeats
the round with the per-p Stim fit of [SURFACE_LAW.md](SURFACE_LAW.md) for every encoded arm
(A = A_X + A_Z, p_th = pΛ; Λ = 7.41 at p = 10⁻³, 89 at 10⁻⁴). Physical-qubit rows must equal
round v1 exactly (1120 checks). Distances rise by 2 at small lattices and up to 6 at the largest;
transitions move about 1.3–1.4× later at p = 10⁻³ and barely at 10⁻⁴. No ordering changes. The
manuscript uses v2. `qsub -v ROUND=v2 pbs/optimization_round.pbs` reproduces it;
`make optimization-round-v2-check` audits it.

## Findings (round v1, reference law)

Headline: time-growing workload τ = L/4, all five changes, one week. Exact tables are in the [summary](../results/optimization_round_v1/summary.md).

1. **The qualitative picture survives.**
   - NISQ plateaus: τ ≤ 1 at p = 10⁻³ and τ ≤ 2 at p = 10⁻⁴, with assumed effective mitigation.
   - STAR and the hybrid saturate: L = 8–12 (τ ≈ 2–3) at p = 10⁻³ and L = 12–20 (τ ≈ 3–5) at p = 10⁻⁴.
   - FT keeps scaling to the sampled ceiling τ = 8.
2. **Crossovers.** Encoded computation overtakes the NISQ effective line at about 2×10⁵ qubits for ε = 0.01 at both p. FT alone leads from about 6×10⁵ qubits (p = 10⁻³, ε = 0.01). At ε = 10⁻³ the FT take-off moves from 1.6×10⁷ to 6×10⁵ qubits at p = 10⁻³, and from 1.3×10⁷ to 3×10⁶ at p = 10⁻⁴.
3. **Which choices change the winner.** Out of 64 cells (both families, both p, both ε, eight budgets), counted against the v1 menu:

   | Change | Cells with a different leader | Mechanism |
   |---|---:|---|
   | FT delivery | 8 | Whole-group HWP (T-count 13.4M → 3.1M at strict 48×48), parallel ports and more factories (105 s → 2 s per shot); factories are about 57% of that footprint |
   | Richardson | 6 | Base r = 8τ with a K₄/r⁴ residual frees the Trotter budget; strict precision only, and it helps every encoded arm |
   | Lanes | 4 | Removes STAR's lone small-Q lead at p = 10⁻⁴ (becomes a tie) and helps the hybrid |
   | Separate distances | 4 | At p = 10⁻⁴ FT fits at d = 9–11 instead of 15; also removes the STAR-only cell |
   | Raw RUS | 1 | Rarely flips a leader, but lowers the STAR/hybrid exponent about 1.5×. For example, STAR 16×16 at τ = 1, ε = 0.01 goes from 10¹⁰·¹ to 10⁸·⁰ qubits. This holds only within our instrument; Toshio et al. App. E finds raw RUS slightly worse than PCEC ([RAW_RUS_NOTE.md](RAW_RUS_NOTE.md)) |
   | All five | 15 | Mostly FT separating earlier; STAR never leads alone |
4. **No STAR-only window remains.** Its only solo win in v1 (p = 10⁻⁴, ε = 0.01, Q = 10⁵) becomes a tie. With all five changes the hybrid is never below STAR (in partial menus it trails only in that one cell), so it is the strongest partially protected option.
5. **Conditional.** Raw-RUS gains rest on first-order coherent cancellation in our restricted
   instrument, which the paper's Appendix E does not reproduce ([RAW_RUS_NOTE.md](RAW_RUS_NOTE.md)).
   Delivery tile charges and concurrency caps are declared, not compiled. The Trotter envelope beyond τ = 4 and all lattices beyond 4×4 are transfers. The surface-code law is the reference A = 0.1, Λ = 10. The Stim calibration
   ([SURFACE_LAW.md](SURFACE_LAW.md)) gives Λ = 7.4 at p = 10⁻³, which means +4–6 in distance
   (about 1.4× area) for both encoded arms; round v2 propagates it. The NISQ effective overhead is assumed, and interacting bulk convergence is open.

## Design

`python/optimization_round.py` reuses the integrated_v1 instruments and ledgers without editing
them. Each change is a switch (`configs/optimization_round.json`):

| Switch | What it changes | Declared assumption |
|---|---|---|
| `raw_rus` | STAR/hybrid channels may use raw repeat-until-success in addition to noisy PCEC | Same restricted generator-fault instrument; raw RUS carries its coherent residue and retry tail in the error budget (`coherent_bias_bound + rus_tail_bound`), amplified exactly as before |
| `lanes` | Analog lanes chosen from 1–128 instead of fixed at 4 | 6 tiles per lane (as before); lanes ≤ L for the Jordan–Wigner termwise circuit (long Z strings limit disjoint supports) and ≤ N for plaquette rotations; more than 4 lanes require the fast data block |
| `ft_delivery` | FT T-consumption ports (1–128), factories and Hamming-weight batches up to 1024 optimized jointly; STAR lane tiles no longer charged to FT | 2 tiles per extra port; ports need the fast block; synthesized rotations are serial T sequences, so consumption is Σ_stages repeat·⌈rotations/ports⌉·T_per_rotation plus parallel fixed and Toffoli T; applies to the hybrid's fixed-T consumption too |
| `distances` | Data distance no longer clamped to the cultivation end distance (15) | Transfer from the d = 15 source into a data patch of distance d is charged 2·max(d, 15) cycles per T in the source service time; STAR receivers keep their measured growth targets |
| `richardson` | Two-branch Richardson extrapolation (r, 2r) with coefficients (−1/3, 4/3) | See below |

`v1_menu` (all switches off) reproduces `results/integrated_v1` (STAR, hybrid, FT; component and
zero-growth; τ = 1) and the NISQ rows of `results/physics_sweep_v1`. `--check` enforces the
replay to 10⁻⁹ in log₁₀ Q wherever v1 is below 10¹⁰ qubits. Above that, the new exposure prune
may drop a design that v1 recorded, which is irrelevant to capacity.

## Richardson accounting

For branches i with coefficients c_i, per-attempt variance multiplier F_i = e^{χ_i} and attempt
time t_i, optimal shot allocation gives a total time proportional to (Σ|c_i|√(F_i t_i))².
With equal F and t_2 = 2t_1 this is (1/3 + (4/3)√2)² = 4.92 base runtimes, as the cross-review
noted. The branches share hardware and code distance; every branch must satisfy the logical,
magic and synthesis budgets divided by ‖c‖₁ = 5/3 after its own Γ amplification, and
instrument biases add as Σ|c_i|·b_i. The residual Trotter error is K₄/r⁴, where K₄ is 1.1 times
the running maximum through τ of the measured |(4E(2r) − E(r))/3|·r⁴ over all signed same-patch
pairs (2×2, 3×3, 4×3, 4×4). At τ = 1 the 4×4 plaquette pairs give 0.0514 at both r = 8 and 16,
so the residual is cleanly r⁻⁴. The finite-step guard stores unsigned errors; its only valid
bound, (4|E(2r)| + |E(r)|)/3, cannot cancel the r⁻² term, so it is excluded. The base step
count keeps the r ≥ 8τ domain. For NISQ each branch pays its own mitigation exponent.

## Time-growing family

τ = L/4 for L = 4, 8, …, 32. At τ ∈ {1, 2, 4} the Trotter coefficient is exactly integrated_v1's
extended envelope. At other times it is the same rule, max(3×3 running maximum, 1.1 × every
measured 4×4/4×3 coefficient through τ), without the restriction to τ ∈ {1, 2, 4}. No patch
larger than 4×4 has been measured, and interacting bulk convergence is not established.

## Exposure prune

A STAR/hybrid candidate whose exposure exceeds χ = 120 is dropped. Such a design needs at least
10²⁴ qubits for one week under any attempt time in this model, so the prune cannot affect
capacity at Q ≤ 10⁹. Noisy PCEC is evaluated only where raw RUS stays below the prune. The
check asserts on every cached angle that PCEC dephasing is at least raw dephasing, which
justifies this gate.

## Reproduction

```bash
qsub pbs/optimization_round.pbs          # 128-core node; resumable
make optimization-round-figures           # figures and summary from saved rows
make optimization-round-check             # replay, budgets, fault ledgers, switch domains
```
