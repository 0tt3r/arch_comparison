# Surface-code law under standard uniform noise (SD6)

The study now uses the standard depolarizing circuit model (SD6) throughout: depolarizing error p after
every one- and two-qubit gate, p on every idle qubit in every layer (including measurement and reset
layers), and flip probability p after resets and before measurements. This is the convention of the
cultivation (Gidney et al.) and distillation (Litinski) estimates the study uses, and one-qubit gates are
no longer ten times better than two-qubit gates. The earlier calibration with one-qubit error p/10 on a
CZ schedule is kept in [SURFACE_LAW.md](SURFACE_LAW.md) and `results/surface_law_v1`.

## Method (`python/surface_law_sd6.py`, `results/surface_law_sd6_v1`)

Stim's rotated memory circuit (native CNOT) with SD6 noise applied layer by layer. Sampling, PyMatching
decoding, the fit, and the independent 100d holdout reuse `surface_law_calibration.py` unchanged, with
the same design as v1 (28 training cells at d = 3–9 for p = 10⁻³ and d = 3–7 for p = 10⁻⁴, 7.8×10⁷
training shots; 12 holdout cells, 1.8×10⁷ shots) and fresh seeds. Noiseless circuits fire no detectors
and every error model decomposes for matching (`validation.json`).

## Result

| p | Λ [95%] | A = A_X + A_Z | v1 (p/10 one-qubit) |
|---|---|---|---|
| 10⁻³ | 6.30 [6.18, 6.43] | 0.103 | Λ = 7.41, A = 0.077 |
| 10⁻⁴ | 43.5 [41.1, 46.1] | 0.186 | Λ = 89, A = 0.125 |

A single law is strongly rejected (deviance 680 on 5 degrees of freedom at 10⁻³; 341 on 3 at 10⁻⁴):
with this many shots, small changes in the distance-to-distance slope are resolved. The local slopes are
5.6, 7.0, and 6.0 for d = 3→5→7→9 at p = 10⁻³, so the fitted 6.30 is about 5% above the slope at the
largest distances and the extrapolated error at d ≈ 25 is optimistic by up to about 2×, roughly one code
distance. At p = 10⁻⁴ the local slopes are 41 and about 110 (few failures at d = 7), so the fit is
conservative there. The fit's own sensitivity (d ≥ 5 only, zero endpoint correction, the 100d holdout)
changes the required distance by at most 2.

Relative to the reference law (A = 0.1, Λ = 10), a protected volume of 10⁸ patch-rounds needs d = 27–33
instead of 21–25 for budgets 10⁻⁴–10⁻⁶ at p = 10⁻³, a factor 1.5–1.8 in patch area.

## STAR growth under SD6 (`python/growth_sd6.py`, `results/growth_sd6_v1`)

The 24 growth cells of `results/frontier_v1` (1 million shots each) repeated with SD6 noise; the circuit
is `star_growth_probe.growth_circuit` with its one-qubit and idle noise raised from p/10 to p and its 0.8p
per-round reset allowance removed (asserted source edits; noiseless circuit identical). Failures roughly
double where any occur: 5 → 25 at p = 10⁻³ gives 119 (X) and 819 (Z) failures, against 51 and 419, and
the union upper bound rises from 5.6×10⁻⁴ to 1.1×10⁻³. Cells from d_i = 9 and 13 at p = 10⁻³ now see a few
failures; all cells at p = 10⁻⁴ still see none.

## Not rerun

The lattice-surgery calibration ([SURGERY_LAW.md](SURGERY_LAW.md)) remains on the p/10 CZ schedule. It
supports a relative statement (isolated surgery fails at 0.73–1.15× the memory law under the same noise)
that is expected to carry over.
