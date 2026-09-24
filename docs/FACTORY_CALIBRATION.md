# Distillation factories under the calibrated surface-code law

The published factory tables of Litinski, Quantum 3, 205 (2019), price every patch, merge, and
storage step with one logical-error law, 0.1(100p)^((d+1)/2): our reference law (A = 0.1,
Λ = 10 at p = 10⁻³). Our Stim calibration gives A = 0.077, Λ = 7.4 at p = 10⁻³ and A = 0.125,
Λ = 89 at 10⁻⁴ ([SURFACE_LAW.md](SURFACE_LAW.md)). This note re-prices the factories on the same
footing as the rest of the study.

## Method (`python/factory_calibration.py`, `results/factory_calibration_v1/`)

Litinski's own model (github.com/litinski/magicstates, commit 8fffa48, Apache 2.0) is vendored in
`python/third_party/litinski_magicstates`. The only change to it is removing a stray `breakpoint()`
that aborts batch runs. The script replaces its `plog` with the calibrated law, keeps its injection
model (faulty π/8 rotations, error p/3 of each Pauli type), and evaluates at 40-digit precision:

- one-level 15-to-1 over all odd (d_X, d_Z, d_m) with d_X ≤ 31 and d_Z, d_m ≤ 19;
- two-level (15-to-1)×(15-to-1) over 8 (p = 10⁻³) or 5 (p = 10⁻⁴) level-1 designs, level-2
  d_X2 = 17–55 (11–35), d_Z2 = 7–25 (5–15), d_m2 = d_Z2 or d_Z2 + 2, and 4 or 6 level-1 factories;
- (15-to-1)×(20-to-4) in 45 steps around each published design (it is about 20 times slower).

That is 9255 designs. The menu keeps, per decade of output error, the Pareto designs with the
smallest footprint and the smallest space-time per output: 26 at p = 10⁻³ and 34 at 10⁻⁴.

**Validation.** With the reference law the vendored code reproduces all 11 published rows used
here (output error, qubits, cycles); the small-footprint (15-to-1)₉,₅,₅×(15-to-1)₂₁,₉,₁₁ row uses a
different layout and is excluded.

## Results

Published designs evaluated with the calibrated law (p = 10⁻³):

| Design | Published error | Calibrated error |
|---|---:|---:|
| (15-to-1)₁₇,₇,₇ | 4.5×10⁻⁸ | 2.4×10⁻⁷ |
| (15-to-1)⁶₁₁,₅,₅×(15-to-1)₂₅,₁₁,₁₁ | 2.7×10⁻¹² | 1.4×10⁻¹⁰ |
| (15-to-1)⁶₁₃,₅,₅×(15-to-1)₂₉,₁₁,₁₃ | 3.3×10⁻¹⁴ | 3.1×10⁻¹² |
| (15-to-1)⁶₁₇,₇,₇×(15-to-1)₄₁,₁₇,₁₇ | 4.5×10⁻²⁰ | 1.9×10⁻¹⁷ |
| (15-to-1)⁴₁₃,₅,₅×(20-to-4)₂₇,₁₃,₁₅ | 2.6×10⁻¹¹ | 9.0×10⁻¹⁰ |

At p = 10⁻⁴ the changes are 1.6–11×. Re-optimizing the distances recovers most of the loss: at the same
output error, two-level factories at p = 10⁻³ need 17–45% more qubits (for example 6.8×10⁻¹³ at
35,902 qubits against the published 2.7×10⁻¹² at 30,700).

**Consequences.** On the Hubbard task (`calibrated_factories` in round v3) every crossover matches the
published-factory run to within 3% ([OPTIMIZATION_ROUND.md](OPTIMIZATION_ROUND.md)). On random
circuits the calibrated factories still remove the cultivation plateau; the RSA-2048 target is
reached with 4×10⁸ qubits at p = 10⁻³ (published factories: 2.5×10⁸) and 10⁷ at 10⁻⁴
([QUOPS_SENSITIVITY.md](QUOPS_SENSITIVITY.md)).

**Scope.** This replaces the logical-error law inside Litinski's model; it is not a circuit-level
simulation of a factory. Injection errors, layouts, and the 15-to-1 protocol structure are his.
