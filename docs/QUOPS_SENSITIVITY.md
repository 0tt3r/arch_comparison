# Angle sweep and delivery sensitivity on random circuits

`python/quops_sensitivity.py` reruns the QUOPS model (`python/quops_study.py`) under the standard
polarization threshold with one change at a time; `python/plot_quops_sensitivity.py` draws
[figures/quops_angle_sweep.pdf](../figures/quops_angle_sweep.pdf). Results:
[results/quops_sensitivity_v1](../results/quops_sensitivity_v1/).

## Angle sweep: the per-radian law under one criterion

STAR's injected angles are drawn uniformly on [−a, a] for a = 0.003, 0.01, 0.03, 0.1, and π/8
(the QUOPS value). Physical-qubit and FT costs do not depend on the angles, so the sweep changes
only STAR's exposure, with the circuit family, success criterion, routing, and delivery held fixed.

| Mean injected angle | STAR saturation, p = 10⁻³ | p = 10⁻⁴ | s·p·⟨\|θ\|⟩ (10⁻³ / 10⁻⁴) |
|---:|---:|---:|---:|
| 0.0015 | 2.7×10⁵ | 2.7×10⁶ | 0.41 / 0.41 |
| 0.005 | 8.1×10⁴ | 8.0×10⁵ | 0.41 / 0.40 |
| 0.015 | 2.8×10⁴ | 2.8×10⁵ | 0.42 / 0.42 |
| 0.05 | 1.5×10⁴ | 9.0×10⁴ | 0.73 / 0.45 |
| 0.196 (QUOPS) | 7.1×10³ | 3.3×10⁴ | 1.4 / 0.64 |

Below about 0.03 rad the saturation follows **s_sat ≈ 0.41/(p⟨|θ|⟩)** to within 3% across two decades
of angle and both error rates: slope −1, as the exposure picture predicts. At larger angles the
injected error grows more slowly than linearly and STAR does up to three times better than the
small-angle law.

The angle sets the height of STAR's plateau, not the crossover. At p = 10⁻³, shrinking the angles
to Hubbard size raises STAR's ceiling about forty times, yet full fault tolerance still overtakes it
for good at 6×10⁴ qubits rather than 3×10⁴: FT wins by scaling past the plateau, and the crossover
sits where FT's fixed area is paid. A narrow STAR-only window (about 1–5×10⁴ qubits) exists at
p = 10⁻³ for every angle scale and at no angle scale at p = 10⁻⁴.

## Delivery sensitivity (p = 10⁻³)

| Change | STAR saturation | FT at 10⁶ qubits | FT ahead to stay | STAR-only window |
|---|---:|---:|---:|---|
| Baseline | 7.1×10³ | 5.5×10⁶ | 3.2×10⁴ | 1.0–2.5×10⁴ |
| Surgery time ×2 | 7.1×10³ | 5.5×10⁶ | 3.2×10⁴ | 1.0–2.5×10⁴ |
| Surgery congestion ∝ √(w/16) | 7.1×10³ | 5.5×10⁶ | 3.2×10⁴ | 1.0–2.5×10⁴ |
| 3 tiles per lane | 7.1×10³ | 5.5×10⁶ | 3.2×10⁴ | 0.8–2.5×10⁴ |
| 12 tiles per lane | 7.1×10³ | 5.5×10⁶ | 3.2×10⁴ | 1.3–2.5×10⁴ |
| 8 tiles per extra port | 7.1×10³ | 4.8×10⁶ | 3.2×10⁴ | 1.0–2.5×10⁴ |

The random-circuit conclusions are insensitive to the declared delivery models. Rotation layers
(about 25 d cycles for FT, and the analog service time for STAR) dominate the lattice-surgery CNOT
layers (2 d), so even a congestion penalty that grows as √w changes no capability at two significant
figures. Lane charges move only the edges of STAR's window; port charges cost FT about 13% at 10⁶
qubits.

## T source: published distillation factories

The FT source menu adds the protocols of Litinski, Quantum 3, 205 (2019), Tables 1–2, at the
matching p (published qubits, cycles per output, and output error, plus one delivery tile) to
cultivation with optional cleanup.

| Q | FT, cultivation (p = 10⁻³) | + factories | FT, cultivation (p = 10⁻⁴) | + factories |
|---:|---:|---:|---:|---:|
| 10⁵ | 8.5×10⁴ | 8.5×10⁴ | 9.4×10⁵ | 2.1×10⁶ |
| 10⁶ | 5.5×10⁶ | 5.5×10⁶ | 1.1×10⁷ | 2.6×10⁷ |
| 10⁷ | 1.1×10⁷ | 1.2×10⁷ | 1.1×10⁷ | 2.6×10⁸ |
| 10⁸ | 1.1×10⁷ | 1.2×10⁸ | 4.8×10⁷ | 9.1×10⁸ |
| 10⁹ | 7.5×10⁷ | 9.7×10⁸ | 3.5×10⁸ | 4.9×10⁹ |

Below 3×10⁶ qubits at p = 10⁻³ cultivation remains optimal and nothing changes. Above it the
factories remove the cultivation plateau, and the RSA-2048 target is reached with 2.5×10⁸ qubits
(p = 10⁻³) and 10⁷ (p = 10⁻⁴). At p = 10⁻⁴ the small single-level factories help at every size.

With the factories re-optimized under the calibrated surface-code law
([FACTORY_CALIBRATION.md](FACTORY_CALIBRATION.md)), the plateau still disappears: FT scores 9.9×10⁷
at 10⁸ qubits and 7.4×10⁸ at 10⁹ (p = 10⁻³), and the RSA-2048 target needs 4×10⁸ qubits instead of
2.5×10⁸. At p = 10⁻⁴ the calibrated and published menus agree to 10⁷ qubits and reach the target
with 10⁷. The shaded band of the capability figure uses the calibrated menu.
