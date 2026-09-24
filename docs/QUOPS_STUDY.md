# QUOPS: a second, fine-grained workload

QUOPS (Proctor et al., [arXiv:2609.12146](../refs/2609.12146v1.pdf)) scores a machine by the
largest random circuit it can run with mean process polarization at least 1/√e. The circuits
alternate a CNOT layer on ⌊w/2⌋ uniformly random pairs with a layer of R_P(θ) rotations,
P ∈ {X, Y, Z} and θ uniform on [0, 2π). Size counts one per rotation and two per CNOT, and the
score is restricted to w² ≤ s ≤ w³ ([specification](QUOPS_SPEC.md)).

The workload complements the Hubbard task in two ways. Its size axis is nearly continuous, so
capacity has no lattice-grid artifacts. Its rotation angles are of order one, so STAR loses the
small-angle discount it enjoys on product formulas. The prediction from the exposure picture is
that STAR's advantage over full fault tolerance should shrink sharply. That prediction is the test.

## Model (`python/quops_study.py`, `configs/quops_study.json`)

For every shape (w, s) and arm we compute the fewest physical qubits Q_min; the capability at a
budget Q is the largest s with Q_min ≤ Q. Two success modes are reported:

- **Standard.** Mean polarization γ ≥ e^{-1/2}, with the mirror test (declared: 100 mirror circuits × 100
  shots at twice the forward size, within the paper's experimental range of 2×10³–2×10⁵) completed within one week.
- **Mitigated.** Our extension of the paper's α threshold: the test's 1/γ² overhead must fit in one
  week, using replicas.

Polarization is the embedded process-infidelity product, γ = Π(1 − a_w ε_i), with
a_w = 4^w/(4^w − 1) and ε_i each operation's process infidelity. This is exact to 0.2% on w = 2–8
([validation](QUOPS_POLARIZATION.md)); the product of local gate polarizations would be up to 24%
pessimistic.

| Arm | Cost model |
|---|---|
| Physical qubits | Measured nearest-neighbor routing per layer ([routing](QUOPS_ROUTING.md)); the lower-loss router of Cirq RouteCQC and a constructive sorting network at each width. The study's Stim convention (as in the surface-code calibration): DEPOLARIZE2(p) per two-qubit gate, DEPOLARIZE1(p/10) per one-qubit gate, DEPOLARIZE1(p) per idle 10 ns slot; process infidelity equals total Pauli error probability. No SPAM (the mirror test excludes it). |
| STAR | The Hubbard study's generator-fault injection instrument averaged over the wrapped angle distribution (uniform on [−π/8, π/8]); raw RUS or noisy PCEC; seeds 9 and 25 with receivers 13 or 25; lanes 1–1024. Cliffords by lattice surgery at two code clocks per CNOT layer (declared). |
| FT | Gridsynth with quadratic polarization charge ε² per rotation (the paper's convention); cultivated T states with optional 15-to-1 cleanup; ports and factories optimized; two code clocks per CNOT layer (declared). |

Both encoded arms use the Stim-calibrated surface-code law ([SURFACE_LAW.md](SURFACE_LAW.md)).

## Calibration against the paper

- The physical arm is anchored by routing measured on the same circuit family, not by fitting
  the measured scores. Willow's 216 was obtained at its own error rates and connectivity.
- For the paper's planar FT projection shapes (Table III), our FT needs 4–7× more physical qubits.
  The calibrated Λ = 7.4 raises the distance relative to a uniform-depolarizing law, and the fast
  data block and 1024-qubit cultivation sources add area. We report the ratio rather than tune it.

## Known limitations

Lattice-surgery routing of random long-range CNOTs is assumed congestion-free at two code clocks
per layer. The isolated-gadget calibration ([SURGERY_LAW.md](SURGERY_LAW.md)) supports the two clocks
(2d rounds plus 300 ns) and the memory-volume law (0.73–1.15× at R_m = d) for adjacent CNOTs; congestion
is not covered. The one-sided 95% confidence rule and SPAM references are not modeled; with 10⁴ executions
the statistical margin is small. The STAR channel is the restricted instrument of the Hubbard
study, with the raw-RUS caveat of [RAW_RUS_NOTE.md](RAW_RUS_NOTE.md).

## Results (p = 10⁻³, standard test)

Physical qubits score about 310 at any size (another workshop implementation gets about 150; Willow
measured 216 at its own error rates). STAR overtakes them at 10⁴ qubits, leads alone between 10⁴ and
2.5×10⁴, and saturates at 7×10³. FT is ahead to stay from 3×10⁴ qubits, reaches 10⁴ at 5×10⁴ and
5.5×10⁶ at 10⁶ (about 800× STAR), then plateaus near 1.1×10⁷ from 3×10⁶ to 1.3×10⁸ qubits, where
cultivated states (2×10⁻⁹) no longer suffice and cleanup factories are needed. With the published
distillation factories of Litinski (2019) added to the source menu
([QUOPS_SENSITIVITY.md](QUOPS_SENSITIVITY.md)) the plateau disappears: the RSA-2048 target
(2.5×10⁸) is reached with 2.5×10⁸ qubits at p = 10⁻³ and 10⁷ at p = 10⁻⁴; with cultivation alone
it is reached only at p = 10⁻⁴ with about 8×10⁸. (A declared 100× better cultivated source,
`ft_better_source`, reaches it with 4×10⁷ at p = 10⁻³.) The dashed curves are a declared
1/γ² throughput scenario, not the cost of certifying polarization. The shaded band in the figure uses the factories re-optimized under the calibrated law
([FACTORY_CALIBRATION.md](FACTORY_CALIBRATION.md)), which reach the RSA-2048 target with 4×10⁸ qubits at p = 10⁻³. Figure:
[figures/quops_capability.pdf](../figures/quops_capability.pdf).
