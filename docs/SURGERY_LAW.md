# Isolated lattice-surgery calibration — 23 September 2026

**At $R_m=d$, the memory-volume law is adequate at the requested factor-of-two level in every well-resolved cell. Shortening the merge can invalidate it.**

**Setup.** 678 million fixed shots, 158 cells; Stim 1.16.dev0, PyMatching 2.3.1, ordinary MWPM. Four CZ slots, H noise, preparation/readout flips and idle noise exactly match [SURFACE_LAW.md](SURFACE_LAW.md); rounds take 240 ns. Rotated patches merge vertically for ZZ and horizontally for XX. CNOT uses an L-shaped control/ancilla/target arrangement: ZZ, XX, then ancilla-Z readout. $R_{\rm pre}=R_{\rm post}=d$. Inputs, including the ancilla, are **ideally encoded**; ideal Bell references and terminal checks probe the channel. Bridge preparation and all surgery readout are noisy. This does not price physical initial encoding.

**Measured / $V p_L$d$$, $R_m=d$** [95% Clopper–Pearson intervals]:

| $p$ | $d$ | ZZ, $l=0$ | ZZ, $l=1$ | ZZ, $l=3$ | Full CNOT probe |
|---|---:|---:|---:|---:|---:|
| 10⁻³ | 3 | .94 [.93,.96] | .94 [.93,.96] | .96 [.95,.97] | .92 [.92,.93] |
| 10⁻³ | 5 | .93 [.91,.96] | .92 [.90,.94] | .91 [.89,.93] | .91 [.90,.92] |
| 10⁻³ | 7 | .87 [.82,.93] | .86 [.81,.91] | .88 [.84,.93] | .84 [.81,.87] |
| 10⁻³ | 9 | .79 [.73,.85] | .78 [.72,.84] | .73 [.68,.78] | — |
| 10⁻⁴ | 3 | .85 [.80,.90] | .91 [.86,.96] | .94 [.90,.99] | .87 [.85,.90] |
| 10⁻⁴ | 5 | 1.13 [.83,1.49] | 1.00 [.74,1.32] | 1.00 [.77,1.28] | 1.14 [.98,1.32] |
| 10⁻⁴ | 7† | <5.45 | <4.78 | 3.03 [.62,8.84] | 1.61 [.44,4.12] |

† Fewer than 30 failures: unresolved; “<” denotes a zero-failure interval upper endpoint. Ratios condition on the saved memory fit, excluding its uncertainty and model error. Its $A_X+A_Z$ is a conservative component sum. $l=0$ means a **one-data-column seam**; $l=1,3$ gaps have $ld$ columns. Area-equivalent volumes are $V_{ZZ}=4d+[2+\max(1,ld)/d]R_m$, $V_{CX}=5d+6R_m+2R_m/d$. Nominal tile counts, split-idle equivalents, all $R_m$, and separate outcome/data failures are in [raw.csv](../results/surgery_v1/raw.csv).

**Time-like component.** Define $q_t$ as wrong reported parity without a data-channel failure; it is an outcome classification, not a microscopic fault-location label. Exactly $q_{ZZ}=q_{\rm space}+q_t$. Fit $q_t=B_{d,l}\Lambda_t^{-(R_m+1)/2}$: shared slopes are 23.16 [22.75,23.57] at 10⁻³ and 264.5 [239.8,292.8] at 10⁻⁴ (95% profile intervals). **The former is rejected** (goodness-of-fit $p=2.5\times10^{-33}$); use geometry-specific fits. For example, $d=5,l=3,p=10^{-3}$ gives $B=5.75$, $\Lambda_t=18.13$ [17.11,19.24]. One sparse local fit has only a lower bound; see [fit interpretation](../results/surgery_v1/fit_interpretation.json).

At 10⁻³, smallest tested $R_m$ with separated $q_t<q_{\rm space}$ intervals, ordered by $l=0,1,3$, are **(3,3,3), (3,5,5), (5,5,7), (7,7,9)** for $d=3,5,7,9$. At 10⁻⁴ they are **(3,3,3), (5,5,5)**; $d=7$ is unresolved. For shorter merges, use $q_{ZZ}\approx Vp_L+B_{d,l}\Lambda_t^{-(R_m+1)/2}$, retaining local-fit uncertainty; memory alone underestimates the worst tested cell by 57×.

**Conclusion.** Retain the memory-volume law for these isolated $d$-round gadgets: no >2× penalty is supported in resolved cells. With pre/post storage made noiseless, the CNOT core gives ratios 0.77–1.15 and $V=6d+2$. Its duration is **$2d\times240\,\mathrm{ns}+300\,\mathrm{ns}$**, so exactly two clocks omit 9–21% at $d=7$–3, assuming a ready ancilla and prestaged gap resets; otherwise add preparation latency. Common multiplier fits show finite-distance deviations; do not replace them with a universal fitted constant. This supports isolated adjacent CNOT accounting, not congestion-free long-range layers. Simultaneous routing, twists/Y surgery, correlated noise, physical encoding/decoder latency, and $d>9$ remain uncovered.

**Checks and reproduction.** All 14 original-memory controls reproduce the saved fit within cell intervals. Zero-noise, logical-flip, single-fault decoding, source/hash and complete circuit-rebuild audits pass. A rejected detector-basis run is archived separately, never pooled. Run `~/soft/mambaforge/envs/stim/bin/python python/surgery_calibration.py audit` (`analyze`/`plot` likewise). [Figure](../figures/surgery_law.png), [fits](../results/surgery_v1/fit.json), [manifest](../results/surgery_v1/manifest.json), and [detailed conventions](../results/surgery_v1/METHOD.md) retain the evidence.
