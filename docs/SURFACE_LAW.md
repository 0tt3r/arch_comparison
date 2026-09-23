# Surface-law calibration — 23 September 2026

**The reference suppression Λ=10 at p=10⁻³ is optimistic at larger distances.** The exponential is an approximate resource model: finite-distance lack of fit prevents treating its statistical intervals as prediction guarantees.

**Experiment.** 78 million fixed training shots (R=10d,20d), plus 18 million independent R=100d holdout shots; d=3,5,7,9 at p=10⁻³ and d=3,5,7 at 10⁻⁴ (holdout d=3,5). [Stim](https://github.com/quantumlib/Stim) 1.16.dev0 and [PyMatching](https://pymatching.readthedocs.io/en/stable/) 2.3.1, ordinary MWPM, one serial process/node. Circuit controls and distance checks passed.

**C6/C7 implementation.** CZ depolarizing p; every H depolarizing p/10; preparation/readout flips p; idle depolarizing p per unused CZ slot and per data-qubit measurement/reset slot. Four CZ slots + 100 ns measurement + 100 ns reset give nominal 240 ns. **Basis changes are fused into these slots; their additional physical duration is unmodeled.** Existing growth controls use different idle conventions and were not pooled.

**Fit.** q_b(d)=A_b Λ^{−(d+1)/2}; P_b(R)=[1−(1−2q_b)^{R+δ_db}]/2. X/Z have separate amplitudes, common Λ, and fitted endpoint corrections δ. AΣ=A_X+A_Z is the conservative Pauli-component sum convention; basis-mean A=AΣ/2. It is not a directly measured any-Pauli rate.

| p | Λ [95% profile interval] | AΣ [95% profile interval] |
|---|---:|---:|
| 0.001 | 7.41 [7.2, 7.62] | 0.0773 [0.0724, 0.0825] |
| 0.0001 | 89.2 [63.5, 146] | 0.125 [0.0618, 0.336] |

Cell intervals are exact Clopper–Pearson; profile intervals are asymptotic and conditional on this law. At 10⁻³ the training-law lack-of-fit p=0.0133; 100d also rejects a single d=3–9 law (p=8.8e-11). A **post hoc d≥5,100d sensitivity fit** gives Λ=7.97 [7.82, 8.11], AΣ=0.097 [0.0914, 0.103] (endpoint correction neglected; fit p=0.17). At 10⁻⁴,d=7 all four cells have zero failures; each million-shot R=140 cell gives q_b<2.64×10⁻⁸ (two-sided 95% interval upper endpoint). Sparse endpoint parameters hit bounds. pΛ is an effective fit threshold, not a measured threshold crossing.

**Conditional resource extrapolation.** Smallest odd d≥3 with VAΣΛ^{−(d+1)/2}≤ε, V=tiles×cycles. Entries are assumed→fitted d (d² ratio); every entry extrapolates beyond sampled distances. JSON also gives statistical intervals, fit-choice sensitivity, and basis-mean comparisons.

| p | V | ε=10⁻⁴ | ε=10⁻⁵ | ε=10⁻⁶ |
|---|---:|---:|---:|---:|
| 0.001 | 1e+08 | 21→25 (1.42) | 23→27 (1.38) | 25→29 (1.35) |
| 0.001 | 1e+11 | 27→31 (1.32) | 29→35 (1.46) | 31→37 (1.42) |
| 0.001 | 1e+14 | 33→39 (1.40) | 35→41 (1.37) | 37→43 (1.35) |
| 0.0001 | 1e+08 | 11→11 (1.00) | 11→13 (1.40) | 13→13 (1.00) |
| 0.0001 | 1e+11 | 13→15 (1.33) | 15→15 (1.00) | 15→17 (1.28) |
| 0.0001 | 1e+14 | 17→17 (1.00) | 17→19 (1.25) | 19→19 (1.00) |

**Scope.** The same law applies to STAR and FT, so it creates no arm-specific advantage; different volumes and rounding can still change their relative costs. It does not cancel against NISQ. Memory only: surgery, correlated/coherent hardware noise, and high-distance performance are unvalidated.

**Reproduce:** `~/soft/mambaforge/envs/stim/bin/python python/surface_law_calibration.py fit` then `audit`. [Evidence and manifests](../results/surface_law_v1/fit.json) retain seeds, circuits, counts, intervals and PBS scripts.
