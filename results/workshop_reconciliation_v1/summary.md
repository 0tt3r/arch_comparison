# Workshop models at a shared operating point

4 × 4, tau = 1 and 2, p = 0.001, 0.01 absolute error, 10 ns elementary gates, 100 ns measurement, 240 ns syndrome cycle, one time point, one-week budget. All runtimes use the same Gaussian doublon planning count at D*=0.1 and statistical allowance 0.004. This is a common estimator assumption, not a validation of native mitigation bias.
The common step count covers the larger of our plaquette and termwise planning coefficients using 25% of the error budget. It is transferred to the supplied cost models to isolate counts; their unspecified/different circuits and observables are not thereby accuracy-certified. Their native step estimates are retained in the CSV.
The table preserves native gate counts, register sizes, layouts, delivery, logical criteria, and source menus. Only task settings, steps, elementary clocks, and the final shot driver are normalized. Our source count is selected for throughput at Q=1e6, matching Martin's objective; IQUS retains its native point construction. IQUS uses the full circuit rather than a physical-cone discount. Its fixed signal is explicitly set. Adapters do not edit either supplied package.

| tau | r | Package | Arm / scenario | G2 or R | T / shot | Q / copy | seconds / shot | chi | log10 seconds at Q=1e6 |
|---:|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | 11 | ours | nisq: effective | 4749.827160493825 | — | 36 | 7.387e-06 | 5.003 | -3.068 |
| 1 | 11 | ours | nisq: two-qubit PEC | 4750 | — | 36 | 7.387e-06 | 19.001 | 3.012 |
| 1 | 11 | ours | star: count | 2443 | — | 58,482 | 0.05026 | 2.607 | 2.937 |
| 1 | 11 | ours | star: angle | 2443 | — | 58,482 | 0.05026 | 0.482 | 2.015 |
| 1 | 11 | ours | ft: individual | 764 | 51244 | 79,066 | 0.5856 | 0.000 | 3.023 |
| 1 | 11 | ours | ft: HWP 8 | 764 | 28602 | 88,768 | 0.3294 | 0.000 | 2.811 |
| 1 | 11 | ours | hybrid: angle | 764 | 1584 | 61,974 | 0.03206 | 0.295 | 1.764 |
| 1 | 11 | martin | nisq: native attenuation rescaling | 3520.0 | — | 32 | 6.7e-06 | 12.674 | 0.170 |
| 1 | 11 | martin | star: native count | 2288 | — | 28,800 | 0.001861 | 0.915 | 0.470 |
| 1 | 11 | martin | star: canonical count | 2288 | — | 28,800 | 0.001861 | 2.441 | 1.133 |
| 1 | 11 | martin | star: angle | 2288 | — | 28,800 | 0.001861 | 0.476 | 0.279 |
| 1 | 11 | martin | ft: native individual synthesis | 2288 | 148720 | 62,256 | 0.05608 | 0.000 | 1.879 |
| 1 | 11 | iqus | nisq: full-circuit two-qubit PEC | 2640.0 | — | 48 | 1.42e-06 | 10.561 | -1.245 |
| 1 | 11 | iqus | star: canonical count | 880.0 | — | 21,632 | 0.0006864 | 0.939 | -0.084 |
| 1 | 11 | iqus | star: angle | 880.0 | — | 21,632 | 0.0006864 | 0.419 | -0.309 |
| 1 | 11 | iqus | ft: native HWP | 880.0 | 13967 | 125,908 | 0.01268 | 0.000 | 1.593 |
| 2 | 23 | ours | nisq: effective | 9931.456790123455 | — | 36 | 1.523e-05 | 10.461 | -0.383 |
| 2 | 23 | ours | nisq: two-qubit PEC | 9932 | — | 36 | 1.523e-05 | 39.731 | 12.328 |
| 2 | 23 | ours | star: count | 5107 | — | 71,442 | 0.1161 | 5.449 | 4.652 |
| 2 | 23 | ours | star: angle | 5107 | — | 58,482 | 0.1051 | 0.972 | 2.548 |
| 2 | 23 | ours | ft: individual | 1580 | 110752 | 79,066 | 1.265 | 0.000 | 3.358 |
| 2 | 23 | ours | ft: HWP 8 | 1580 | 61641 | 88,768 | 0.7095 | 0.000 | 3.144 |
| 2 | 23 | ours | hybrid: angle | 1580 | 3312 | 61,974 | 0.06666 | 0.591 | 2.211 |
| 2 | 23 | martin | nisq: native attenuation rescaling | 7360.0 | — | 32 | 1.39e-05 | 26.498 | 6.491 |
| 2 | 23 | martin | star: native count | 4784 | — | 36,992 | 0.00441 | 1.914 | 1.379 |
| 2 | 23 | martin | star: canonical count | 4784 | — | 36,992 | 0.00441 | 5.104 | 2.765 |
| 2 | 23 | martin | star: angle | 4784 | — | 36,992 | 0.00441 | 0.959 | 0.964 |
| 2 | 23 | martin | ft: native individual synthesis | 4784 | 325312 | 62,256 | 0.1226 | 0.000 | 2.219 |
| 2 | 23 | iqus | nisq: full-circuit two-qubit PEC | 5520.0 | — | 48 | 2.86e-06 | 22.081 | 4.062 |
| 2 | 23 | iqus | star: canonical count | 1840.0 | — | 28,800 | 0.001656 | 1.963 | 0.875 |
| 2 | 23 | iqus | star: angle | 1840.0 | — | 28,800 | 0.001656 | 0.842 | 0.388 |
| 2 | 23 | iqus | ft: native HWP | 1840.0 | 30672 | 151,648 | 0.0313 | 0.000 | 2.052 |

## Identical STAR ledger, different laws

| tau | rotations | total angle | law | log variance multiplier |
|---:|---:|---:|---|---:|
| 1 | 2443 | 64.0 | Martin native: chi = 0.4 p R | 0.9772 |
| 1 | 2443 | 64.0 | IQUS native: chi = 1.066 p R | 2.6042 |
| 1 | 2443 | 64.0 | Canonical 2024 inverse | 2.6066 |
| 1 | 2443 | 64.0 | Angle + floor | 0.4817 |
| 2 | 5107 | 128.0 | Martin native: chi = 0.4 p R | 2.0428 |
| 2 | 5107 | 128.0 | IQUS native: chi = 1.066 p R | 5.4441 |
| 2 | 5107 | 128.0 | Canonical 2024 inverse | 5.4489 |
| 2 | 5107 | 128.0 | Angle + floor | 0.9723 |

The two native count-law coefficients are themselves different: 0.4 p R in Martin versus approximately 1.066 p R in IQUS and the canonical Akahoshi inverse. Therefore choosing “count law” alone does not reconcile the implementations.
Native rotations per step and injection schedules also differ. The supplied models are useful alternative scenarios, not independent replications of one numerical result. The shared-ledger table is an algebra check, not an empirical collapse plot.
Martin retains its single-attenuation NISQ estimator and T2 convention. Its FT logical criterion omits factory wait exposure where supply is slower. IQUS retains its unvalidated HWP workspace/synthesis accounting and factory table; its quoted FT runtime is diagnostic. No qLDPC result is imported into the nearest-neighbor hardware comparison.
The IQUS native Trotter calibration belongs to the spin-correlator/dimer task and is not independent calibration of our doublon quench. Its reconciliation preset retains an old four-step workload and does not control all the switches it labels.
Classical ED capacity is not a classical impossibility boundary. No frontier line or Lieb-Robinson certificate is imported without a compatible initial state, boundary-error theorem, and classical-method scope.

Source hashes are recorded in comparison.json; source_snapshot.zip preserves only the supplied Python modules used for reproduction. Run `python python/workshop_reconciliation.py`; when the supplied folders are absent, or with `--snapshot`, the adapter restores that snapshot under build/. The independent repository inside fh_iqus_idears is left untouched.
