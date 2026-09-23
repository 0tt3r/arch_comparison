# Published benchmarks: what is reproduced and what remains assumed

[TopQAD figure](../../figures/published_topqad.pdf) · [Experimental figure](../../figures/published_phasecraft.pdf) · [Methods](../../docs/PUBLISHED_BENCHMARKS.md)

The published examples are compatible with the resource laws once task, estimator, parallelism and accuracy are kept explicit. This pass does not fit an experimental mitigation exponent or replace the integrated capacity curves.

## TopQAD / STAR: Table IV

All 16 STAR footprints reproduce the published two-significant-digit values from (4L²+1)2d². 14 of 16 PEC-free runtimes agree within 2.3%. The 16 FT time-optimal entries agree within 2.8% with one T consumed per d syndrome cycles.

The replay supplies the paper’s distances, RUS clocks and T counts. It reconstructs the remaining area/time ledger, including both Hadamard-test quadratures, the multilevel time grid, integer steps and control clocks. The growth inferred from the total/PEC-free ratio is a diagnostic, not a forward validation of the growth fit.

| U | L | Hardware | PEC-free replay / published | Total time with anchor growth / published | Implied growth per state | Flag |
|---:|---:|---|---:|---:|---:|---|
| 4 | 4 | target | 1.001 | 1.001 | 3.2e-06 | none |
| 4 | 4 | desired | 0.996 | 0.996 | 1.6e-07 | none |
| 4 | 6 | target | 0.982 | 0.997 | 3.1e-06 | none |
| 4 | 6 | desired | 0.997 | 1.002 | 1.2e-07 | none |
| 4 | 8 | target | 0.977 | 1.049 | 2.8e-06 | none |
| 4 | 8 | desired | 1.003 | 0.994 | 2.1e-07 | none |
| 4 | 10 | target | 0.995 | 1.044 | 2.9e-06 | none |
| 4 | 10 | desired | 0.992 | 0.984 | 2.1e-07 | none |
| 8 | 4 | target | 0.988 | 0.971 | 3.2e-06 | none |
| 8 | 4 | desired | 0.993 | 1.003 | 1.1e-07 | none |
| 8 | 6 | target | 0.984 | 1.092 | 2.8e-06 | none |
| 8 | 6 | desired | 0.925 | 0.999 | no nonnegative solution | PEC-free replay differs by more than 5 percent; Published PEC multiplier is below the stated small-angle contribution |
| 8 | 8 | target | 1.005 | 1.091 | 2.9e-06 | none |
| 8 | 8 | desired | 1.561 | 1.551 | 1.8e-07 | Published clock and single-shot columns disagree; PEC-free replay differs by more than 5 percent |
| 8 | 10 | target | 0.995 | 1.096 | 2.8e-06 | none |
| 8 | 10 | desired | 0.998 | 0.990 | 1.9e-07 | none |

For target hardware, the inferred growth contribution clusters near 3×10⁻⁶. That numerical agreement is suggestive, but the inference uses the published runtime and cannot validate it independently. Fitting only U=4,L=4 and transferring the same growth parameter reproduces all seven other target-hardware total runtimes within 10%. This is a calibrated workload-transfer check, not independent gadget validation. Two desired-hardware rows remain inconsistent under the declared replay: U=8,L=8 lists 750 clocks, which predicts 470 ms rather than 297 ms; U=8,L=6 has a PEC factor below the angle-only estimate. Neither row is silently repaired. The preceding STAR supplementary-code endpoint returned HTTP 403 during this pass; a complete independent acceptance/growth/schedule reproduction remains open.

At L=8,U=4 under target hardware, STAR uses 230,000 qubits and 8.6 minutes. FT uses 570,000 qubits and 1.9 minutes in the space-optimal layout or 890,000 qubits and 21.2 seconds in the time-optimal layout. This is an area-versus-time tradeoff on the paper’s GSEE task, not an ordering for our local dynamics observable.

### Change one assumption on the fixed 4×4 GSEE workload

| Single change | Runtime [minutes] | Ratio to anchor |
|---|---:|---:|
| Published anchor | 21.31 | 1.000 |
| Clock only: 240 ns | 14.62 | 0.686 |
| Angle coefficient only: doubled | 24.06 | 1.129 |
| Growth only: our d9-to-d13 upper bound | 123.79 | 5.808 |

These isolate model sensitivities. The substituted growth bound has a different noise model, preparation and growth protocol; it is not a calibrated migration of the published hardware. All other clocks, distances and schedule choices stay fixed.

## Phasecraft / Quantinuum: direct data reduction

The pinned archive contains 3,200 shots per interaction value, 160 shots at each of 20 times, on 56 qubits. The mean doublon reconstructed from raw shots agrees with every saved raw mean to floating-point precision. The following RMS differences use all 20 hardware times; local entries also use all 28 sites.

| U | Observable | Method | Reference | RMS difference | Median published standard error |
|---:|---|---|---|---:|---:|
| 0 | mean_doublon | Raw | FLO | 0.0259 | 0.0054 |
| 0 | mean_doublon | TFLO | FLO | 0.0115 | 0.0059 |
| 0 | mean_doublon | Mitigated | FLO | 0.0064 | 0.0027 |
| 0 | local_doublon | Raw | FLO | 0.0433 | 0.0282 |
| 0 | local_doublon | TFLO | FLO | 0.0361 | 0.0312 |
| 0 | local_doublon | Mitigated | FLO | 0.0200 | 0.0144 |
| 4 | mean_doublon | Raw | TDVP chi=2048 | 0.0584 | 0.0054 |
| 4 | mean_doublon | TFLO | TDVP chi=2048 | 0.0378 | 0.0060 |
| 4 | mean_doublon | Mitigated | TDVP chi=2048 | 0.0362 | 0.0028 |
| 4 | local_doublon | Raw | TDVP chi=2048 | 0.0672 | 0.0279 |
| 4 | local_doublon | TFLO | TDVP chi=2048 | 0.0506 | 0.0313 |
| 4 | local_doublon | Mitigated | TDVP chi=2048 | 0.0412 | 0.0146 |

For U=4 the raw mean-doublon variance per shot is about 0.0035, compared with about 0.149 for an average individual site and the discussion model’s Bernoulli reference 0.09. Spatial averaging and correlations already change the shot budget substantially; smoothing also shares information across time. Thus a single-point gatewise-PEC estimate is not a prediction of the experimental shot count. The U=0 comparison uses training information. U=4 differences combine experimental, Trotter and finite-bond reference errors, and cannot be interpreted as held-out residual bias or as proof of failure of mitigation.

## Confidence after this pass

| Item | Evidence now | Remaining assumption |
|---|---|---|
| NISQ experimental footprint and shot count | Direct archive census and exact reconstruction of raw means | Different state, geometry and observable from our baseline |
| NISQ useful observable accuracy | Published processed curves and reference comparisons | U=4 exact truth and transfer to a chosen local target |
| STAR published area and PEC-free ledger | 16 areas and 14 consistent runtimes | Supplied RUS clocks and distances; full noisy injection/growth |
| STAR mitigation overhead | Inferred growth scale from published PEC ratio | Independent forward channel prediction |
| FT time-optimal throughput | 16 one-port timing checks | Published synthesis counts and private factory/layout optimizer |

Willow shot data are also public: [Zenodo 17881822](https://zenodo.org/records/17881822), an 8.8 GB archive with raw shots and processed expectations. Its README and metadata were checked; that archive was not reduced here. This corrects the earlier uncertainty about availability and makes a hardware-matched follow-up possible.
