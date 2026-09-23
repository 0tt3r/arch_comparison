# Willow with a matched interacting reference

[Figure](../../figures/willow_matched_reference.pdf) · [Methods](../../docs/WILLOW_REFERENCE.md)

We independently compute continuous dynamics and ideal experimental circuits for the 4×4, U=4, zero-flux center-hole input. All 24 nonzero saved times and all three step counts are included. This adds the missing interacting reference to the preceding archive analysis; it does not refit the published mitigation pipeline or extrapolate its performance.

## What the circuit reconstruction establishes

The Appendix D2 swap-network order reproduces all 288 archived U=0 circuit-reference points at both 4×4 and 6×6 and both fluxes, with maximum local-Z discrepancy 5.6e-15. This requires diagonal symmetry averaging: the unsymmetrized Z discrepancy reaches 0.345. The onsite opposite-spin ZZ-derived doublons reproduce within 7.5e-9, the precision of the archived quantities. The U>0 onsite evolution and complete small-system state were separately checked against a full Fock-space matrix construction.

## Observed differences for the processed experimental estimates

All values below are local-doublon RMS differences across 24 times and 16 sites. The experimental discrepancy contains finite sampling and processing error; it is not an estimate of ensemble residual bias.

| Steps | Processed vs ideal circuit | Ideal circuit vs continuous | Processed vs continuous | Fraction of local points within 0.01 of continuous |
|---:|---:|---:|---:|---:|
| 1 | 0.004688 | 0.044808 | 0.044288 | 0.221 |
| 2 | 0.007641 | 0.019731 | 0.015813 | 0.435 |
| 3 | 0.009462 | 0.013150 | 0.010301 | 0.714 |

The two contributions respond oppositely to refinement: the ideal dynamics improves from one to three steps, while the processed-data discrepancy from the ideal circuit grows. At three steps, the total RMS is smaller than the deterministic circuit error alone because the signed contributions partially cancel. This is a measured depth-versus-accuracy tradeoff for this pipeline, not a projection beyond three steps.

The signed identity is

    observed − continuous = (observed − ideal circuit) + (ideal circuit − continuous).

Consequently, the squared RMS total includes twice the mean product of the two signed contributions. The saved decomposition retains this cross term, so accidental cancellation is not credited as a generally reliable error budget.

## Approximate classical references

| MPS bond dimension | Local doublon RMS vs continuous | Maximum local difference |
|---:|---:|---:|
| 1024 | 0.00399118 | 0.0178296 |
| 2048 | 0.00213554 | 0.00908684 |

## Numerical scope

The (N_up,N_down)=(7,8) sector has 147,232,800 amplitudes. The 97 records comprise 24 continuous time points, 72 ideal circuit points and one independent direct-time check. Maximum norm error is 2.9e-14; the accumulated analytical Chebyshev state-tail bound is 2.9e-13. The independent direct evolution to Jt=1 agrees with incremental evolution within 1.4e-14 in the reported observables. The tail bound excludes floating-point rounding; independent small-state comparisons and this direct-time check control the numerical implementation.

This finite-system calculation addresses the actual experimental input and ordering, which differ from the half-filled Néel/plaquette resource baseline. It does not certify every point at 0.01 or 0.001, determine a universal mitigation exponent, or validate larger interacting lattices.
