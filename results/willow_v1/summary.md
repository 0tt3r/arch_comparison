# Willow: measured operating points and what the sampling budget means

[Operating-point figure](../../figures/willow_operating_points.pdf) · [Sampling figure](../../figures/willow_sampling_design.pdf) · [Methods](../../docs/WILLOW_BENCHMARK.md)

This pass compares the public experiment with the questions in our resource model. It does not fit a universal mitigation exponent or change the integrated capacity curves.

The reduction covers 36 declared center-hole workloads: 4×4 and 6×6, both fluxes, U=0,4,8, and one to three steps. It includes 864 time points and 172,800,000 attempted shots. Each time point has 20 randomized compilations with 10,000 shots each. A U=0 training time series costs additional shots, and the published smoothing uses shared time-series information.

## Resource scale

| L | Steps | Physical qubits | Two-qubit gates | Strict acceptance, median | Hypothetical two-qubit PEC variance factor, p=0.001 |
|---:|---:|---:|---:|---:|---:|
| 4 | 1 | 34 | 420 | 0.250 | 5.37 |
| 4 | 2 | 34 | 820 | 0.126 | 26.6 |
| 4 | 3 | 34 | 1220 | 0.082 | 132 |
| 6 | 1 | 74 | 1492 | 0.037 | 391 |
| 6 | 2 | 74 | 2932 | 0.023 | 1.24e+05 |
| 6 | 3 | 74 | 4372 | 0.020 | 3.94e+07 |

The paper counts 34 and 74 physical qubits, including ancillas; the models have 32 and 72 fermionic modes. Each 24-time series uses 4.8 million shots. The PEC column is a specified inverse-channel calculation at our reference p, not a prediction of the experimental pipeline. The experiment uses different hardware noise and a trained estimator, and its useful accuracy must be assessed separately.

## Accuracy scale, without promoting approximate references to truth

The table below uses three steps and zero flux, all 24 nonzero times and all sites. Derived doublon uncertainties are covariance-agnostic upper bounds from the supplied Z/ZZ errors; they are not measured total errors.

| L | U | Reference | Local doublon RMS difference | Mean doublon RMS difference | Local Z median reported SE |
|---:|---:|---|---:|---:|---:|
| 4 | 0 | continuous FLO | 0.0050 | 0.0027 | 0.0152 |
| 4 | 4 | TDVP chi=2048 | 0.0096 | 0.0068 | 0.0156 |
| 4 | 8 | TDVP chi=2048 | 0.0135 | 0.0094 | 0.0160 |
| 6 | 0 | continuous FLO | 0.0183 | 0.0143 | 0.0803 |
| 6 | 4 | TDVP chi=2048 | 0.0479 | 0.0451 | 0.0810 |
| 6 | 8 | TDVP chi=2048 | 0.0706 | 0.0668 | 0.0826 |

U=0 is a training consistency check. The archived `exact` values change with Trotter depth and describe the free circuit. We independently compute continuous U=0 dynamics by orbital evolution, checked against the existing many-body sector solver. At U=4 and 8, the exact groups are empty: the exact interacting comparisons shown in the paper are not supplied in this pinned archive. The MPS comparisons therefore include reference and Trotter errors, not just experimental error.

## Sampling design matters before an exponential cost law

Shots from one randomized circuit share its noise realization. With K circuit groups, n retained shots per group and conditional observable variance v, a useful planning decomposition is Var(mean) ≈ (a + v/n)/K, where a is variation between randomized-circuit expectations. Increasing n at fixed K does not remove the a/K contribution; redrawing independent circuits can. The fitted nonnegative a is a diagnostic from only 20 groups, not a certified irreducible floor.

| L | Central doublon SE now | Same circuits, 100× shots | 100× independent circuits |
|---:|---:|---:|---:|
| 4 | 0.0043 | 0.0032 | 0.00043 |
| 6 | 0.0022 | 0.0010 | 0.00022 |

These examples use U=4, zero flux, three steps and Jt=1, before TFLO and smoothing. They account for number-selection rejection (exact sector at 4×4; at most two total occupation-number deviations at 6×6). They forecast sampling variance only, assuming new circuit groups have the same distribution and are independent; they do not forecast mitigation accuracy.

## Archive and pipeline checks

Raw published Z and onsite ZZ-derived doublons reproduce from equal-weight circuit groups, number selection and diagonal symmetry in 36 of 36 cases using the declared strict/relaxed/no-selection menu. The recovered number-deviation thresholds are zero at 4×4, two at 6×6 with zero flux, and one at 6×6 with pi flux. Thus the archive label `raw` already includes processing; it must not be compared to unfiltered bits as though they were the same estimator.

The initial-state JSON repeats a singlet configuration under the center-hole labels, including out-of-range coordinates for 4×4. It is retained and flagged, not used. The consistent product-state occupations are read from the zero-time TDVP observables and checked for one missing particle. All selected ZIP members have checked sizes/CRCs and SHA-256 hashes. About 1.06 GB was fetched by byte ranges; the full 8.8 GB archive MD5 was not verified.

## Consequence for the architecture comparison

The physical-qubit experiment is a real low-footprint operating point with a shared, trained estimator and observable-dependent accuracy. Our gate-count scale can describe its circuit size, but a gatewise PEC budget cannot describe its estimator. Neither the experimental curves nor their reference differences certify every local observable at 0.01 or 0.001. The useful next accuracy check is the missing matched interacting circuit/continuous reference for this pipeline, not a larger arbitrary held-out training campaign. For STAR/FT, full injection/growth and delivery remain separate unresolved physical inputs.
