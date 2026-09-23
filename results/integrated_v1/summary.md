# Integrated comparison and NISQ accuracy requirements

This is the current integrated comparison. Historical sweeps remain available for replay. The encoded bands are scenario envelopes; they do not establish statistical error bars or an unconditional architecture ordering.

[Capacity figure](../../figures/integrated_capacity.pdf) · [NISQ requirements](../../figures/nisq_requirements.pdf) · [Methods](../../docs/INTEGRATED_COMPARISON.md)

## Capacity at tau=1 and one million qubits

Entries are central lattice sides with the full scenario range in brackets. The shared accuracy-transfer and delivery assumptions still apply.

| p | Error | STAR | Hybrid | FTQC | Unique leader across scenarios |
|---:|---:|---:|---:|---:|---|
| 0.001 | 0.01 | 10 [10, 10] | 12 [12, 12] | 12 [12, 12] | unresolved |
| 0.001 | 0.001 | 4 [4, 6] | 6 [6, 6] | 4 [4, 4] | unresolved |
| 0.0001 | 0.01 | 16 [16, 24] | 24 [20, 24] | 24 [24, 24] | unresolved |
| 0.0001 | 0.001 | 8 [8, 12] | 8 [8, 12] | 6 [6, 6] | unresolved |

The table leaves close rankings unresolved. The full capacity table also identifies regions with a unique leader across this menu; that statement remains conditional on the shared component-composition and layout assumptions.

## Current 4×4, tau=1 operating points

| p | Error | Arm | Scenario | Q for one week | Q per copy | seconds/attempt | steps | injection seed / data d | log variance multiplier |
|---:|---:|---|---|---:|---:|---:|---:|---|---:|
| 0.001 | 0.01 | star | component | 60,648 | 60,648 | 0.06362 | 10 | 9 / 19 | 0.685 |
| 0.001 | 0.01 | hybrid | component | 62,394 | 62,394 | 0.1038 | 13 | 9 / 19 | 0.4058 |
| 0.001 | 0.01 | ft | component | 75,994 | 75,994 | 2.696 | 12 | — / 21 | 0 |
| 0.001 | 0.01 | star | zero_growth | 60,648 | 60,648 | 0.06981 | 11 | 9 / 19 | 0.4755 |
| 0.001 | 0.01 | hybrid | zero_growth | 62,394 | 62,394 | 0.1038 | 13 | 9 / 19 | 0.3169 |
| 0.001 | 0.01 | ft | zero_growth | 75,994 | 75,994 | 2.696 | 12 | — / 21 | 0 |
| 0.001 | 0.01 | star | slow | 60,648 | 60,648 | 0.06733 | 10 | 9 / 19 | 0.685 |
| 0.001 | 0.01 | hybrid | slow | 62,394 | 62,394 | 0.1972 | 13 | 9 / 19 | 0.4058 |
| 0.001 | 0.01 | ft | slow | 75,994 | 75,994 | 5.553 | 10 | — / 21 | 0 |
| 0.001 | 0.001 | star | component | 88,872 | 88,872 | 0.2081 | 28 | 9 / 23 | 1.123 |
| 0.001 | 0.001 | hybrid | component | 90,954 | 90,954 | 0.3613 | 38 | 9 / 23 | 0.5938 |
| 0.001 | 0.001 | ft | component | 177,472 | 177,472 | 0.4844 | 37 | — / 23 | 0 |
| 0.001 | 0.001 | star | zero_growth | 74,088 | 74,088 | 0.2098 | 31 | 9 / 21 | 0.5158 |
| 0.001 | 0.001 | hybrid | zero_growth | 90,954 | 90,954 | 0.3894 | 41 | 9 / 23 | 0.3293 |
| 0.001 | 0.001 | ft | zero_growth | 177,472 | 177,472 | 0.4844 | 37 | — / 23 | 0 |
| 0.001 | 0.001 | star | slow | 105,000 | 105,000 | 0.2026 | 28 | 9 / 25 | 1.123 |
| 0.001 | 0.001 | hybrid | slow | 93,036 | 93,036 | 0.3657 | 38 | 9 / 23 | 0.5938 |
| 0.001 | 0.001 | ft | slow | 226,384 | 226,384 | 0.5569 | 36 | — / 25 | 0 |
| 0.0001 | 0.01 | star | component | 28,392 | 28,392 | 0.03669 | 10 | 5 / 13 | 0.2517 |
| 0.0001 | 0.01 | hybrid | component | 39,274 | 39,274 | 0.07381 | 12 | 5 / 15 | 0.1034 |
| 0.0001 | 0.01 | ft | component | 39,274 | 39,274 | 1.929 | 12 | — / 15 | 0 |
| 0.0001 | 0.01 | star | zero_growth | 28,392 | 28,392 | 0.04336 | 11 | 9 / 13 | 0.04755 |
| 0.0001 | 0.01 | hybrid | zero_growth | 39,274 | 39,274 | 0.08167 | 13 | 9 / 15 | 0.0317 |
| 0.0001 | 0.01 | ft | zero_growth | 39,274 | 39,274 | 1.929 | 12 | — / 15 | 0 |
| 0.0001 | 0.01 | star | slow | 28,392 | 28,392 | 0.03672 | 10 | 5 / 13 | 0.2517 |
| 0.0001 | 0.01 | hybrid | slow | 39,274 | 39,274 | 0.1787 | 12 | 5 / 15 | 0.1034 |
| 0.0001 | 0.01 | ft | slow | 39,274 | 39,274 | 6.115 | 12 | — / 15 | 0 |
| 0.0001 | 0.001 | star | component | 28,392 | 28,392 | 0.102 | 28 | 5 / 13 | 0.6514 |
| 0.0001 | 0.001 | hybrid | component | 39,274 | 39,274 | 0.2263 | 37 | 5 / 15 | 0.2744 |
| 0.0001 | 0.001 | ft | component | 80,192 | 80,192 | 0.3285 | 36 | — / 15 | 0 |
| 0.0001 | 0.001 | star | zero_growth | 28,392 | 28,392 | 0.1328 | 35 | 9 / 13 | 0.05132 |
| 0.0001 | 0.001 | hybrid | zero_growth | 39,274 | 39,274 | 0.2476 | 40 | 9 / 15 | 0.03311 |
| 0.0001 | 0.001 | ft | zero_growth | 80,192 | 80,192 | 0.3285 | 36 | — / 15 | 0 |
| 0.0001 | 0.001 | star | slow | 28,392 | 28,392 | 0.102 | 28 | 5 / 13 | 0.6514 |
| 0.0001 | 0.001 | hybrid | slow | 39,274 | 39,274 | 0.551 | 37 | 5 / 15 | 0.2744 |
| 0.0001 | 0.001 | ft | slow | 91,984 | 91,984 | 0.5329 | 36 | — / 15 | 0 |

STAR and hybrid use the exact existing generator-fault submodel at the workload angles, including noisy PCEC, repeat-until-success service time, detector rejection, coherent residue and retry tails. Seed 5 uses its measured growth-campaign upper bound; no claim is made about untested seed 7. The zero-growth sensitivity is available for seed 9, whose zero-failure result is an upper limit rather than a zero-error measurement. Seed 25 avoids growth but pays its preparation acceptance. Receivers use distance-13 or distance-25 component anchors. Greater data distances add protected handoff time and volume; their injection channels are not newly calibrated.

FT selects individual synthesis or Hamming-weight batches 4, 8, 16, and compact or fast one-port blocks. These are checked arithmetic counts with an aggregate delivery model, not full routed schedules. Hybrid fixed-angle supply uses the same cultivation quality and source menu. The slow scenario increases detector rejection and multiplies source preparation cycles by four.

## NISQ requirements at Q=10^6, p=10^-3, tau=1

The remaining offset is set to 10% of epsilon in this table. Training sensitivity charges six million additional attempts of the same duration; it is not a proposed experimental training budget. Comparison is to the fastest of the selected encoded implementations that fits Q, under the component scenario. Those implementations were optimized for one-week area, not reoptimized for runtime at Q.

| Error | L | Required speedup | Training shots | Encoded comparator | Largest F | Steps | kappa threshold at that step count |
|---:|---:|---:|---:|---|---:|---:|---:|
| 0.01 | 4 | 1 | 0 | star | 3.26e+07 | 13 | 3.08 |
| 0.01 | 4 | 1 | 6000000 | star | 3.26e+07 | 13 | 3.08 |
| 0.01 | 4 | 10 | 0 | star | 3.26e+06 | 13 | 2.67 |
| 0.01 | 4 | 10 | 6000000 | star | 3.26e+06 | 13 | 2.67 |
| 0.01 | 6 | 1 | 0 | hybrid | 1.38e+08 | 13 | 0.989 |
| 0.01 | 6 | 1 | 6000000 | hybrid | 1.38e+08 | 13 | 0.989 |
| 0.01 | 6 | 10 | 0 | hybrid | 1.38e+07 | 13 | 0.868 |
| 0.01 | 6 | 10 | 6000000 | hybrid | 1.38e+07 | 13 | 0.868 |
| 0.01 | 8 | 1 | 0 | hybrid | 3.68e+08 | 13 | 0.439 |
| 0.01 | 8 | 1 | 6000000 | hybrid | 3.68e+08 | 13 | 0.439 |
| 0.01 | 8 | 10 | 0 | hybrid | 3.68e+07 | 13 | 0.388 |
| 0.01 | 8 | 10 | 6000000 | hybrid | 3.68e+07 | 13 | 0.388 |
| 0.001 | 4 | 1 | 0 | hybrid | 8.08e+07 | 40 | 1.05 |
| 0.001 | 4 | 1 | 6000000 | hybrid | 8.08e+07 | 40 | 1.05 |
| 0.001 | 4 | 10 | 0 | hybrid | 8.08e+06 | 40 | 0.921 |
| 0.001 | 4 | 10 | 6000000 | hybrid | 8.08e+06 | 40 | 0.921 |
| 0.001 | 6 | 1 | 0 | hybrid | 1.26e+08 | 40 | 0.32 |
| 0.001 | 6 | 1 | 6000000 | hybrid | 1.26e+08 | 40 | 0.32 |
| 0.001 | 6 | 10 | 0 | hybrid | 1.26e+07 | 40 | 0.28 |
| 0.001 | 6 | 10 | 6000000 | hybrid | 1.26e+07 | 40 | 0.28 |
| 0.001 | 8 | 1 | 0 | hybrid | 1.19e+08 | 40 | 0.135 |
| 0.001 | 8 | 1 | 6000000 | hybrid | 1.19e+08 | 40 | 0.135 |
| 0.001 | 8 | 10 | 0 | hybrid | 6.97e+07 | 40 | 0.131 |
| 0.001 | 8 | 10 | 6000000 | hybrid | 6.97e+07 | 40 | 0.131 |

F includes postselection through variance per attempted circuit. Unknown residual bias is an axis, not a fitted constant. Smoothing, training and other information-sharing strategies must meet the displayed variance and offset requirements with their own resource accounting. A larger affordable F at greater size or stricter accuracy can arise from a slower encoded comparator; it does not establish easier mitigation at the deeper physical circuit. No large noisy-dynamics calculation or experimental pipeline fit was performed.

## Interpretation

A capacity leader is marked robust only when its smallest capacity across the declared scenarios exceeds every other encoded arm’s largest capacity. This is robustness to this menu, not a theorem about architectures. Ties and overlaps remain unresolved. NISQ has no unconditional size boundary in these results. The data-distance floor, transferred finite-lattice Trotter coefficients, growth composition, source quality and uncompiled delivery remain named assumptions.

Run `make integrated` to resume/recompute the new version and regenerate its figures. `make integrated-check` audits source contracts, component composition and accounting.
