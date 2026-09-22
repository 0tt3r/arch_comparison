# Optimized discussion example

One week, p = 0.001; measured finite-patch time dependence and per-arm error allocation. The total precision is unchanged.
The planning coefficient is 1.1 times the running maximum of converged 3 × 3 coefficients for each ordering through the requested time. This avoids relying on isolated finite-patch cancellations. Pointwise coefficients are a sensitivity; transfer to larger lattices remains an assumption. No global time power is fitted. L = 2 uses its own measured coefficients.
The coefficient also covers measured finite-step corrections at the minimum allowed r and twice r; selected optimized counts are validated directly.
Statistics receive the budget remaining after Trotter, logical, synthesis, magic, and residual-bias costs. Mitigated arms reserve 10% for uncalibrated residual bias. Its affordability is assumed.
Optimization searches the declared allocation grid, integer-step neighbors, and the shared cultivated-state menu with zero or one cleanup stage. It is not a proof of the global optimum.

## Central areas at tau = 1

| Lattice | Precision | Arm | Steps | Q | Cleanup | Statistical share | Trotter share | Magic share |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| 4 × 4 | 0.001 | nisq | 20 | 36 | 0 | 0.293 | 0.607 | 0.000 |
| 4 × 4 | 0.001 | star | 21 | 71,442 | 0 | 0.653 | 0.197 | 0.000 |
| 4 × 4 | 0.001 | hybrid | 37 | 73,348 | 0 | 0.610 | 0.177 | 0.013 |
| 4 × 4 | 0.001 | ft | 31 | 1,434,498 | 0 | 0.336 | 0.253 | 0.362 |
| 4 × 4 | 0.01 | nisq | 8 | 36 | 0 | 0.521 | 0.379 | 0.000 |
| 4 × 4 | 0.01 | star | 8 | 46,818 | 0 | 0.564 | 0.136 | 0.000 |
| 4 × 4 | 0.01 | hybrid | 10 | 48,420 | 0 | 0.457 | 0.243 | 0.000 |
| 4 × 4 | 0.01 | ft | 8 | 61,974 | 0 | 0.264 | 0.379 | 0.007 |
| 6 × 6 | 0.001 | nisq | 18 | ~10^10.6 | 0 | 0.151 | 0.749 | 0.000 |
| 6 × 6 | 0.001 | star | 20 | 298,356 | 0 | 0.633 | 0.217 | 0.000 |
| 6 × 6 | 0.001 | hybrid | 34 | 180,798 | 0 | 0.649 | 0.210 | 0.036 |
| 6 × 6 | 0.001 | ft | 36 | 11,554,972 | 1 | 0.733 | 0.187 | 0.000 |
| 6 × 6 | 0.01 | nisq | 8 | 81 | 0 | 0.521 | 0.379 | 0.000 |
| 6 × 6 | 0.01 | star | 8 | 101,802 | 0 | 0.664 | 0.136 | 0.000 |
| 6 × 6 | 0.01 | hybrid | 10 | 103,548 | 0 | 0.556 | 0.243 | 0.001 |
| 6 × 6 | 0.01 | ft | 8 | 126,268 | 0 | 0.372 | 0.379 | 0.019 |

## Largest lattice at tau = 1, strict precision

| Q | NISQ effective | NISQ 2q PEC | NISQ compiled PEC | STAR | Hybrid | FTQC |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 4 × 4 | 2 × 2 | 2 × 2 | — | — | — |
| 1,000 | 4 × 4 | 2 × 2 | 2 × 2 | — | — | — |
| 10,000 | 4 × 4 | 2 × 2 | 2 × 2 | — | — | — |
| 100,000 | 4 × 4 | 2 × 2 | 2 × 2 | 4 × 4 | 4 × 4 | 2 × 2 |
| 1,000,000 | 4 × 4 | 2 × 2 | 2 × 2 | 6 × 6 | 6 × 6 | 2 × 2 |
| 10,000,000 | 4 × 4 | 2 × 2 | 2 × 2 | 8 × 8 | 8 × 8 | 4 × 4 |
| 100,000,000 | 4 × 4 | 2 × 2 | 2 × 2 | 10 × 10 | 12 × 12 | 10 × 10 |

The [ablation table](ablation.csv) separates budget optimization under the old coefficient from replacing the coefficient. The [budget table](budgets.csv) records every selected split.
The time-growing panel stops at L = 32 because tau = L/4 then reaches the measured endpoint tau = 8. The fixed-time sweeps still reach L = 48.
STAR uses the termwise ordering; FT/hybrid use plaquettes. The NISQ gate proxy transfers the plaquette coefficient, without claiming a match to the experimental ordering.
FT and hybrid use the same source quality, service, and cleanup menu. Source-quality/throughput tradeoffs beyond this menu remain uncalibrated.
PEC references retain Hoeffding counts. Other curves use Gaussian planning counts. Bands show delivery/mitigation scenarios, not confidence intervals or bulk-calibration uncertainty.

Figures: [strict](../../figures/discussion_capacity.pdf), [loose](../../figures/discussion_capacity_loose.pdf), [time calibration](../../figures/time_calibration.pdf).
Full story: [main manuscript](../../tex/main.pdf). Previous [fixed-allocation example](../discussion_v2/summary.md).
