# Archived fixed-allocation discussion example

Superseded by [discussion v3](../discussion_v3/summary.md). Source at commit `26b9a07`; fixed error allocation and assumed cubic time law.

One week, p = 0.001, 10 ns gates, 240 ns syndrome cycles. Four architectures; effective mitigation and PEC are shown separately.
The chosen step policy is small_patch_fit; the transferred fit is 0.22 tau^3/r^2, with Trotter allowance 0.25 epsilon.
At tau = 1 this gives r = 30 for epsilon = 0.001 and r = 10 for epsilon = 0.01.
Larger sizes, other product-formula orderings, and the time law are assumptions, not calibrated accuracy.
With coefficient 0.22, a 25% allocation gives 30/10 steps (strict/loose); a 12.5% allocation gives 42/14.

## Largest lattice at tau = 1, absolute error 0.001

| Physical qubits | NISQ effective | NISQ 2q PEC | NISQ compiled PEC | STAR | Hybrid | FTQC |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 2 × 2 | 2 × 2 | 2 × 2 | — | — | — |
| 10,000 | 4 × 4 | 2 × 2 | 2 × 2 | — | — | — |
| 100,000 | 4 × 4 | 2 × 2 | 2 × 2 | 4 × 4 | 4 × 4 | 2 × 2 |
| 1,000,000 | 4 × 4 | 2 × 2 | 2 × 2 | 6 × 6 | 6 × 6 | 2 × 2 |
| 10,000,000 | 4 × 4 | 2 × 2 | 2 × 2 | 6 × 6 | 8 × 8 | 4 × 4 |

## Central area required at tau = 1

| Lattice | Absolute precision | Steps | NISQ effective | STAR | Hybrid | FTQC |
|---:|---:|---:|---:|---:|---:|---:|
| 4 × 4 | 0.001 | 30 | 1,404 | 71,442 | 73,348 | 5,664,510 |
| 4 × 4 | 0.01 | 10 | 36 | 58,482 | 60,228 | 73,348 |
| 6 × 6 | 0.001 | 30 | ~10^17.7 | 745,890 | 302,520 | 18,908,136 |
| 6 × 6 | 0.01 | 10 | 162 | 124,362 | 103,548 | 131,986 |

## What step refinement changes at 4 × 4, tau = 1, strict precision

| Arm | Q at r = 4 | Q at the chosen refined r |
|---|---:|---:|
| nisq | 36 | 1,404 |
| star | 58,482 | 71,442 |
| hybrid | 60,228 | 73,348 |
| ft | 79,066 | 5,664,510 |

The FT jump includes crossing the chosen magic-state cleanup threshold; it is not purely a linear step-count penalty.
The hybrid does not automatically beat termwise STAR, which already avoids fixed T gates.
The effective NISQ line still assumes affordable residual bias and is not a measured performance claim.
The two-qubit PEC reference uses exact inverse-channel overhead on the gate proxy. The compiled PEC edge transfers saved full-noise exposure and timing, linearly in step count.
Both PEC references use Hoeffding shot counts; central effective, STAR, hybrid, and FT curves use Gaussian planning counts.
The hybrid injects variable angles and supplies high-fidelity fixed T states. It shares the FT factory and applies the same Gamma-amplified logical/magic budget as STAR.

The fourth panel follows tau = 0.25 L. With the chosen small-patch-fit time power, r grows as L^(1.5), so the NISQ exponent grows as L^(4.5).
The L^4 law applies instead when the step count grows linearly in time; choose trotter_time_power = 2 for that sensitivity.

Details: [discussion guide](../../docs/DISCUSSION_MODEL.md).
Figures: [strict precision](figures/discussion_capacity.pdf), [0.01 sensitivity](figures/discussion_capacity_loose.pdf).
The [original fixed-step example](../discussion_v1/summary.md) is archived.
