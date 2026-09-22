# Back-of-the-envelope+ discussion example

Central scenarios, one week, p = 0.001, 10 ns gates, 240 ns syndrome cycles.
These are analytical planning estimates; mitigation and Trotter biases are assumed affordable.
NISQ gate/depth counts are anchored to Willow. Its assumed sampling penalty is 100 at the anchor, not a measured fit.

## Largest lattice at tau = 1, absolute error 0.001

| Physical qubits | NISQ | STAR | FTQC |
|---:|---:|---:|---:|
| 100 | 6 × 6 | — | — |
| 10,000 | 8 × 8 | — | — |
| 100,000 | 8 × 8 | 4 × 4 | 4 × 4 |
| 1,000,000 | 8 × 8 | 8 × 8 | 4 × 4 |
| 10,000,000 | 8 × 8 | 10 × 10 | 8 × 8 |

## Area required for 6 × 6 at tau = 1

| Absolute precision | NISQ | STAR | FTQC |
|---:|---:|---:|---:|
| 0.001 | 81 | 124,362 | 2,263,170 |
| 0.01 | 81 | 101,802 | 126,268 |

The strict FT case uses an extra magic-state cleanup stage; the loose case does not.
At fixed circuit size, sampling scales as inverse precision squared. Integer replicas and code distance add thresholds.

NISQ: exp(kappa p G2) effective sampling cost; kappa is fixed by the chosen anchor overhead.
STAR: exp(4 alpha p Theta + 4 f R), with four effective concurrent rotation lanes.
FTQC: synthesized T count, with preparation and consumption charged together.

Bands vary NISQ anchor overhead 10–1000, STAR floor 0–1e-4, delivery lanes 8–1, and factory service time 0.5–4 times central.
They are scenarios, not confidence intervals.

Details: [discussion guide](../../docs/DISCUSSION_MODEL.md).
Figures: [strict precision](../../figures/discussion_capacity.pdf), [0.01 sensitivity](../../figures/discussion_capacity_loose.pdf).
