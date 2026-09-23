# Physical error, precision, and algorithm alternatives

The p = 0.001, held-source rows reproduce the frozen added-patch baseline. The main figure leads with 0.01 absolute error; 0.001 remains equally visible.
[Capacity figure](../../figures/physical_error_capacity.pdf) · [Hamming-weight comparison](../../figures/hamming_weight_comparison.pdf) · [Full manuscript](../../tex/main.pdf)

The source footprint (1024 qubits), nominal preparation cycles (100), and minimum handoff distance (15) stay fixed. Held scenarios retain source error 2e-9 and STAR floor 1e-5. The scaled scenario assumes source error proportional to p³ and floor proportional to p. Slow supply multiplies the source preparation cycles by four. These are scenarios, not confidence bounds or new gadget calibrations.
The original STAR count law uses the exact dephasing inverse with q = 4p/15 per rotation. Its delivery model is held fixed to isolate the injection law. It is not a full reconstruction of the original protocol.
Cleanup has ideal input-error suppression 35 q³; faults in its protected circuitry and stored states are charged separately to the logical tile-time budget. No zero-fault-per-campaign requirement is imposed.

## Capacity at Q = one million, tau = 1

| p | absolute error | NISQ effective | NISQ 2q PEC | STAR angle | STAR count | Hybrid | FTQC |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.001 | 0.01 | 6 | 4 | 10 | 6 | 12 | 12 |
| 0.001 | 0.001 | 4 | 2 | 6 | 2 | 6 | 2 |
| 0.0003 | 0.01 | 10 | 6 | 16 | 12 | 20 | 20 |
| 0.0003 | 0.001 | 6 | 4 | 8 | 4 | 10 | 4 |
| 0.0001 | 0.01 | 12 | 8 | 20 | 16 | 24 | 24 |
| 0.0001 | 0.001 | 10 | 6 | 10 | 8 | 10 | 4 |

Entries are lattice sides L. A value of 48 reaches the sampled size ceiling, not a hardware limit. Full discrete staircases and ties are retained in capacity.csv.

## FT repeated-angle synthesis

| L | epsilon | Individual T / shot | Selected HWP T / shot | HWP batch | Live ancillas | Q ratio (individual / HWP) |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 0.01 | 34,752 | 34,752 | 1 | 0 | 1.00 |
| 8 | 0.01 | 182,464 | 182,464 | 1 | 0 | 1.00 |
| 16 | 0.01 | 820,688 | 443,643 | 8 | 7 | 1.16 |
| 4 | 0.001 | 192,332 | 80,123 | 16 | 15 | 6.74 |
| 8 | 0.001 | 949,112 | 155,105 | 64 | 63 | 5.78 |
| 16 | 0.001 | 4,374,668 | 708,325 | 64 | 63 | 6.04 |

HWP batches only equal-angle rotations within the same stage. A b-bit batch uses b-popcount(b) clean ancillas/Toffolis and floor(log2 b)+1 synthesized rotations. Ancillas are reused between batches but charged for the entire shot. T counts include four T per measurement-assisted Toffoli; runtime additionally charges two logical Clifford clocks per Toffoli. Routing and multiple delivery ports remain the reference layout assumptions.
At strict 4x4 neither selected circuit uses cleanup. HWP reduces the magic allocation from 38.5% to 16.0% and increases the statistical allocation from 29.8% to 58.0%. The reduced shot count makes the area gain larger than the T-count reduction.

## Local encoding and ordering

Compact encoding uses a 3N-qubit register estimate. The parity-gadget scenario charges four CNOTs per weight-three hopping Pauli and two per onsite Pauli, plus 8N preparation gates. A two-CNOT hopping scenario assumes fusion; a routing scenario doubles gate and placement costs. Native-grid layout, preparation, and 50% gate utilization are uncompiled assumptions. Both the effective-mitigation and two-qubit PEC costs are shown in algorithms.csv.
The compact scenario uses termwise evolution and the existing termwise accuracy envelope. Thus the comparison with the plaquette-calibrated native proxy changes algorithm as well as mapping; counts are separately exposed. Compact STAR keeps that same envelope and compares four with sixteen assumed analog lanes. No parallel speedup is inferred from encoding locality alone.
The separate ordering.json experiment changes only bond ordering, then changes onsite grouping as a second variant, on 2x2 and 3x3 at tau=1. These diagnostic coefficients do not replace the frozen planning envelope.

The workshop-model comparison is in [../workshop_reconciliation_v1/summary.md](../workshop_reconciliation_v1/summary.md).
Reproduce with `make physics`; audit with `python python/physics_sweep.py --check`.
