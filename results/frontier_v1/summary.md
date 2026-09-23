# Same-task algorithms and measured components

This pass keeps the previous capacity baseline frozen. New costs are separately identified algorithm and layout scenarios; new channel measurements constrain individual mechanisms.

## Qubitized dynamics on the same open-boundary Hamiltonian

The centered Hamiltonian has alpha = 4L(L−1) + L² at U/J=4: alpha=64 on 4×4. Repeating each onsite Pauli twice makes all coefficients magnitude 1/2. At 4×4 the 128-term PREPARE is seven Hadamards. A controlled block encoding uses 508 T and no arbitrary rotations, versus 508 T plus 256 rotations with generic amplitude preparation. The local adapter repairs SELECT signs/zero padding without modifying pyLIQTR.

The cosine/sine branches and their LCU give accepted amplitude 1/4 before amplification. Counts optimize zero through three amplification steps, simulation/synthesis/logical budgets, source count, and cleanup. Conditional fault bias and rejected attempts are charged. Large-degree phase generation and physical SELECT routing remain unvalidated. The two timing choices bracket extra Clifford serialization only; they are not rigorous hardware bounds.

| 4×4, tau=1 | Absolute error | T per attempted circuit | Q for one week |
|---|---:|---:|---:|
| split_uniform, Clifford serialization=0 | 0.01 | 296,556 | 159,622 |
| split_uniform, Clifford serialization=1 | 0.01 | 296,556 | 181,024 |
| trotter_individual | 0.01 | 45,984 | 92,902 |
| trotter_hwp | 0.01 | 17,885 | 113,236 |
| split_uniform, Clifford serialization=0 | 0.001 | 551,403 | 9,506,556 |
| split_uniform, Clifford serialization=1 | 0.001 | 551,403 | 12,147,266 |
| trotter_individual | 0.001 | 192,332 | 1,335,536 |
| trotter_hwp | 0.001 | 82,336 | 177,472 |

QSVT beats the best of these two Trotter scenarios at 0/32 tested size/time/precision points. This does not rule out better structured block encodings or coherent observable estimation. At tau=8 the Trotter coefficient still comes only from 3×3; its finite-size transfer is explicit.

## Exposure diagnostics

The measured NISQ variance multipliers are much smaller than whole-circuit PEC inverse-norm costs, but many held-out bias tests fail. Passing points and failures remain on the same axes. Training adds six million shots per fit and is retained in the table. These data do not collapse to a universal variance-plus-bias law in one exposure coordinate. STAR q/(p theta) is plotted from the enumerated generator-fault channels, with ideal and noisy PCEC separated; it is a submodel test, not full injection calibration.

## Larger decoded growth components

PBS smoke 7790942 and production 7790943 both completed F/0. The production census has 24 cases of one million shots (24 million total), all with pinned seeds, circuit hashes, ideal controls, and simultaneous one-sided 95% upper bounds over all 24 components.

| Initial→final distance | p | X failures | Z failures | Pauli-union upper bound |
|---|---:|---:|---:|---:|
| 5→13 | 0.001 | 60 | 175 | 0.000302 |
| 5→13 | 0.0001 | 0 | 0 | 1.23e-05 |
| 9→13 | 0.001 | 0 | 0 | 1.23e-05 |
| 9→13 | 0.0001 | 0 | 0 | 1.23e-05 |
| 13→13 | 0.001 | 0 | 0 | 1.23e-05 |
| 13→13 | 0.0001 | 0 | 0 | 1.23e-05 |
| 5→25 | 0.001 | 51 | 419 | 0.000556 |
| 5→25 | 0.0001 | 0 | 0 | 1.23e-05 |
| 9→25 | 0.001 | 0 | 0 | 1.23e-05 |
| 9→25 | 0.0001 | 0 | 0 | 1.23e-05 |
| 13→25 | 0.001 | 0 | 0 | 1.23e-05 |
| 13→25 | 0.0001 | 0 | 0 | 1.23e-05 |

The X/Z union assumes both experiments describe one Pauli logical channel. For zero-failure pairs its bound is 1.235e-5, so this shot count does not certify a 1e-5 total floor. The distance-5 seed is visibly worse at p=1e-3. The starting state is ideally encoded; growth, target memory, and readout are included. No subtraction of the memory control or claim about a full injected-state channel is made.

## Explicit magic delivery

Magic states enter through perimeter ports, walk by nearest-neighbor SWAPs, are consumed, and return the data permutation by inverse SWAPs. Traversed data patches and surgery mediators remain reserved until return. One bounded output buffer per port and the entire square footprint (225 or 289 tiles) are charged, with the reference source area added. This is a valid conservative logical construction, not an optimized factory layout or physical surgery simulation.

| Batch | Ports | Abstract delivery (ms) | Explicit delivery (ms) | Array tiles |
|---:|---:|---:|---:|---:|
| 1 | 1 | 46.824 | 581.118 | 225 |
| 1 | 4 | 12.282 | 100.518 | 225 |
| 4 | 1 | 43.518 | 691.434 | 225 |
| 4 | 4 | 17.292 | 145.830 | 225 |
| 8 | 1 | 33.222 | 525.618 | 289 |
| 8 | 4 | 15.666 | 149.454 | 289 |
| 16 | 1 | 27.510 | 358.890 | 289 |
| 16 | 4 | 16.398 | 84.102 | 289 |

The batch ranking changes with placement, transport, and port count. An aggregate HWP T-count saving therefore remains conditional until the complete chosen layout is scheduled. The previous capacity figure is not silently replaced by this isolated onsite stage.

## Centered locality certificate

After removing onsite U by an interaction picture, each two-spin hopping bond has norm 2J. The centered doublon has norm 1/2. Counting four first bonds and at most seven successors gives a sufficient omitted-boundary bound min(1, (2/7) sum[n≥R+1] (28 tau)^n/n!). At tau=1, the geometric majorant first meets 0.0025 at R=78 (a 157×157 enclosing square). This very loose certificate cannot justify replacing the full lattice by a 6×6 or 8×8 patch. It also does not say the physical light cone has radius 78. Finite-reference differences remain the useful necessary check, with outer-boundary error unbounded.

## What is settled and what remains scientific work

All six calculation deliverables in docs/FRONTIER_WORK.md are implemented and checked. The remaining evidence gaps are end-to-end noisy arbitrary-angle injection/growth, transfer of trained mitigation beyond eight qubits, tight state-dependent bulk truncation, and full routed SELECT/HWP schedules. These require additional protocol or physics input; passing this software audit does not settle them.
