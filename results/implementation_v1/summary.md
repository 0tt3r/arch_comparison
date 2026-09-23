# Bounded implementation checks

These experiments test estimator bias, a specified injection submodel, an explicit FT layer, and finite-patch truncation. The capacity curves remain conditional; these data are not silently substituted for their assumptions.

## NISQ: trained mitigation with independent interacting tests

Exact eight-qubit channels include native routing, depolarizing gates, idle noise and SPAM. Six U=0 circuits train each affine correction; U=2,4 circuits at other times are held out. Zero-angle onsite slots retain noise. Each training circuit receives one million attempted shots in the cost model.

| p | Test noise | Estimator | Bias ≤ 0.005 | Total error plan ≤ 0.01 |
|---|---|---|---:|---:|
| 0.001 | matched | Raw | 0/24 | 0/24 |
| 0.001 | matched | Number selection | 8/24 | 6/24 |
| 0.001 | matched | Trained affine | 10/24 | 10/24 |
| 0.001 | matched | Selection + trained affine | 12/24 | 10/24 |
| 0.001 | mismatch | Raw | 0/24 | 0/24 |
| 0.001 | mismatch | Number selection | 8/24 | 6/24 |
| 0.001 | mismatch | Trained affine | 11/24 | 9/24 |
| 0.001 | mismatch | Selection + trained affine | 9/24 | 7/24 |
| 0.0001 | matched | Raw | 10/24 | 9/24 |
| 0.0001 | matched | Number selection | 23/24 | 19/24 |
| 0.0001 | matched | Trained affine | 24/24 | 19/24 |
| 0.0001 | matched | Selection + trained affine | 24/24 | 19/24 |
| 0.0001 | mismatch | Raw | 9/24 | 8/24 |
| 0.0001 | mismatch | Number selection | 23/24 | 19/24 |
| 0.0001 | mismatch | Trained affine | 23/24 | 19/24 |
| 0.0001 | mismatch | Selection + trained affine | 24/24 | 18/24 |

The total plan adds absolute mitigation bias, absolute Trotter bias and a 0.0025 Gaussian statistical allowance; it does not exploit bias cancellation. The binomial training bootstrap checks the delta-method fit variance. It does not establish hardware coverage or transfer to larger lattices.

## STAR: channel, acceptance, and measured growth components

The generator-Pauli fault channel is exact under perfect syndrome projection. Ideal PCEC makes that channel Pauli; one level of noisy correction leaves a nonzero coherent component, which is reported separately. RUS correction angles, their preparation times, and their growth errors are included. Workload truncation tails consume less than 1.1% of epsilon. Detector rejection and high-distance growth remain declared scenarios.

Twenty circuit-level X/Z memory probes use 500,000 shots each, ideal small-patch preparation, noisy direct growth, target-distance stabilization, and final readout. At p=0.001, d=5→7 gives 42 and 41 failures: 8.4e-5 and 8.2e-5. These include memory/readout and are not the error of an injected arbitrary-angle state. At p=0.0001, zero failures in those two rows give two-sided 95% upper limits 7.38e-6, not zero error. No high-distance floor is calibrated by these probes.

## FT: compiled arithmetic and constructive nearest-neighbor scheduling

CUDA-Q Logical compiles the explicit four-T clean-AND arithmetic, X-measurement uncomputation and conditional CZs. All two-patch gates are routed to adjacent square-grid patches. For one 4×4 onsite layer, batch 16 reduces 1280 T to 460 T but increases routed CX from 32 to 1524. With one source port and one buffered state, time falls from 46.824 ms to 27.510 ms. With four ports it rises from 12.282 ms to 16.398 ms. Unlimited buffering helps time but can require hundreds of extra patches.

This is a valid logical dependency schedule under declared CX, measurement and feedforward macros. Magic-state transport paths are abstracted. It is not a physical lattice-surgery layout, decoded error calibration, optimized router, or full dynamics campaign.

## Locality: finite-reference differences

At tau=1, the 3×3 versus 4×4 doublon difference is 0.04015; 4×3 versus 4×4 is 0.01488. Every tested proper subpatch of 4×4 fails the 0.0025 allowance at tau=1 and tau=2. Similar Trotter coefficients therefore do not certify a small truncation error. These references share coordinate (1,1), retain nearby outer boundaries, and do not establish a bulk limit.

See [the manuscript](../../tex/main.pdf), [methods and reruns](../../docs/IMPLEMENTATION_CHECKS.md), and the four `implementation_*.pdf` figures.
