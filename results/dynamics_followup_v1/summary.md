# Direct allocation checks and longer-time 4 × 4 dynamics

All 16 original headline allocations pass their saved Trotter reserves, before reoptimization.
The reference precision remains epsilon=0.001 (1% of D*=0.1), with 0.01 as the loose sensitivity.
This checks finite-lattice Trotter error only. NISQ mitigation bias and experimental-ordering transfer remain assumptions.

| Error | Arm | Layout | Saved r | Error / saved reserve | Saved Q | Q with measured-error budget at saved r | Best measured r | Best measured Q |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 0.001 | nisq | effective | 21 | 0.8997 | 72 | 72 | 21 | 72 |
| 0.001 | nisq | pec_two_qubit | 19 | 0.9000 | 10^13.46 | 10^12.92 | 19 | 10^12.92 |
| 0.001 | nisq | compiled_pec | 19 | 0.9000 | 10^40.84 | 10^40.30 | 19 | 10^40.30 |
| 0.001 | star | reference | 33 | 0.8933 | 71442 | 71442 | 32 | 71442 |
| 0.001 | hybrid | reference | 41 | 0.8985 | 73348 | 73348 | 41 | 73348 |
| 0.001 | ft | reference | 33 | 0.8987 | 1.8759e+06 | 1.6552e+06 | 32 | 1.6552e+06 |
| 0.001 | ft | fast_b1 | 29 | 0.8989 | 1.1555e+06 | 9.8811e+05 | 32 | 9.8811e+05 |
| 0.001 | ft | fast_four_ports | 33 | 0.8987 | 3.8505e+05 | 3.8505e+05 | 32 | 3.8505e+05 |
| 0.01 | nisq | effective | 8 | 0.9091 | 36 | 36 | 8 | 36 |
| 0.01 | nisq | pec_two_qubit | 8 | 0.9091 | 36 | 36 | 8 | 36 |
| 0.01 | nisq | compiled_pec | 8 | 0.9091 | 10^13.36 | 10^13.28 | 8 | 10^13.28 |
| 0.01 | star | reference | 11 | 0.9012 | 46818 | 46818 | 11 | 46818 |
| 0.01 | hybrid | reference | 10 | 0.9051 | 48420 | 48420 | 10 | 48420 |
| 0.01 | ft | reference | 8 | 0.9091 | 61974 | 61974 | 8 | 61974 |
| 0.01 | ft | fast_b1 | 10 | 0.9051 | 79906 | 79906 | 10 | 79906 |
| 0.01 | ft | fast_four_ports | 10 | 0.9051 | 92902 | 92902 | 10 | 92902 |

The original Q column is validated with its unchanged allocation. Subsequent columns reoptimize the allocation using measured error plus 1e-10. The last column searches only the measured step menu; it is not a global optimum.

| Patch | Time | Ordering | Exact central doublon | Tail coefficient | Tail spread | Largest sampled coefficient |
|---|---:|---|---:|---:|---:|---:|
| 4x3 | 1 | plaquette | 0.1303479408 | 0.245411 | 0.057% | 0.248182 |
| 4x3 | 1 | termwise | 0.1303479408 | 0.139358 | 0.093% | 0.142035 |
| 4x4 | 1 | plaquette | 0.1452313819 | 0.261816 | 0.025% | 0.26484 |
| 4x4 | 1 | termwise | 0.1452313819 | 0.187732 | 0.007% | 0.191033 |
| 4x4 | 2 | plaquette | 0.1436930909 | 1.0959 | 0.253% | 1.10688 |
| 4x4 | 2 | termwise | 0.1436930909 | 0.775179 | 0.225% | 0.782097 |
| 4x4 | 4 | plaquette | 0.1454409046 | 4.19314 | 0.290% | 4.24238 |
| 4x4 | 4 | termwise | 0.1454409046 | 3.02911 | 0.353% | 3.07199 |

## Capacity sensitivity

The new figure holds delivery and source choices fixed while changing the accuracy coefficient. Each curve reoptimizes error allocations and steps. The added-patch envelope retains the old 3x3 running maximum and adds every measured 4x3/4x4 coefficient through the requested time, with a 10% margin. L=2 keeps its own calibration.
This envelope is a planning sensitivity, not a bulk bound, a continuously sampled time maximum, or a validation at every selected r. The main discussion-v3 curves are preserved, and no 4x4 inference at tau=8 is made.

[Capacity figure](../../figures/patch_capacity.pdf) · [Patch/time figure](../../figures/patch_time_extension.pdf)

| Time | Ordering | Old planning coefficient | Added-patch coefficient |
|---:|---|---:|---:|
| 1 | plaquette | 0.24282 | 0.291324 |
| 1 | termwise | 0.0869206 | 0.210136 |
| 2 | plaquette | 1.06719 | 1.25359 |
| 2 | termwise | 0.424585 | 0.860307 |
| 4 | plaquette | 5.79728 | 5.79728 |
| 4 | termwise | 2.6591 | 3.37919 |
