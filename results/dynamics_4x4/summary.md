# Direct 4 × 4 calibration at tau = 1

Central coordinate (2,2); OBC Neel quench at U/J=4. All rows retain the reference signal D*=0.1 for sampling.
This is a measured-step comparison; it does not replace the larger-lattice planning envelope.
NISQ still transfers plaquette dynamics to its experimental circuit-count proxy, and residual mitigation bias remains assumed.

| Ordering | r | Exact doublon | Signed Trotter error | abs(error) r² |
|---|---:|---:|---:|---:|
| plaquette | 8 | 0.1452313819 | 0.0041381287 | 0.264840 |
| plaquette | 16 | 0.1452313819 | 0.0010251201 | 0.262431 |
| plaquette | 32 | 0.1452313819 | 0.00025569174 | 0.261828 |
| termwise | 8 | 0.1452313819 | -0.0029848895 | 0.191033 |
| termwise | 16 | 0.1452313819 | -0.00073588352 | 0.188386 |
| termwise | 32 | 0.1452313819 | -0.00018333239 | 0.187732 |

Resource allocations below optimize only over directly measured r values and the existing budget grid. A 1e-10 numerical allowance is added to the measured error.

| Error | Arm | Layout | r | Physical qubits |
|---:|---|---|---:|---:|
| 0.001 | nisq | effective | 32 | 2196 |
| 0.001 | nisq | pec_two_qubit | 32 | 5.6588e+21 |
| 0.001 | nisq | compiled_pec | 32 | 2.7745e+67 |
| 0.001 | star | reference | 32 | 71442 |
| 0.001 | star | lanes_2 | 32 | 73002 |
| 0.001 | hybrid | reference | 32 | 73348 |
| 0.001 | hybrid | lanes_2 | 32 | 62764 |
| 0.001 | hybrid | fast_b1 | 32 | 90988 |
| 0.001 | hybrid | fast_four_ports | 32 | 1.0686e+05 |
| 0.001 | ft | reference | 32 | 1.6552e+06 |
| 0.001 | ft | fast_b1 | 32 | 9.8811e+05 |
| 0.001 | ft | fast_four_ports | 32 | 3.8505e+05 |
| 0.01 | nisq | effective | 8 | 36 |
| 0.01 | nisq | pec_two_qubit | 8 | 36 |
| 0.01 | nisq | compiled_pec | 8 | 1.9231e+13 |
| 0.01 | star | reference | 8 | 46818 |
| 0.01 | star | lanes_2 | 8 | 49818 |
| 0.01 | hybrid | reference | 8 | 48420 |
| 0.01 | hybrid | lanes_2 | 8 | 41484 |
| 0.01 | hybrid | fast_b1 | 16 | 74668 |
| 0.01 | hybrid | fast_four_ports | 8 | 71986 |
| 0.01 | ft | reference | 8 | 61974 |
| 0.01 | ft | fast_b1 | 8 | 79906 |
| 0.01 | ft | fast_four_ports | 8 | 92902 |

The doubling grid can penalize an arm that strongly prefers intermediate r. The following **interpolated sensitivity** retains the original continuous-in-r error model, now with a 4x4 coefficient and 10% margin. Its selected r values are not all directly checked.

| Error | Arm | Layout | r | 4x4 planning coefficient | Physical qubits |
|---:|---|---|---:|---:|---:|
| 0.001 | nisq | effective | 21 | 0.291324 | 72 |
| 0.001 | nisq | pec_two_qubit | 19 | 0.291324 | 2.8764e+13 |
| 0.001 | nisq | compiled_pec | 19 | 0.291324 | 6.8807e+40 |
| 0.001 | star | reference | 33 | 0.210136 | 71442 |
| 0.001 | hybrid | reference | 41 | 0.291324 | 73348 |
| 0.001 | hybrid | fast_b1 | 43 | 0.291324 | 90988 |
| 0.001 | hybrid | fast_four_ports | 34 | 0.291324 | 1.0686e+05 |
| 0.001 | ft | reference | 33 | 0.291324 | 1.8759e+06 |
| 0.001 | ft | fast_b1 | 29 | 0.291324 | 1.1555e+06 |
| 0.001 | ft | fast_four_ports | 33 | 0.291324 | 3.8505e+05 |
| 0.01 | nisq | effective | 8 | 0.291324 | 36 |
| 0.01 | nisq | pec_two_qubit | 8 | 0.291324 | 36 |
| 0.01 | nisq | compiled_pec | 8 | 0.291324 | 2.2976e+13 |
| 0.01 | star | reference | 11 | 0.210136 | 46818 |
| 0.01 | hybrid | reference | 10 | 0.291324 | 48420 |
| 0.01 | hybrid | fast_b1 | 13 | 0.291324 | 74668 |
| 0.01 | hybrid | fast_four_ports | 14 | 0.291324 | 71986 |
| 0.01 | ft | reference | 8 | 0.291324 | 61974 |
| 0.01 | ft | fast_b1 | 10 | 0.291324 | 79906 |
| 0.01 | ft | fast_four_ports | 10 | 0.291324 | 92902 |
