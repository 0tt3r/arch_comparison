# Scientific evidence

Versioned result directories retain the inputs, compiler artifacts, counts,
model outputs, and source hashes used at that stage:

| Directory | Content |
|---|---|
| `published_benchmarks_v1` | TopQAD Table IV area/time replay, calibrated growth transfer, reduced Phasecraft raw shots and processed observables, source availability and confidence table |
| `integrated_v1` | Current component-informed STAR/hybrid and synthesis-selected FT capacity; explicit NISQ offset/variance requirements and training sensitivities |
| `frontier_v1` | Same-Hamiltonian QSVT costs, 24M growth shots at distances 13/25, geometric magic delivery, independent exposure diagnostics, centered locality certificate |
| `implementation_v1` | Held-out NISQ bias/variance, STAR projected channels and decoded growth probes, explicit CUDA-Q FT arithmetic, nested finite-patch checks |
| `physics_sweep_v1` | Frozen-baseline p/precision sweep, source/injection scenarios, Hamming-weight and compact-encoding costs, small-patch ordering experiment |
| `workshop_reconciliation_v1` | Shared-task adapters, native count/delivery differences, identical STAR ledger, captured workshop Python sources |
| `analytical_v1` | Initial analytic curves and sensitivities |
| `cudaq_v1` | Analog programs, native routes, synthesized words, backend diagnostics |
| `nisq_v2` | Native Rzz compilation, exact cones, parity acceptance, factory sensitivities |
| `refinement_v3` | Parallel compilation, patch dynamics, noisy observables, cultivation, STAR retries |
| `presentation_v1` | Capacity inversion and the repository/presentation provenance record |
| `discussion_v1` | Archived fixed-step analytical example and its original figures (source at commit `7212660`) |
| `discussion_v2` | Archived cubic-time, fixed-allocation example and its figures (source at commit `26b9a07`) |
| `delivery_v1` | Reoptimized data-block, analog-lane, and shared-source delivery sensitivities |
| `dynamics_4x4` | Matrix-free exact reference, directly measured Trotter steps, PBS provenance, and finite-lattice resource sensitivity |
| `dynamics_followup_v1` | Direct checks of saved headline allocations, 4 × 4 evolution at tau = 2,4, and isolated patch-calibration capacity sensitivity |
| `discussion_v3` | Measured time dependence by ordering, joint budget/source optimization, finite-patch and grid checks |

The principal current figures are in [../figures](../figures), and the full
story is in [../tex/main.pdf](../tex/main.pdf). Presentation changes do not
turn conditional circuit costs into accuracy-qualified dynamics results.

The v3 manifest remains unchanged by the repository reorganization. Its source
hashes resolve to original snapshots where the current source has moved or
changed; the audit also replays the calculations using current code. Historical
root filenames map to their new locations in
[../configs/file_locations.json](../configs/file_locations.json).
