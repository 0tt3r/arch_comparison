# Scientific evidence

Versioned result directories retain the inputs, compiler artifacts, counts,
model outputs, and source hashes used at that stage:

| Directory | Content |
|---|---|
| `analytical_v1` | Initial analytic curves and sensitivities |
| `cudaq_v1` | Analog programs, native routes, synthesized words, backend diagnostics |
| `nisq_v2` | Native Rzz compilation, exact cones, parity acceptance, factory sensitivities |
| `refinement_v3` | Parallel compilation, patch dynamics, noisy observables, cultivation, STAR retries |
| `presentation_v1` | Capacity inversion and the repository/presentation provenance record |
| `discussion_v1` | Archived fixed-step analytical example and its original figures (source at commit `7212660`) |
| `discussion_v2` | Archived cubic-time, fixed-allocation example and its figures (source at commit `26b9a07`) |
| `discussion_v3` | Measured time dependence by ordering, joint budget/source optimization, finite-patch and grid checks |

The principal current figures are in [../figures](../figures), and the full
story is in [../tex/main.pdf](../tex/main.pdf). Presentation changes do not
turn conditional circuit costs into accuracy-qualified dynamics results.

The v3 manifest remains unchanged by the repository reorganization. Its source
hashes resolve to original snapshots where the current source has moved or
changed; the audit also replays the calculations using current code. Historical
root filenames map to their new locations in
[../configs/file_locations.json](../configs/file_locations.json).
