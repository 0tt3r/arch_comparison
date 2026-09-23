# Next comparison layer

The accepted objective remains a fair resource comparison and supported
improvements, not agreement with the other workshop implementations. This pass
closes the following concrete calculation tasks and reports what they do and
do not settle:

1. Count qubitized dynamics for the same open-boundary, half-filled Hamiltonian;
   include polynomial approximation, synthesis, state-selection/amplification,
   and the existing shot and source accounting. Keep periodic specialized
   encodings and coherent estimators separately identifiable.
2. Test the proposed exposure variables against the independently computed
   NISQ bias/variance and STAR channels. Failed bias targets must stay visible.
3. Add geometric magic-state delivery and storage/surgery area to the explicit
   FT stage. Compare source rates using one layout contract and report routing
   cost rather than assigning free nonlocal delivery.
4. Extend decoded growth components to distances 13 and 25, with fixed shot
   counts, simultaneous upper bounds, complete PBS accounting, and no conversion
   of zero failures into zero error. These remain components of injection.
5. Construct a centered, state-independent locality certificate and compare it
   with finite-patch evidence. If the sufficient patch is impractically large,
   report that result rather than interpreting a loose bound as physical range.
6. Assemble these results in the main manuscript, with reproducible figures,
   independent tests and source checks; commit the completed work.

An end-to-end hardware-validated arbitrary-angle source, a universal bound on
trained NISQ bias, and a classical solution of bulk interacting dynamics are
not finite software deliverables. The task is to quantify their effect on the
comparison and avoid presenting them as established by smaller checks.

## Completed calculation layer

| Deliverable | Implementation and evidence | Claim boundary |
|---|---|---|
| Same-task qubitization | `python/qubitization_study.py`; 8 compiled inventories and 128 cost rows in `results/frontier_v1/qubitization/` | Polynomial/query construction and shared source ledger; large-degree phases and routed SELECT unvalidated |
| Independent exposure diagnostics | `python/exposure_diagnostics.py`; 384 held-out estimators and 1,640 generator-fault channel points | Eight-qubit density matrices and restricted injection channels; no universal mitigation law |
| Geometric magic delivery | `python/geometric_delivery.py`; 8 schedules with complete resource reservations | Constructive SWAP transport for one onsite stage, not an optimized full algorithm |
| Data-distance growth components | `python/frontier_growth.py`; 24 million fixed shots, 24 cells, complete PBS accounting | Ideal encoded input; growth, target memory, and readout included; full injection excluded |
| Centered locality certificate | `python/locality_certificate.py`; analytic commutator tail and sufficient-radius table | Deliberately loose state-independent bound, not measured propagation speed |
| Synthesis and verification | `python/frontier_analysis.py`, `tests/test_frontier.py`, `tex/main.tex` | Numerical/source/circuit controls, not hardware validation |

The [result summary](../results/frontier_v1/summary.md) reports the numbers and
remaining scientific dependencies. The baseline datasets and supplied workshop
folders are preserved.

## Reproduction

From the repository root, with the `qb_gsee` interpreter and single-threaded
BLAS:

```bash
make PYTHON=/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python frontier
make PYTHON=/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python frontier-check
```

`frontier` regenerates the local calculations using the completed growth data.
It does not submit jobs. `frontier-figures` regenerates figures and the summary
without recounting pyLIQTR gates. The package versions are saved with the counts.
Uniform-state preparation is independently simulated in the science tests; the
Hamiltonian, Bessel approximation, amplification normalization, conditional
budgets, magic-transport channel, and resource reservations have separate checks.

The pyLIQTR adapter is local to this repository: SELECT receives only signed
unitaries and identity padding; PREPARE carries nonnegative magnitudes. The
split-weight encoding uses Qualtran's uniform preparation. No installed package
or source checkout was modified.

## Growth campaign contract

The frozen inputs are listed in
`results/frontier_v1/growth_inputs.sha256`. Production uses one Improv node and
eight independent workers, each with bounded detector batches. Every case has
a source/circuit contract, a fixed seed, exclusive writer lock, and atomic final
record. Sampling stops at one million shots even if failures are absent.
Confidence intervals are Clopper–Pearson; the one-sided upper bounds spend
0.05/24 per component. The sum of X/Z upper bounds is conditional on a common
Pauli-channel interpretation and is not a direct bound on arbitrary-angle
injection.

The local distance-25 topology check and PBS smoke preceded production. The
normal `QC_NA_Sim` allocation had no hours remaining, so the completed jobs used
the permitted backfill queue:

```bash
qsub -r y -q backfill -v FRONTIER_MODE=smoke \
  -o results/frontier_v1/growth_smoke.pbs.log pbs/frontier_growth.pbs
qsub -r y -q backfill -l walltime=02:00:00 -v FRONTIER_MODE=production \
  -o results/frontier_v1/growth_production.pbs.log pbs/frontier_growth.pbs
```

Smoke `7790942.imgt1` completed F/0 in 40 seconds; production `7790943.imgt1`
completed F/0 in 4 minutes 53 seconds. Both final `qstat -xf` records and logs
are saved. Validation requires every expected case and its full shot count,
not only the scheduler exit code. The source hashes are unchanged after the
production run. Re-running the original contract resumes completed cases; for
an independent campaign or changed inputs, use a new output directory and
contract rather than overwriting these records. The archived audit intentionally
checks the canonical job identifiers and accounting.

## Interpretation

The qubitization timing comparison uses the same fast-block and source menu
as its Trotter comparators. Its no-extra-Clifford and serial-Clifford options
are sensitivities, not lower and upper bounds on all possible architectures.
Hamming-weight costs in that comparison remain aggregate scenarios; the new
geometric onsite stage shows why their full routing still matters.

The locality certificate uses graph distance around the observed site, not a
corner-anchored square. Consistent initial reduced states and even fermionic
observables are required. The large sufficient radius is a limitation of this
particular uniform bound. It cannot be used as a necessary patch size or an
empirical bulk-convergence estimate.
