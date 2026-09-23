# Reproducing the calculations

Run commands from the repository root. Python 3.12 with NumPy, SciPy,
Matplotlib, Numba, pytest, and packaging is sufficient for the models, tests, plots,
and manuscript helper. The existing `qb_gsee` environment provides these.
Set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` for the bounded examples.

## Current integrated comparison

`make integrated` evaluates the existing STAR instrument at actual workload
angles, integrates the measured growth components as explicit scenarios, and
selects FT synthesis/layout alternatives. It also generates the NISQ accuracy
and sampling requirements. It launches no noisy-dynamics simulation or PBS job.
Channel results and completed resource cases resume under pinned inputs.
`make integrated-figures` uses saved results; `make integrated-check` audits
their contracts and accounting. See [INTEGRATED_COMPARISON.md](INTEGRATED_COMPARISON.md)
for the physical and statistical assumptions.

## Published resource estimates and experimental data

`make published-benchmarks` replays the transcribed TopQAD inputs and saved
experimental reduction, then regenerates the report and two figures.
`make published-check` verifies the numerical replay and hashes offline.
`make published-figures` only redraws the saved results. The
[methods](PUBLISHED_BENCHMARKS.md) distinguish supplied clocks, inferred growth,
and experimental reference agreement from independent validation.

To repeat the raw-data reduction, run `python python/phasecraft_data.py` first.
It needs h5py in addition to NumPy, downloads the pinned 111.9 MB archive,
checks its checksum and uses approximately 1.1 GB under ignored `build/`.
No mitigation pipeline is refitted. The 8.8 GB Willow archive is documented
as available but is not downloaded by these commands. The original archives
are not required to regenerate the committed report and figures.

The subsequent [Willow benchmark](WILLOW_BENCHMARK.md) uses
`make willow-benchmark` and `make willow-check` for offline replay.
To repeat its retrieval, run `python python/willow_archive.py`, then
`python python/willow_data.py`. Selected HTTP ranges transfer about 1.06 GB;
the extracted cache needs 2.24 GB. Member CRCs and SHA-256 hashes are checked,
while the full 8.8 GB archive MD5 is explicitly not verified. The reduction
preserves randomized-compilation groups and uses h5py. `make willow-figures`
redraws the saved results without downloading or reducing shots again.

`make willow-reference` analyzes the saved 4×4 interacting sector results and
regenerates the matched-reference figure. `make willow-reference-check` replays
the numerical accounting and free-circuit ordering controls offline. The
[methods](WILLOW_REFERENCE.md) describe the separate PBS calculation, its
initial-state and ordering contract, and its numerical checks. These Make
targets do not submit jobs or rerun the 147-million-amplitude solver.

## Same-task algorithms and larger components

See [FRONTIER_WORK.md](FRONTIER_WORK.md). `make frontier` regenerates the
open-boundary pyLIQTR/Qualtran counts, geometric delivery, exposure diagnostics,
and locality bounds, using the completed growth campaign. `make frontier-figures`
only replays saved arithmetic and plots. `make frontier-check` checks saved
accounting plus the Stim circuit hashes. The science environment additionally
needs pyLIQTR 1.4.2, Qualtran, and Cirq for the compiler/uniform-preparation test;
exact installed versions are recorded in the qubitization results. No PBS job
is submitted by these Make targets. A fresh growth campaign requires its own
output contract and the supplied PBS runner, as described in the methods.

## Estimator and implementation checks

The four new bounded studies and exact rerun commands are described in
[IMPLEMENTATION_CHECKS.md](IMPLEMENTATION_CHECKS.md). `make implementation` runs
them; `make implementation-figures` only replays and plots saved data. The
growth check uses Stim/PyMatching in `QEC_PYTHON`, and the stage compilation uses
`CUDAQ_PYTHON`. Both can be supplied as Make variables. `make audit` includes
the growth circuit controls; the ordinary science tests do not need Stim.

## Current figures and manuscript

The physical-error/precision study and bounded algorithm alternatives run with
`make physics`. This recalculates the small 2 × 2 / 3 × 3 ordering diagnostic,
864 resource scenarios, and 104 algorithm alternatives. It uses the saved
larger-patch calibration and does not launch a 4 × 4 simulation or PBS job.
Edit `configs/physics_sweep.json` for the declared source/floor scenarios.
`python python/physics_sweep.py --plot` uses saved numbers; `--check` audits
the frozen baseline replay and error budgets.

The five model/audit files changed since commit `1503a89` are preserved in
`results/physics_sweep_v1/baseline_sources.zip`. Historical audits check the
original source hashes there and replay saved costs with the current model.
They still require unchanged live simulation inputs and data; checking an
archive does not substitute for replay or rewrite the old validation records.

Run `make reconcile` for the common-point workshop comparison. It reads the
supplied folders without editing them. If they are absent, it uses the captured
Python sources in `results/workshop_reconciliation_v1/source_snapshot.zip`,
restored under `build/workshop_sources/`. The archive excludes repository
metadata, results, and environments. `--check` validates the saved arithmetic
and captured source hashes without executing the supplied packages.
Use `python python/workshop_reconciliation.py --snapshot` to reproduce from
the archive even when the original folders are present.

For the fast analytical discussion example alone, run `make discussion`.
Edit `configs/discussion_model.json` for delivery assumptions and
`configs/budget_optimization.json` for the coefficient policy and search grid.
This regenerates `results/discussion_v3/` costs and the capacity/calibration
figures from saved data, without calling CUDA-Q or a dynamics solver. See [DISCUSSION_MODEL.md](DISCUSSION_MODEL.md).

```bash
conda activate qb_gsee
make test
make paper
make audit
```

Alternatively pass the interpreter explicitly:

```bash
make PYTHON=/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python paper
```

`make paper` regenerates the capacity plot and architecture drawings, compiles
`tex/main.tex`, and records a presentation manifest. It does not rerun the
scientific simulations or change the old result directories. A Tectonic
executable on PATH is preferred; the helper also supports the existing local
conda package cache. With a standard TeX installation, `cd tex && latexmk -pdf main.tex`
builds the same main document.

## Current time calibration and optimization checks

Calibration resumes completed (patch, time, ordering) cases and writes each
completed case atomically. Use a new `--out` directory for a fresh independent
run; point the model at the intended saved calibration before comparing results.
The default 2 × 2 and 3 × 3 sweep runs on one CPU thread. The 4 × 3 checks
use about 1.4 GB and several minutes per time point on the workshop node.

```bash
python python/time_calibration.py
python python/time_calibration.py --patches 4x3 --times 2 4 --out results/discussion_v3/rectangle_checks
python python/check_discussion_v3.py --guard
make discussion
python python/check_discussion_v3.py --supplementary
python python/check_discussion_v3.py
```

`--supplementary` checks the actual error at selected optimized step counts
and searches a finer allocation grid at tau = 1. Run it again after changing
the model assumptions. `make audit` verifies saved model budgets, capacity
inversion, calibration convergence/source hashes, and the older frozen
scientific artifacts. The two old discussion directories and their figures
are archived; `make discussion` does not overwrite them.

## Latest scientific calculations

These commands intentionally regenerate the named result files. The current
presentation pass only reads those saved calculations.

```bash
.venv-cudaq/bin/python python/run_parallel_refinement.py --full-L4
.venv-cudaq/bin/python python/run_parallel_refinement.py --sizes 2 4 --steps 16
.venv-cudaq/bin/python python/run_parallel_refinement.py --sizes 2 4 --steps 64
.venv-cudaq/bin/python python/cudaq_star_gadget.py
.venv-cudaq/bin/python python/cudaq_schedule_reproducer.py
.venv-cudaq/bin/python python/cudaq_minimal_timing.py

# qb_gsee / science Python:
python python/small_dynamics.py
python python/noisy_small_observable.py
python python/cultivation_model.py
python python/star_parallel.py
python python/plot_refinement_v3.py
python -m pytest -q --junitxml=results/refinement_v3/tests.xml
python python/audit_refinement_v3.py --write-manifest
```

The normal v3 audit without `--write-manifest` verifies the frozen evidence;
the write flag is for deliberately replacing that evidence after a scientific
rerun. Current tests may be more numerous than the archived test report.

## Earlier calculations and compiler environment

```bash
python python/make_plots.py

# Create the compiler environment only if it does not exist:
python -m venv .venv-cudaq
.venv-cudaq/bin/python -m pip install -r configs/requirements-cudaq.lock

.venv-cudaq/bin/python python/cudaq_compile.py --sizes 2 4 --full-synthesis
.venv-cudaq/bin/python python/cudaq_compile.py --sizes 6 8
python python/plot_compiled.py
.venv-cudaq/bin/python python/cudaq_compile.py --manifest-only
```

The compiler environment runs on CPU; no CUDA-Q source build or GPU is needed
for these paths. Encoded backend diagnostics retain rejected statuses and are
not used as physical timing evidence. See the historical
[CUDA-Q review](archive/CUDA_Q_RESULTS.md) and
[v3 review](archive/REFINEMENT_V3.md) for their scope.

Full v1 MLIR dumps are stored as `.mlir.gz`; the small synthesis fixtures stay
plain text. `python python/package_results.py --restore` verifies and restores
the originals using `results/cudaq_v1/storage.json`. It preserves different
local data by refusing to overwrite it. After an intentional compiler rerun
and manifest refresh, `python python/package_results.py` refreshes the compressed
copies. Tests and the current figures do not require restoring the full dumps.

## Repository reorganization

Python modules and command-line scripts are in `python/`; tests are in `tests/`.
Imports remain simple script imports, with `pytest.ini` adding `python/` for
tests. All data/config paths resolve from the repository location, independent
of the caller's working directory, except explicit user-supplied relative paths.
`configs/file_locations.json` maps old root filenames for historical provenance.
`results/presentation_v1/reorganization.json` records the pre-move source hashes.

`PLAN.md` remains at the root because it had an active editor lock during the
move. Other planning drafts are in `docs/plans/`. The accepted decisions are in
[WORKING_DECISIONS.md](WORKING_DECISIONS.md); chronological reviews are in
`docs/archive/`. Earlier TeX source and PDFs are in `tex/archive/`.

## Delivery menus and the large finite patch

`make delivery` reoptimizes the encoded layouts and writes
`results/delivery_v1/` plus the delivery figures. It uses the saved v3 physics
coefficients so the layout effect is isolated. The source-count menu, block
rates, and analog-lane choices are in `configs/delivery_sensitivity.json`.
`python python/delivery_sensitivity.py --check` audits the saved error budgets,
logical/magic bounds, source-menu endpoint, replicas, and throughput conventions.

The large solver needs GCC with C++17/OpenMP, plus the science Python environment.
Build and run its independent controls before submitting a new calculation:

```bash
mkdir -p build
g++ -O3 -std=c++17 -fopenmp cpp/hubbard_sector.cpp -o build/hubbard_sector
python -m pytest tests/test_large_dynamics.py -q
sha256sum cpp/hubbard_sector.cpp python/large_dynamics.py build/hubbard_sector python/run_checked_dynamics.py > build/dynamics_inputs.sha256
```

The tests compile a separate temporary executable and do not replace a binary
used by an active PBS job. The production driver rejects full 4 × 4 runs outside
a PBS allocation. Query current limits, reservations, job counts, and allocation
balance before submitting; scripts contain the workspace path and pinned Python.

```bash
qsub -r y -q backfill -o results/dynamics_4x4/smoke.pbs.log pbs/dynamics_smoke.pbs
# Inspect terminal exit status and all four full-size smoke records first.
qsub -r y -q backfill -o results/dynamics_4x4/checked_steps.pbs.log pbs/dynamics_checked_steps.pbs
```

The checked-step driver reuses only contract-matching completed cases and runs
r = 8,16,32 for both orderings. It adds r = 64 only if the two finest coefficients
disagree by more than 5%. `numactl --interleave=all` avoids placing the entire
state vector in one memory domain. The saved same-node benchmark compares
identical one-step circuits and records the speedup and numerical agreement.
The initial run was stopped after a converged saved plaquette r = 32 point;
its intentional nonzero PBS exit is recorded separately from successful resumes.
Do not run two drivers into the same case directory concurrently.

For a fresh run with changed inputs, use a new output directory (and adapt the
checked-step runner path). Existing results deliberately reject a different
source/binary contract. Each completed case is atomic; an interrupted case
restarts from its initial state. No state-vector checkpoint is retained.

After the step menu is complete:

```bash
python python/large_patch_analysis.py
make paper
make audit
```

This produces direct measured-step costs and a separately labeled interpolation
sensitivity at L = 4, tau = 1. The main capacity curves retain their original
finite-patch transfer assumptions. The separate follow-up below extends time.

## Direct headline checks and the time extension

`make patch-study` regenerates `results/dynamics_followup_v1/` analysis and the
patch/time and capacity figures from saved raw cases. It does not launch jobs.
The [method note](DIRECT_CHECKS_AND_TIME_EXTENSION.md) distinguishes validation
of the original allocations from subsequent measured-error reoptimization.

`configs/dynamics_followup.json` fixes the selected row set, observed site,
thread count, and convergence rule. `configs/patch_sensitivity.json` fixes the
capacity sweep. The three completed submissions are recorded in
`results/dynamics_followup_v1/submissions.json`, with terminal accounting in
`jobs.json`. Each task owns a separate case directory.

For a new scientific run, use a new result directory and update the runner/PBS
paths and input manifest together. Follow the cluster preflight and smoke
requirements above; do not duplicate active writers or overwrite the frozen
production executable. The existing submission pattern is:

```bash
qsub -r y -q backfill -v TASK=selected -o results/dynamics_followup_v1/selected.pbs.log pbs/dynamics_followup.pbs
qsub -r y -q backfill -v TASK=tau2 -o results/dynamics_followup_v1/tau2.pbs.log pbs/dynamics_followup.pbs
qsub -r y -q backfill -v TASK=tau4 -o results/dynamics_followup_v1/tau4.pbs.log pbs/dynamics_followup.pbs
```

The scripts check `inputs.sha256` before execution. That file includes the
source, executable, runner, configuration, and original selected-row input.
Rebuilding on another system requires a new binary contract and result
directory; the saved production binary itself is not distributed.

After terminal completion and the case census:

```bash
make patch-study
python python/dynamics_followup_analysis.py --check
make test
make paper
make audit
```
