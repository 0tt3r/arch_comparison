# Reproducing the calculations

Run commands from the repository root. Python 3.12 with NumPy, SciPy,
Matplotlib, pytest, and packaging is sufficient for the models, tests, plots,
and manuscript helper. The existing `qb_gsee` environment provides these.
Set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` for the bounded examples.

## Current figures and manuscript

For the fast analytical discussion example alone, run `make discussion`.
Edit `configs/discussion_model.json` for the scenario knobs. This regenerates
`results/discussion_v2/` and the two `discussion_capacity` figures without
calling CUDA-Q or a dynamics solver. See [DISCUSSION_MODEL.md](DISCUSSION_MODEL.md).

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
