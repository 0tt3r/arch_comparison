# Hubbard dynamics on a fixed qubit array

For the team discussion, start with the **[physical-error and algorithm comparison](results/physics_sweep_v1/summary.md)**
and the **[main manuscript](tex/main.pdf)**. The analytical
[assumptions and formulas](docs/DISCUSSION_MODEL.md) remain the common model.

![Achievable lattice size versus physical qubits, physical error, and precision](figures/physical_error_capacity.png)

The new sweep varies p = 0.001, 0.0003, and 0.0001 at 0.01 and 0.001 absolute
error. All 144 matching p = 0.001 cases reproduce the frozen baseline.
Source quality and STAR floors at lower p are explicit scenarios, not new
gadget calibrations. Both precision targets remain prominent.

The [algorithm comparison](figures/hamming_weight_comparison.pdf) finds that
Hamming-weight phasing lowers the strict 4 × 4 FT estimate from 1.88 million
to 278,192 qubits, including workspace and arithmetic. At loose precision,
individual synthesis remains cheaper for 4 × 4 and 8 × 8. Compact-encoding
counts and an independent small-patch ordering diagnostic are retained in the
[same results directory](results/physics_sweep_v1/summary.md).

The **[workshop reconciliation](results/workshop_reconciliation_v1/summary.md)**
evaluates all three codebases at 4 × 4, tau = 1,2, common steps, 0.01 absolute
error, and a 240 ns syndrome cycle. Compilation, estimator, and delivery
differences remain explicit; similar headline capacities do not establish an
independent validation. A small source archive makes the comparison reproducible
without modifying the supplied repositories.

The **[direct-check summary](results/dynamics_followup_v1/summary.md)** and
[patch-capacity figure](figures/patch_capacity.pdf) retain the finite-size study.

The current comparison includes exact 4 × 4 calibration at tau = 1,2,4.
It shows the effect of adding larger patches while holding reference delivery and source
menus fixed. Both curves optimize the error allocation and step count. All
16 saved tau = 1 headline allocations now pass their directly measured
Trotter allowances. See the [method and scope](docs/DIRECT_CHECKS_AND_TIME_EXTENSION.md).
Larger-lattice transfer and effective NISQ residual bias remain assumptions.

The earlier **[discussion-v3 reference](results/discussion_v3/summary.md)**
retains the 3 × 3 planning envelope and the longer time grid:

Its curves optimize each arm's error allocation, step count, and resource-state
supply. Small-patch Trotter coefficients are measured through tau = 8 for both termwise
and plaquette orderings; the planning envelope avoids isolated small-patch
cancellations. See the [calibration figure](figures/time_calibration.pdf).
Larger-lattice transfer and effective NISQ residual bias remain assumptions.
The [0.01 absolute-precision sensitivity](figures/discussion_capacity_loose.pdf)
is also available. The fourth panel follows tau = L/4, stopping at L = 32 to
stay within the measured time interval. Run `make discussion` to regenerate
these estimates from saved calibration.

The new [encoded-delivery comparison](figures/delivery_capacity.pdf) tests
whole-block consumption rates, analog lanes, and shared source supply.
[Details and the large-patch method](docs/DELIVERY_AND_LARGE_PATCH.md) separate
published block rates from uncompiled parallelism assumptions. The
[completed 4 × 4 calculation](results/dynamics_4x4/summary.md) now checks
finite-lattice transfer at tau = 1; [this figure](figures/large_patch_calibration.pdf)
shows the ordering-dependent difference from the smaller patches. The
[time extension](figures/patch_time_extension.pdf) adds the tau = 2,4 checks.

Start with the **[main manuscript](tex/main.pdf)** ([LaTeX](tex/main.tex),
[bibliography](tex/references.bib)). It brings together the analytical resource
laws, CUDA-Q compilation, mitigation, STAR, FTQC and hybrid supply, and the
small-system accuracy checks.

The principal comparison is **maximum modeled lattice size versus available
physical qubits**, subject to a one-week runtime cutoff.

The separate [compiled-reference capacity plot](figures/capacity_vs_qubits.pdf)
retains the more detailed previous implementations.

These are conditional four-step circuit costs. Four steps fail the accuracy
target on the small square patches checked; the current 4 × 4 calibration uses different step counts.
The staircases show only tested sizes, and reaching 8 × 8 marks the end of the
current sweep. The [underlying CSV](results/presentation_v1/capacity_thresholds.csv)
retains the original area thresholds.

The NISQ curve is a **sufficient gatewise-PEC sampling budget**, not an
optimized NISQ limit. Its 4 × 4 circuit uses 32 active qubits (36 including
square packing); the strict-precision estimate buys 59,239 replicas to
complete the prescribed shot budget within a week. See the
[experimental comparison and cost breakdown](docs/NISQ_EXPERIMENT_COMPARISON.md)
for why this does not predict the resources used in the Phasecraft experiments.

The [architecture overview](figures/architecture_overview.pdf) and
[patch/supply detail](figures/architecture_detail.pdf) explain resource roles.
The [rotation protocols](figures/rotation_protocols.pdf) show the distinction
in time: STAR injects angle-specific states with retries, while FTQC consumes
high-fidelity T states in a synthesized sequence. Both encoded approaches
share protected Clifford operations; STAR preparation can reuse surgery
ancillas. These drawings are schematic and available as editable SVG,
vector PDF, and PNG in `figures/`.

| Location | Purpose |
|---|---|
| `python/` | Models, compilers, plotting, build and audit scripts |
| `cpp/` | Matrix-free fixed-spin-sector dynamics kernel |
| `pbs/` | Bounded Improv smoke and calibration jobs |
| `tests/` | Independent mathematical and implementation checks |
| `configs/` | Assumptions, compiler dependency lock, historical filename map |
| `tex/main.tex` | One complete current manuscript |
| `tex/references.bib` | Dedicated bibliography for the main manuscript |
| `tex/archive/` | Earlier complete notes and their original PDFs/bibliographies |
| `figures/` | Current capacity plots, algorithm comparisons, and architecture cartoons |
| `results/` | Versioned numerical evidence, compiler artifacts, and manifests |
| `docs/` | Accepted decisions, reproduction guide, older plans and reviews |
| `refs/` | Reference papers |

Use the existing science environment, then:

```bash
conda activate qb_gsee
make test
make paper
make audit
```

`make paper` regenerates the current figures and builds `tex/main.pdf`.
It does not rerun the scientific calculations. To specify an interpreter,
use `make PYTHON=/path/to/python paper`. The build helper uses Tectonic on
PATH or the existing local package cache; `latexmk -pdf main.tex` from `tex/`
also works with a standard TeX installation.

See [REPRODUCING.md](docs/REPRODUCING.md) for scientific reruns, the isolated
CPU CUDA-Q environment, archived-note builds, and compressed MLIR handling.
[WORKING_DECISIONS.md](docs/WORKING_DECISIONS.md) records the accepted conventions.
`PLAN.md` remains at the root because it had an active editor lock during the
reorganization. The latest chronological calculation report is
[REFINEMENT_V3.md](docs/archive/REFINEMENT_V3.md); its source datasets and original
hash manifest remain intact.
