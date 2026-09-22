# Parallel compilation, accuracy, cultivation, and STAR retries

The important change is that accuracy now constrains the architecture comparison.
At tau = 1, the shared four-step plaquette formula has local doublon errors
0.01364 on 2 × 2 and 0.01439 on 3 × 3. Both exceed even the loose 0.01 total
error target. Thus the r = 4 architecture curves remain useful conditional cost
comparisons, but cannot establish a crossover for an accurate dynamics task.
The smaller patches do not determine the error at L = 4 or in the bulk.

The [derivation note](../../tex/archive/refinement_v3.pdf) uses the dedicated
[bibliography](../../tex/archive/refinement_v3_references.bib), together with the earlier bib files.
All new evidence is in [results/refinement_v3](../../results/refinement_v3).
The earlier result directories are preserved.

## Parallel native compilation

`parallel_nisq.py` preserves the ordered plaquette evolution, gathers independent
fermionic blocks using parallel odd-even transposition sorting, and fuses Gaussian
operations on adjacent pairs. Every swap is a physical fermionic swap on a snake
embedded in the nearest-neighbor square. The compiler chooses among six initial
mode-order/local-routing candidates. Arbitrary Rzz and one-qubit pulses retain
the declared 10 ns duration, p2 = 0.001, p1 = p2/10, and idle p per 100 ns.

The observable-specific circuit additionally removes the first onsite phase
(scalar on the initial Fock state), the final onsite phase (commutes with doublon),
and blocks outside the exact backward fermionic cone. Changing initial mode order
is free here because the initial state is a product of occupations; final mode
permutations are retained explicitly. These optimizations are available to every
architecture. The present encoded curves retain their original full workload,
so they have not yet received these reductions.

| L = 4, tau = 1, r = 4 | Native entanglers | 10 ns layers | Local PEC exposure | Week-cutoff Q |
|---|---:|---:|---:|---:|
| v2 greedy route + cone | 3,149 | 2,187 | 35.297 | 2.07e13 |
| v3 parallel route, full evolution | 3,204 | 1,074 | 22.214 | 2.13e7 |
| v3 parallel route + observable reductions | 2,648 | 964 | 20.017 | 2.13e6 |

The last exposure is 10.5927 (two-qubit), 4.1201 (one-qubit), 5.1681 (idle),
and 0.1361 (preparation/readout). The zero-idle value for this same circuit is
Q = 12,168. These changes combine scheduling, mode placement, fusion, and
observable algebra; they are not solely a scheduling speedup. The route still
uses only the snake edges, so it is not an optimized two-dimensional compiler.

Nine new native circuits were compiled with CUDA-Q Logical: L = 2, 4, 6, 8 at
r = 4 with observable reductions; L = 4 without them; and L = 2, 4 at r = 16, 64.
The [compiler figure](../../results/refinement_v3/parallel_compilation.png) and
[CSV](../../results/refinement_v3/parallel_summary.csv) retain the idle sensitivities.
All costs here describe a specified local-PEC estimator, not a lower bound on
NISQ computation or a no-go result for other mitigation strategies.

## Accuracy and finite spatial patches

`small_dynamics.py` evolves the sparse Hubbard Hamiltonian in fixed spin-number
sectors, independently of the compiled gate decomposition. It also evolves the
same second-order pink/gold/onsite product formula. The largest sector has 15,876
basis states (nine sites). The saved results contain 112 Trotter comparisons,
164 points on exact patch dynamics curves, and four nested spatial comparisons.

The refinement scenario splits the existing 25% systematic reserve equally
between Trotter and spatial truncation. At epsilon = 0.001, each gets 0.000125;
this replaces the old reserve rather than enlarging the total error budget.
Signal-relative budgets using each exact patch value are also recorded.

| Patch at tau = 1 | Exact D | Error, r = 4 | Error, r = 16 | Error, r = 64 |
|---|---:|---:|---:|---:|
| 2 × 2 | 0.1179790 | 0.0136397 | 0.00087075 | 0.00005447 |
| 3 × 3 | 0.1050849 | 0.0143939 | 0.00085350 | 0.00005317 |

For these patches, r = 64 is the first tested power of two passing the strict
split Trotter budget; r = 16 passes its epsilon = 0.01 counterpart. This does
not establish the minimum r, and cannot be extrapolated to L = 4 without further
analysis. The asymptotic local error is consistent with r^-2; coarse steps can
be nonmonotonic. Compiled L = 2 PEC cost increases from one nine-qubit copy at
r = 16 to about 4.73e6 qubits at r = 64. This is the accuracy/noise tradeoff for
the specified estimator, with Trotter error independently checked at L = 2.

For physical truncation, the five-site cross is nested inside 3 × 3, with the
same Neel pattern and central measured site. Their exact D difference is
0.000475 at tau = 0.1 and 0.07935 at tau = 1, above the strict truncation
allowance. These are finite-reference differences; the 3 × 3 outer-boundary
error is not bounded. In particular, this does not certify a physical light
cone or establish that 3 × 3 is sufficient. The [accuracy figure](../../results/refinement_v3/accuracy_calibration.png)
keeps these limitations visible.

## Actual post-selection bias on eight qubits

`noisy_small_observable.py` propagates the eight-qubit density matrix through
the saved native circuit with primitive gate, idle, preparation, and readout
channels. It measures both acceptance and doublon bias, rather than inferring
bias from acceptance. At L = 2, r = 4:

| Selection | Acceptance | Bias relative to noiseless Trotter circuit |
|---|---:|---:|
| None | 1 | 0.020371 |
| Both spin parities | 0.81958 | 0.005808 |
| Exact numbers N_up = N_down = 2 | 0.80943 | 0.004212 |

Number selection helps but fails the strict precision target. Including the
separately measured Trotter error, its bias relative to exact dynamics is
0.01785. At r = 64 its noise bias is 0.1070. The [post-selection figure](../../results/refinement_v3/postselection_bias.png)
shows all three step counts. These calculations concern one eight-qubit noise
model; they do not rule out more advanced mitigation, different compilation,
or selection at L = 4. Parity acceptance independently agrees with the analytic
Pauli-character calculation, and density trace/positivity checks pass.

## Cultivation is available to both FT and hybrid

The new source model uses Gidney, Shutty, and Jones,
[arXiv:2409.17595v1](https://arxiv.org/abs/2409.17595v1), Table 1 and Figs. 1, 15,
16: output error about 2e-9 at p = 0.001, and a Fig. 15 volume of 6.6e4
qubit-rounds with early rejection and sharing. The source ends in a grafted
matchable code; transfer, storage, and consumption are separate obligations.

`cultivation_assumptions.json` distinguishes sourced values from our supply
assumptions: rounded packed volume 1e5 qubit-cycles, 1,024 physical qubits,
50-cycle successful latency, and a 2d-cycle handoff. An unshared fixed-slot
scenario costs 5.12e6 qubit-cycles; additional handoff error is swept over 0,
1e-9, and 1e-8. Mapping the source paper's cycles and uniform noise to our
240 ns cycle and p1 = p/10 is conditional. No physical cultivation factory
has been compiled for this device.

Both arms may consume cultivated states directly or add one 15-to-1 stage
when their actual state count violates the magic-error budget. With zero added
handoff error, the unamplified strict budget supports 37,500 states at 2e-9
per state. The compiled L = 4 FT circuit uses 19,584; L = 6 needs the extra
stage. A handoff error of 1e-8 already forces it at L = 4. The cutoff is therefore
task- and budget-dependent. TopQAD Fig. 16's roughly 1e5 cutoff counts
small-angle rotations (about 100 T states each), not 1e5 T states.

| L = 4, r = 4 reference | Q at epsilon = 0.001 | Q at epsilon = 0.01 |
|---|---:|---:|
| Serial STAR | 166,698 | 45,486 |
| FT, conventional serial factory | 13,791,250 | 216,890 |
| FT, packed cultivation, zero handoff error | 748,656 | 52,180 |
| Hybrid, conventional serial factory | 669,714 | 54,150 |
| Hybrid, packed cultivation, zero handoff error | 57,472 | 47,232 |

These are conditional r = 4 resource numbers, not successful dynamics runs.
The loose-precision FT calculation retains the strict-precision synthesized
words, so it can still be improved by resynthesis. The
[encoded figure](../../results/refinement_v3/cultivation_parallel.png) shows supply
and handoff sensitivities; the shading is a scenario range, not uncertainty
estimated from data.

## Explicit logical STAR retries and parallel scheduling

`cudaq_star_gadget.py` compiles ancilla preparation, CNOT, measurement, conditional
angle doubling, and a returned exhaustion flag at caps 1, 2, 4, 8, and 24.
This is actual logical P0 control flow. It is not a surface-code injection,
patch-growth, or nearest-neighbor lattice-surgery implementation.

For the stated logical instrument, every single two-qubit Pauli fault before
or after CNOT leaves each outcome at probability 1/2 and each conditional map
Z-covariant. The calculation supports cancellation of the Z-dephasing PEC
inverse across this local sandwich under branch-independent Pauli noise.
It does not certify those assumptions for a physical STAR implementation.
The old worst-case Gamma*q_L bound remains the default, with a separately
labeled covariance-based scenario.

`star_parallel.py` schedules disjoint **full JW supports**, retains dependencies,
charges six injection tiles per concurrent lane, and uses the expected maximum
of the geometric retry counts in each batch. Nine assumed clocks per attempt
recover the previous 18-clock serial mean. At L = 4 and strict precision, two
lanes minimize Q at 121,716, versus 166,698 for the old serial reference. Four
and eight lanes reduce time further but lose on area. At loose precision one
lane still wins. This is a useful counterexample to simply dividing the
18-clock time by the number of rotations that look independent at fermion level.
Asynchronous dispatch and alternative logical mappings remain open refinements.

## CUDA-Q backend timing and reproducibility

`cudaq_schedule_reproducer.py` reduces the encoded scheduling failure to four
fixed-angle rotations repeated four times (16 T states). The installed backend
rejects it with `nested schedule entry must lie within its parent interval`.
The 4-, 8-, and 32-T cases pass a loose timing sanity bound. Separate two-resource
template-call controls all return the expected durations. The root cause is
unresolved; these results do not validate the passing macro schedules.
Exact input operations, P3 IR, and the rejected diagnostic are retained.

The 73 tests pass in `qb_gsee`, including independent dense/sector dynamics,
fermionic permutation signs, observable pruning, explicit Pauli-mixture noise,
STAR instrument algebra, and magic-state budget checks. `audit_refinement_v3.py`
replays the saved native counts/exposures and verifies source/artifact hashes.

Reproduce from the repository root, setting `OPENBLAS_NUM_THREADS=1` and
`OMP_NUM_THREADS=1` for these small calculations:

```bash
.venv-cudaq/bin/python run_parallel_refinement.py --full-L4
.venv-cudaq/bin/python run_parallel_refinement.py --sizes 2 4 --steps 16
.venv-cudaq/bin/python run_parallel_refinement.py --sizes 2 4 --steps 64
.venv-cudaq/bin/python cudaq_star_gadget.py
.venv-cudaq/bin/python cudaq_schedule_reproducer.py
.venv-cudaq/bin/python cudaq_minimal_timing.py

# Python from the qb_gsee environment:
python small_dynamics.py
python noisy_small_observable.py
python cultivation_model.py
python star_parallel.py
python plot_refinement_v3.py
python -m pytest -q --junitxml=results/refinement_v3/tests.xml
tectonic refinement_v3.tex
python audit_refinement_v3.py --write-manifest
python audit_refinement_v3.py
```

The next scientifically decisive steps are a larger-patch or controlled local
truncation estimate, a common accuracy-qualified workload for all arms, and
physical STAR growth/handoff schedules with a calibrated residual channel.
The new results identify those dependencies; they do not settle the ultimate
NISQ limit, the STAR crossover, or the optimal encoded layout.
