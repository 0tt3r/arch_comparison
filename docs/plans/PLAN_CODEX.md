**Plan: resource scaling for Hubbard dynamics on a common qubit array**

Written 2026-09-22 from the workshop discussion. This records the proposed study and the agreed baseline; no crossover estimates or numerical results have been produced. It is an independent plan for comparison with the other agent's proposal.

**Scientific objective**

Determine how much two-dimensional Fermi–Hubbard dynamics can be computed, at a specified observable accuracy, using three different approaches on the same physical hardware:

1. Physical-qubit computation with error suppression, post-selection, and error mitigation.
2. Partially fault-tolerant STAR computation, with protected Clifford operations and analog rotations.
3. Fully fault-tolerant surface-code computation, including the resources needed for non-Clifford operations.

The aim is to understand scaling laws, controlling quantities, and crossovers, with sufficiently calibrated prefactors to make the comparison meaningful. Detailed engineering should be added where it changes those conclusions. The existence and ordering of three winning regimes are hypotheses to test; STAR may or may not have a winning window.

**Common hardware and comparison variables**

- Physical hardware: a square, two-dimensional nearest-neighbor array.
- Gate duration: `t_g = 10 ns`.
- Measurement duration: `t_meas = 100 ns`.
- Nominal physical error probability: `p = 1e-3`.
- Runtime cutoff: one week, or `604800 s`, to obtain the requested answer, including the required sampling.

Keep application size and hardware size distinct:

- `N = L^2`: number of Hubbard sites.
- `Q`: number of physical qubits available in the hardware array.
- Spinful Hubbard has `2N` occupation modes before ancillary and encoding overhead.

At fixed `Q`, allow every approach to allocate spare area usefully. In particular, physical-qubit computation can run independent circuit replicas in parallel, while encoded approaches can trade area between data, routing, factories, and replicas.

The hardware specification needs a minimal common noise convention. Before numerical comparisons, declare the native gate set, how `p` applies to gates, preparation and measurement, and an idle-noise model tied to elapsed time. Declare reset duration and feedforward latency as well. Treat unresolved choices as sensitivity parameters. The meaning of idle noise during a 100 ns measurement can materially affect the comparison.

A surface-code round must include an explicit syndrome-extraction schedule. Its duration cannot simply be set to the 10 ns physical gate time. Use the same physical timing assumptions in all three approaches.

**Agreed baseline and proposed concrete defaults**

The user accepted a simple initial state, a local density or doublon observable, and `0.01` absolute error, with lattice size and evolution time swept. The following choices make that baseline concrete:

- Hamiltonian: nearest-neighbor spinful Hubbard,
  `H = -J sum_<ij>,sigma (c†_i,sigma c_j,sigma + h.c.) + U sum_i n_i,up n_i,down`.
- Geometry: even `L x L` lattices, initially with open boundaries and zero flux.
- Initial state: a half-filled Néel occupation product state.
- Interaction: `U/J = 4` initially.
- Observable: local doublon density `D_i(t) = <n_i,up n_i,down>`, initially at a specified central site. Prefer this over a total local charge density that could be fixed by symmetry for the chosen quench.
- Accuracy: absolute error at most `0.01` at 95% confidence for a specified observable and time point.
- Evolution variable: `tau = J t / hbar`.

Separate statistical uncertainty from systematic error: product-formula error, synthesis error where applicable, residual hardware/mitigation bias, and state-preparation/readout error. Allocate and optimize the total error budget consistently. A narrow statistical interval alone does not meet the total-error target.

If the requested output becomes an entire curve or many observables, explicitly change the success criterion and account for simultaneous coverage and measurement reuse. Do not silently compare one architecture's single-point estimate with another's full time series.

Sweep both size and evolution time. Include fixed-time slices and a family of times that grows with system size, such as a propagation time across the lattice. Fixed-time local observables may only require a bounded region to the desired accuracy, so exponential cost in the full lattice size is not automatic. Any locality-based reduction should be available to all three approaches.

**Resource objective and accounting**

The headline quantity is

`Q_min(N, tau, epsilon; T_answer <= 604800 s)`.

Also report the maximum reachable evolution time at fixed hardware size, and retain the qubit/runtime tradeoff rather than reducing everything to one space-time product.

A useful schematic runtime model is

`T_answer ≈ N_0(epsilon_stat, delta) * Gamma^2 * t_attempt / (a * m)`.

Here `N_0` is the ideal sampling requirement, `Gamma^2` is the mitigation variance multiplier, `a` is whole-circuit acceptance probability, and `m` is the number of concurrent replicas. For bounded observables, ordinary sampling gives an inverse-square dependence on statistical precision, up to variance and confidence factors.

Use this factorization only where justified. Estimate acceptance, conditional bias, and estimator variance together when post-selection and mitigation interact. Include state-factory rejection and repeat-until-success delays in factory throughput and execution schedules; do not also charge them as whole-circuit rejection. Include startup costs where relevant and respect shared-factory bottlenecks.

At each problem point, optimize Trotter step count, error-budget allocation, code distances, factory capacity, and parallel replicas subject to the physical-qubit budget. For expectation estimation, logical faults must satisfy the observable-error criterion; demanding that no fault occur anywhere in the entire sampling campaign is a stronger, different requirement.

**Mechanisms and hypotheses to test**

| Approach | Candidate controlling quantity | Main issue to resolve |
| --- | --- | --- |
| Physical qubits | Effective noisy circuit volume relevant to the observable | Mitigation variance and residual bias after exploiting structure |
| STAR | Total absolute rotation angle plus any per-rotation error floor | When angle-dependent errors dominate, and when preparation/growth costs restore a gate-count penalty |
| Full FT | Protected circuit volume, synthesis demand, and factory throughput | Size of the initial overhead and subsequent scalable error suppression |

For physical-qubit computation, compare an explicit, limited menu: unmitigated execution where accurate enough, symmetry verification/post-selection, and a probabilistic error cancellation (PEC) reference model. Add extrapolation or training-based mitigation when it provides an identifiable benefit with assessable bias. The result is a comparison of specified strategies and favorable envelopes, not an impossibility claim about every future NISQ method.

Under suitable local stochastic noise and PEC assumptions, a useful starting model is `log(Gamma^2) ~ c p V_eff`. The effective volume depends on the observable, circuit, and mitigation strategy. General mitigation lower bounds support exponential penalties under stated assumptions; they do not directly establish this Hubbard task's crossover. Post-selection only detects some faults: acceptance probability cannot be equated with probability of an error-free result.

For STAR, distinguish the original injection protocol from improved small-angle protocols. In the improved effective dephasing model, the residual rotation error can behave as `q_rot(theta) ≈ alpha p |theta|`, giving

`Gamma_STAR^2 ≈ exp(4 alpha p Theta)`, with `Theta = sum_j |theta_j|`.

Use one explicit angle convention throughout. For rotations written as `exp(-i theta P)`, a straightforward positive-time product formula has `Theta ≈ t sum_l |h_l|` for `H = sum_l h_l P_l` with `hbar = 1`. Thus, at fixed order, increasing the Trotter step count need not increase the leading angle-dependent mitigation exponent. For a full-system Hubbard decomposition at fixed `U/J`, this angle budget is extensive in `N tau`. Clifford extraction, cancellations, and alternative decompositions can change the actual angle cost.

A more realistic effective model may include

`q_rot(theta) ≈ alpha p |theta| + q_floor`,

so that

`log(Gamma_STAR^2) ≈ 4 alpha p Theta + 4 N_rot q_floor`.

The TopQAD study explicitly includes an angle-independent patch-growth contribution. The floor depends on the protocol and code choices; it is not a universal fixed constant. Optimize preparation distance, growth, acceptance, RUS corrections, buffering, and delivery rate. Account for protected Clifford errors too. A central scientific question is whether total angle or the accumulated floor controls the practical turnover.

For full FT, begin with a declared below-threshold surface-code model such as

`p_L(d) ≈ A (p / p_th)^((d+1)/2)`.

Use a consistent unit for `p_L`, such as patch-rounds, when accumulating protected circuit errors. This leads to logarithmic distance growth with protected volume and inverse error budget under the model, with data-patch area scaling approximately as `d^2` and many logical operations taking order `d` syndrome rounds. Treat the prefactor and threshold as model inputs requiring calibration or sensitivity analysis.

Independent rotation synthesis gives a useful reference scaling of roughly `N_rot log(N_rot / epsilon_synth)` T gates. Include distillation or another explicitly selected FT state-production protocol, routing, storage, and supply-limited execution. Full FT still pays ordinary measurement sampling costs. Its scalability does not imply feasibility at every fixed `Q`.

**Algorithm and compilation strategy**

First compare a matched second-order Trotter workload across the three approaches. Match the Hamiltonian, state, observable, and accuracy target; allow architecture-specific step counts and compilation choices within those requirements.

Retain width, scheduled depth, Clifford work, rotation-angle histograms, total absolute rotation angle, routing work, and ancillary requirements. Square-lattice fermions do not automatically produce a geometrically local qubit circuit. Compare a straightforward Jordan–Wigner implementation with a locality-conscious fermionic compilation so that a mapping artifact is not mistaken for an architectural scaling law.

Then form an algorithm-optimized comparison. In particular, add qubitized dynamics for FT and evaluate whether it helps the other approaches. Include PREPARE/SELECT, ancillas, synthesis, success overhead where applicable, and physical routing. Keep the matched-algorithm and optimized-algorithm results separately identifiable: STAR's small-angle advantage and FT's ability to use a different algorithm are both scientifically relevant.

**Software roles and inspected starting points**

- CUDA-Q Logical is a candidate circuit/compilation and resource-analysis backbone. The local tree is `~/soft/cuda-quantum`, with the toolkit under `preview/logical/`.
- `preview/logical/examples/03_fermi_hubbard.py` is a simple chain example, not the desired square-lattice benchmark.
- `preview/logical/python/cudaq/logical/algorithms/fermi_hubbard.py` supplies a square-lattice Jordan–Wigner workload, but its current workload is controlled evolution with first-order Trotterization. Adapt it to ordinary observable dynamics and the selected product formula.
- `surface_physical_target` in `preview/logical/python/cudaq/logical/targets/__init__.py` is explicitly a compact reference realization with native Clifford actions and multi-Pauli measurements. Its output is not automatically a cost for our nearest-neighbor gate array. Supply an appropriate device/gadget realization or connect its logical output to a separately calibrated architecture model.
- pyLIQTR at `~/soft/pyLIQTR` includes `Examples/ApplicationInstances/FermiHubbard/fermi_hubbard-dynamics-qubitized.ipynb`. Use it for the qubitization comparison when needed; there is no initial requirement to move all workloads into one package.
- Keep a small, inspectable cost model that can consume workload summaries independently of either package. Detailed physical compilation should resolve important prefactors and scheduling effects without becoming a prerequisite for the first scaling plots.

These are source-inspection findings, not a claim that the software has been installed, run, or validated for this study.

**Experimental reference: arXiv:2510.26300v2**

The Phasecraft/Quantinuum experiment provides a useful executed dynamics workload. It uses a `7 x 4` doubly periodic nearest-neighbor Hubbard lattice, a threaded pi flux, an initial entangled dimer state with holon/doublon defects, and `U/J = 0, 4`, with evolution through `Jt = 2`. Its circuits use four second-order Trotter steps and at most 2,415 two-qubit gates.

Three distinctions matter for using it here:

1. Increasing the long dimension at fixed width four produces a fixed-width scaling family. Our principal two-dimensional study should grow both dimensions at fixed aspect ratio. Boundary conditions, flux, and odd/even periodic dimensions must be declared when extrapolating the experiment.
2. The experiment exploits effective all-to-all ion-trap connectivity and native variable-angle entangling gates. Some SWAP operations are implemented by relabeling and omitted from its gate counts. Recompile the workload for the common physical array rather than transferring the published count or runtime directly.
3. Training-based mitigation and Gaussian-process smoothing contribute to the reported results, while interacting-case Trotter error is assessed indirectly. Experimental execution does not establish our `0.01` total-error criterion. Reassess approximation error, residual mitigation bias, and sampling uncertainty under the common criterion, including training costs and assumptions.

Use three related workloads:

| Purpose | Workload |
| --- | --- |
| Main scaling study | Agreed square-lattice baseline with local doublon density |
| Experimental reference point | Exact `7 x 4` Hamiltonian, initial state, and selected observables, recompiled for our hardware |
| Initial-state sensitivity | A few square-lattice cases with a similar dimer/holon–doublon preparation |

The dimer/defect extension creates a localized disturbance whose charge and spin dynamics can be followed. Neighboring dimers can be prepared with depth bounded independently of lattice size. Keep this as a controlled extension of the agreed baseline.

**Execution sequence and deliverables**

1. **Analytical prototype.** Implement the transparent runtime/area model, identify dimensionless noise and rotation exposures, and produce provisional crossover diagrams with assumptions visible. These plots are model predictions.
2. **Matched workload.** Generate the square-lattice circuits and verify conventions, initial state, observables, and product-formula behavior against exact dynamics at small sizes. Cross-check compiled resource summaries against direct counts. Distinguish rigorous error bounds from empirical small-system calibration.
3. **Calibrated architecture costs.** Check representative nearest-neighbor schedules, syndrome rounds, logical operations, STAR preparation/RUS, and FT state supply. Use small noisy examples to assess mitigation bias, variance, and acceptance. Calibrate the constants that can move the crossover, then sweep the model without simulating every large instance.
4. **Experimental reference point.** Reconstruct the selected `7 x 4` task and compare costs under the common hardware and error criterion. Add a few dimer-state square-lattice cases.
5. **Algorithm and sensitivity comparison.** Add qubitized dynamics and test which conclusions survive changes in noise assumptions, observable, time scaling, precision, mapping, and factory choices.

The initial hackathon result should center on three figures:

- Minimum physical qubits versus Hubbard size under the one-week cutoff, for selected times and precisions.
- Reachable size–time regions at fixed physical-qubit budgets.
- Mitigation/runtime overhead versus effective noisy volume or total rotation angle, testing whether results collapse onto a few controlling quantities.

Retain a compact table of chosen configurations and uncertainty bands around sensitive boundaries. Pin source versions, workload parameters, angle conventions, and reference-model provenance so the estimates can be reproduced. The initial deliverable is an interpretable resource study; a claim of quantum advantage over the best classical methods would require a separate comparison.

**References and their roles**

- [Original STAR proposal](../../refs/PRXQuantum.5.010337.pdf): error-corrected Clifford operations with analog rotations. [DOI](https://doi.org/10.1103/PRXQuantum.5.010337).
- [Improved STAR protocols](../../refs/PhysRevX.15.021057.pdf): small-angle error behavior and total-angle mitigation cost. [DOI](https://doi.org/10.1103/PhysRevX.15.021057).
- [TopQAD partial-FT study](../../refs/2603.13093v1.pdf): preparation/growth models and partial/full-FT comparison. Its Hubbard application estimates ground-state energies; rebuild the dynamics task rather than transferring its runtime. [arXiv](https://arxiv.org/abs/2603.13093).
- [CUDA-Q Logical paper](<refs/CUDA-Q Logical Retargetable Compilation for Fault Tolerant Quantum Computing.pdf>): compilation-linked resource analysis. [arXiv](https://arxiv.org/abs/2609.13388).
- [Executed Hubbard dynamics experiment](../../refs/2510.26300v2.pdf): experimental reference workload and compilation/mitigation ideas. [arXiv v2](https://arxiv.org/abs/2510.26300v2).
- [Locality and Error Mitigation of Quantum Circuits](https://arxiv.org/abs/2303.06496): observable-dependent mitigation cost.
- [Fundamental limits of quantum error mitigation](https://www.nature.com/articles/s41534-022-00618-z): limitations under explicit noise and protocol assumptions.

**Main open scientific questions**

Where do observable-aware mitigation and post-selection cease to buy useful dynamics at fixed hardware size? Does STAR's total-angle advantage survive realistic preparation, growth, and scheduling costs at `p = 1e-3`? Does FT become favorable first through reliable long-time evolution, through a better simulation algorithm, or through a precision requirement? Which crossovers persist when we change the benchmark details while holding the scientific accuracy criterion fixed?
