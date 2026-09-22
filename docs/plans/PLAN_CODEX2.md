**Research plan: protection, mitigation, and algorithm choice for Hubbard dynamics**

2026-09-22. Synthesis of `PLAN_CODEX.md`, `PLAN_CLAUDE.md`, their technical comparison, and the subsequent workshop discussion. This document records a proposed research program; no resource curves or crossover locations have yet been established.

**Purpose and scientific ambition**

The workshop brings together experts in applications, mitigation, partial fault tolerance, and full fault tolerance. The opportunity is to combine those perspectives to identify what governs the resource cost of useful dynamics and to discover better ways of spending a common physical-qubit budget.

The study should explain when errors can be tolerated or mitigated, when information must be protected, and how application structure changes that balance. Fermi–Hubbard dynamics supplies a concrete setting in which these mechanisms can be calculated and tested. A successful outcome would include resource predictions, an explanation of their scaling, and specific algorithmic or architectural changes that improve the accessible physics.

Three approaches provide reference points:

1. Physical-qubit computation with error suppression, post-selection, and mitigation.
2. Surface-code STAR computation with protected Clifford operations and analog rotations.
3. Fully fault-tolerant surface-code computation with reliable non-Clifford operations.

Their ordering and the existence of a useful STAR window are open questions. The study should also allow intermediate choices suggested by the analysis: different protection levels for different operations, hybrid analog/distilled rotations, and observable-aware allocation of correction and mitigation. Such extensions should follow an identified mechanism and have their implementation costs included.

The scope of a calculation is determined by its scientific value. Detailed compilation, an experimental reconstruction, or an alternative simulation algorithm belongs in the study when it resolves an uncertainty that affects the result. Software effort and elapsed workshop time do not define the scientific question.

**Questions that organize the collaboration**

| Question | Quantity or comparison that can resolve it |
| --- | --- |
| How much of the simulated system must remain accurate for the requested answer? | Observable support, physical propagation, truncation error, and the effective region affected by relevant noise |
| Where does error correction buy the most useful computation? | Reduction in estimator cost versus additional area, execution time, and lost replica throughput |
| What limits STAR for dynamics? | Total analog rotation angle, angle-independent residual errors, preparation acceptance, and resource-state delivery rate |
| What makes full FT favorable? | Protected circuit volume, synthesis and factory demand, longer accessible evolution, or a change of simulation/estimation algorithm |
| Which apparent scaling laws survive changes in implementation? | Controlled changes of fermion mapping, circuit ordering, boundary conditions, observable, and error model |
| Can a tailored intermediate strategy improve on all three reference approaches? | Optimize selected protection and rotation choices while preserving the same physical task and error target |

The experts' contributions meet at the interfaces: the application determines which errors matter; mitigation determines the attainable bias–variance tradeoff; QEC determines the cost of suppressing particular errors; compilation determines whether the assumed operations and parallelism can be realized. These interfaces are likely to contain the most consequential assumptions.

**Common device and application variables**

The stipulated physical device is a square two-dimensional nearest-neighbor array with gate time `t_g = 10 ns`, measurement time `t_meas = 100 ns`, nominal physical error probability `p = 1e-3`, and a one-week answer-time cutoff `B = 604800 s`.

Use distinct variables throughout:

- `N = L^2`: Hubbard sites in the main square-lattice family.
- `Q = M^2`: physical qubits in the device.
- `2N`: occupation modes in a direct spinful encoding, before ancillary or QEC overhead.
- `tau = J t / hbar`: dimensionless simulated evolution time.
- `T_answer`: physical elapsed time to obtain the specified answer.
- `epsilon`: absolute observable-error target; `delta`: allowed failure probability of the estimator.

At fixed `Q`, optimize allocation to data, routing, preparation, storage, factories, and independent replicas. A larger physical-qubit array can increase sampling throughput even when a single unencoded circuit already fits. For independent identical copies with no shared bottleneck, the replica count is approximately `floor(Q / Q_copy)`, subject to an actual placement.

Complete the hardware model with a common native gate set and explicit conventions for gate errors, state preparation, readout, reset, idle noise, and feedback. The scalar `p` alone does not specify these channels. A convenient stochastic reference model can be used first, with coherent or correlated components introduced where they could change the inference. Preparation and measurement errors should not silently inherit an unrelated gate-error convention.

Derive the syndrome-round duration from a realizable schedule using the stated physical times. Measurement, reset, idle exposure, and feedback are potentially important at `t_meas / t_g = 10`. Use consistent surface-code suppression assumptions in STAR and full FT; differences must be attributable to actual protocols or schedules.

Sweep `p` around the stipulated point, including `1e-4` as a useful reference. Change suppression-law assumptions consistently across encoded approaches. The one-week cutoff is an operational contour through the resource curves; retain the underlying curves so its role is visible.

**Physical task and benchmark families**

The agreed starting task is a simple initial state, a local density or doublon observable, and absolute error `0.01`, with lattice size and evolution time swept. Adopt the following concrete baseline:

- Standard spinful nearest-neighbor Hubbard Hamiltonian:
  `H = -J sum_<ij>,sigma (c†_i,sigma c_j,sigma + h.c.) + U sum_i n_i,up n_i,down`.
- Even `L x L` squares with open boundaries and zero flux.
- Half-filled Néel occupation product state, initially at `U/J = 4`.
- Local doublon observable `D_i(t) = <n_i,up n_i,down>` at a specified central site.
- Absolute error `epsilon = 0.01` at 95% confidence for a specified time point.

Doublon density gives nontrivial dynamics in this setting. Total local charge density can be pinned by the particle-hole and spin symmetries of the baseline quench and is therefore a poor default diagnostic.

Distinguish three scientifically different size/time questions:

1. **Fixed-time local physics.** For a local Hamiltonian and the chosen simple initial state, a sufficiently large region around the observable may suffice at fixed accuracy. Implement and budget such a truncation if used. Optimized cost can then saturate with `L`; a full-lattice circuit does not acquire this saving automatically.
2. **Dynamics on the scale of the system.** Increase evolution time with linear size, initially considering a propagation-time family `tau = c L`. Other time scalings should be tied to the transport process being studied. Here the relevant region grows and a size-dependent resource comparison can remain meaningful.
3. **A larger information demand.** If needed to test a mechanism, consider a separated correlator, spatial profile, or extended observable. Define its normalization and error criterion explicitly. Spatial averages, one local value, and a whole profile have different statistical costs and may respond differently to noise.

Keep the scientific task fixed when comparing approaches. Any locality reduction, symmetry, or measurement reuse used to improve an answer should be available wherever it can be implemented. State clearly whether a curve represents the straightforward full-system implementation or an optimized implementation exploiting that structure.

Initial-state extensions should test an identified dependence. A short-depth dimer preparation with a holon–doublon disturbance is particularly useful for following charge and spin propagation. Ground-state preparation is a separate resource problem and should be included only if the selected scientific task requires it.

**Accuracy and the object being costed**

Cost the production of an observable estimate. Separate contributions from state preparation, locality truncation, simulation approximation, gate synthesis, residual hardware/mitigation bias, and statistical uncertainty. Assign a total error budget and optimize its allocation. The confidence statement must include the relevant systematic uncertainties; statistical error bars or smoothing confidence bands alone are insufficient.

For logical faults, bound their contribution to the requested observable error. Requiring zero faults throughout every shot in an entire campaign is a different and generally stronger criterion. Use one consistent convention for the event or time interval to which a reported logical error probability applies.

For a curve or several observables, specify pointwise versus simultaneous coverage and count the circuits, training data, and shared measurements accordingly. In the initial matched comparison, all approaches use ordinary sampling. If precision becomes a main driver, coherent estimation can be a separately costed algorithmic extension, including controlled operations, state preparation/inversion, and additional ancillas.

The principal resource outputs are

`Q_min(N, tau, epsilon, delta; T_answer <= B)`

and

`tau_max(N, Q, epsilon, delta; T_answer <= B)`.

Retain the area/runtime Pareto frontier and its bottlenecks. An area-time product alone hides the distinction between sample throughput, circuit latency, and resource-state supply.

Where acceptance and estimator variance factorize, a useful model is

`T_answer ≈ N_0(epsilon_stat, delta) * Gamma^2 * t_attempt / (a * m)`.

`N_0` is the ideal sample requirement, `Gamma^2` is the variance multiplier, `a` is whole-circuit acceptance, and `m` is the replica count. Use the actual conditional estimator when post-selection and mitigation interact. State-factory rejection and RUS retries belong in their supply and execution model; charging them again as whole-circuit rejection would double count the cost. Include setup, calibration, training, and online classical latency where they contribute materially, and state the accounting boundary.

**Candidate scaling variables and what they would mean**

Start with observable-relevant physical noise exposure, total analog rotation exposure, accumulated residual floors, and protected logical volume. Test whether these variables organize results across workloads and implementations; do not assume a universal collapse.

For a specified PEC construction under suitable stochastic noise, a reference form is

`log(Gamma_NISQ^2) ~ sum_j c_j p_j = chi_phys`.

The included locations and coefficients depend on the estimator, circuit, and noise model. Full-circuit PEC and locality-aware mitigation should be distinguishable. A general error-mitigation lower bound only applies when its assumptions match the task. Symmetry post-selection changes acceptance and conditional errors; acceptance is not an error-free probability.

Under this exponential model, with `S_1` attempted shots per replica and `m` replicas, the tolerable exposure scales schematically as `chi_phys <= log(m S_1 / N_0)`. This explains a logarithmic improvement with runtime or replicated hardware once the mitigation term dominates. It does not supply a universal gate ceiling, and `S_1` itself depends on circuit depth and scheduling.

For STAR, use the effective error of a completed rotation, including its RUS sequence. Under a phase-flip model with error `q_j`, the PEC factor is

`Gamma_STAR^2 = product_j (1 - 2 q_j)^(-2) ≈ exp(4 sum_j q_j)`.

This expression presumes the modeled channels and cancellation operations; other errors require their own treatment. Improved small-angle protocols motivate

`q_j ≈ alpha_j p |theta_j| + q_floor,j`,

and hence the candidate exposures

`chi_angle = 4 p sum_j alpha_j |theta_j|`,

`chi_floor = 4 sum_j q_floor,j`.

The coefficients, floors, acceptance probabilities, and preparation times vary with the selected protocol and code distances. A model with fixed `alpha` and fixed `q_floor` is a useful slice through this optimization, not an intrinsic property of STAR.

Use the convention `U_P(theta) = exp(-i theta P)`. For ordinary evolution and a straightforward positive-time product formula, `Theta = sum_j |theta_j| ≈ (t/hbar) sum_l |h_l|`, with `H = sum_l h_l P_l`, before Clifford extraction and other circuit simplifications. First- and second-order positive-time formulas can preserve this total angle as the step count changes. Higher-order signed stages, different decompositions, and optimizations require recounting.

Thus finer Trotter steps can leave the leading angle-dependent exponent unchanged while increasing Clifford work, preparation demand, execution time, and accumulated floors. Additional precision can also increase ordinary sampling cost. Derive any factor of two from the actual gate convention and circuit: phase-estimation constructions using controlled directions can have different angle budgets from ordinary dynamics. Translating notation must leave the physical prediction invariant.

For FT, use a suppression model such as `p_L(d) = A (p/p_th)^((d+1)/2)` in its calibrated regime and interval convention. Below threshold, the required distance can grow logarithmically with protected volume and inverse error budget, with patch area roughly proportional to `d^2`. Logical schedules, routing, and state supply set the accompanying time overhead. Independent deterministic rotation synthesis gives approximately `N_rot log(N_rot / epsilon_synth,total)` T gates; actual counts depend on angle classes and synthesis strategy.

These mechanisms give testable questions: does STAR fail first from total angle, a residual floor, or poor preparation throughput? Does FT gain most from reliable longer evolution, precision, or a different algorithm? Does physical-qubit computation remain competitive through observable structure even after a full-circuit estimate predicts failure?

**Published constants as useful starting cases**

Populate the model immediately with named reference cases, keeping their conventions and validity ranges visible. Replace or refine a reference case when the implemented protocol differs.

| Reference case | Starting information | How to use it |
| --- | --- | --- |
| Original STAR injection [R1] | Leading phase-error probability `2p/15` per injection attempt; approximately `4p/15` per completed rotation after RUS | Reference exponent `log(Gamma^2) ≈ 16 p N_rot / 15`; higher-order errors and the source noise model remain relevant |
| Original compact-block timing [R1] | Approximately 18 logical clocks per analog Pauli rotation on average; a clock is `d` syndrome rounds | A particular layout/consumption model, to be compared with the selected placement and parallel schedule |
| Original surface-code fit [R1] | `C_Z = 0.0679`, `p_th,Z = 0.00385`; `C_X = 0.0819`, `p_th,X = 0.00416`, with reported uncertainties | Fit from the paper's circuit-level experiment; retain its measurement interval and channel convention when translating to rates |
| Common illustrative surface-code model | `A = 0.1`, `p_th = 0.01` | Sensitivity case, used consistently across encoded approaches; not a hardware calibration |
| Improved STAR [R2, R3] | Effective small-angle error approximately `alpha_RUS p |theta|`; a fixed-weight protocol gives a coefficient of order `0.4 k` in the source convention | Translate the angle convention and choose the appropriate protocol, rather than assigning one coefficient to all gadgets |
| Growth model [R3] | Effective rotation error includes `2 I_growth(d_i, d_f)` in addition to its angle-dependent part | Optimize initial/final distances jointly with acceptance, area, and throughput |
| TopQAD comparison [R3] | A favorable STAR region around `10^5–10^6` small-angle rotations under its chosen hardware, width, and error targets | A reference result to explain or reproduce under matching assumptions; not our crossover prediction |

At `p = 1e-3`, the original STAR model suggests a limited per-rotation improvement and motivates the hypothesis of a narrow advantage window. Whether that hypothesis survives protected routing, the ratio of physical gates to rotations, improved injection, and area allocation is part of the calculation. The surface-code fit can strongly affect area, but which parameter matters most should be established through sensitivity analysis.

Count the actual circuit. For an open `L x L` lattice, `E = 2L(L-1)` spatial bonds. The uncombined first-order Pauli decomposition has `4E + 3N = 11N - 8L` nonidentity terms. Second-order ordering and merging change executed multiplicities. In a fixed particle-number sector, some onsite terms can be removed as a sector-dependent phase. Record these choices before using a quoted `11N` count or a Hamiltonian coefficient norm. Physical entangling-gate and idle counts require the fermionic mapping and schedule as well.

**Algorithm–architecture comparison and routes to improvement**

First produce a matched second-order product-formula comparison. Share the Hamiltonian, state, observable, target error, and simulation method. A fixed-circuit slice isolates implementation differences; subsequently allow each approach to optimize the step count and error allocation within the same scientific target.

Compare a direct Jordan–Wigner implementation with at least one locality-conscious fermionic compilation. Include any extra encoding qubits, initialization, constraints, and measurements. Keep routing and idle exposure explicit. Mapping-dependent changes in gate count or depth should be distinguished from the cost of protection itself.

Determine product-formula error for the observable. Locality can remove explicit size dependence at fixed time, but does not by itself establish `r proportional to tau^(3/2) / sqrt(epsilon)`: time-dependent coefficients and support growth remain. Use exact small-system checks, available analytical bounds, and larger-system convergence evidence to justify the extrapolation actually used. Identify empirical calibrations separately from bounds.

Then compare optimized algorithm choices. Qubitized dynamics is a natural full-FT candidate; include PREPARE/SELECT, ancillas, synthesis, routing, and any success overhead. Larger-angle rotations can make it unattractive for STAR, but do not make it incompatible. Evaluate further combinations when a plausible cost argument warrants them. Report matched and optimized comparisons separately so that an algorithm change has an identifiable contribution.

Use the first results to propose targeted improvements. Candidate directions include assigning different rotation classes to analog or distilled implementations, improving small-patch preparation and growth, combining modest correction with mitigation, and optimizing observable support before physical compilation. For fixed ideal sample requirements and a valid factorized cost model, an elementary decision criterion is whether

`Delta log(t_attempt) + Delta log(Gamma^2) - Delta log(a) - Delta log(m) < 0`.

This makes the marginal benefit of additional protection or a better gadget explicit. Include conversion, routing, and synchronization costs when testing a hybrid proposal. A promising extension should predict a measurable shift in the resource frontier and identify its cause.

**Role of the demonstrated Hubbard experiment**

The Phasecraft/Quantinuum experiment [R5] supplies an executed dynamics workload: a `7 x 4` doubly periodic nearest-neighbor lattice with threaded pi flux, an entangled dimer state with holon/doublon defects, `U/J = 0, 4`, and evolution through `Jt = 2`. It uses four second-order Trotter steps and at most 2,415 two-qubit gates, with training-based mitigation and Gaussian-process smoothing.

Use it to ask concrete questions about observable robustness, compilation efficiency, and how much apparent performance comes from postprocessing and exploitable physical structure. Its effective all-to-all ion-trap connectivity, native entanglers, and relabeling of some SWAPs require recompilation for the stipulated device. Its demonstrated execution does not establish the chosen `0.01` total-error target.

Reconstruct the exact instance if doing so resolves one of those questions. A related square-lattice dimer/defect quench can test transferability while separating geometry and initial-state effects. Growing only the length at fixed width four does not produce the same two-dimensional scaling family as growing both dimensions; flux and boundary conditions also change the physics.

The reconstruction is a scientifically motivated workstream, with a defined question and owner when undertaken. Its priority should follow the information it provides about the central mechanisms.

**Calculation strategy and tool roles**

Build an inspectable analytical model that accepts workload summaries and architecture models independently of any compiler. Store rotation angles and multiplicities, dependencies, physical placements where available, error-budget choices, and the provenance of fitted constants. Let the model return feasible configurations and bottlenecks, not only a single cost number.

Use tools where they make a calculation more reliable or reveal a missing cost:

- **CUDA-Q Logical:** candidate compilation and scheduling framework in `~/soft/cuda-quantum/preview/logical`. The square-lattice helper currently implements first-order controlled evolution and defaults to periodic boundaries; adapt it to the actual task. The short Hubbard example is a chain. The built-in physical surface target is a reference model with native Clifford actions and multi-Pauli measurements; provide a nearest-neighbor realization or an explicitly separate cost layer. The inspected synthesis pass uses GridSynth. RUS synthesis formulas describe a different implementation and require their associated scheduling costs.
- **pyLIQTR:** use its Hubbard block encodings and `qsvt_dynamics` for the qubitization comparison. Check boundary conditions and Hamiltonian conventions. In the inspected `FermiHubbardSquare` implementation, the internal variable `N` counts modes: with our site count its normalization is `alpha = 4N|J| + N|U|`. Derive the normalization for the actual encoding used, particularly if the Hamiltonian is shifted or restricted.
- **Q# resource estimation:** a useful independent FT calculation and an initial route while more detailed compilation is being prepared. Compare error budgets, connectivity, synthesis, and factory assumptions before interpreting agreement or disagreement.
- **General circuit simulation:** use exact statevector, trajectories, or density matrices for small non-Clifford Hubbard instances and mitigation checks. Ordinary Stim is appropriate for stabilizer/QEC subproblems, not generic Hubbard rotations without an additional method.

As of this discussion, package metadata was checked in `~/soft/mambaforge/envs/qb_gsee`: `pyLIQTR 1.4.2`, `qualtran 0.4.0`, `cirq-core 1.4.0`, and `qsharp 0.28.302812`. This establishes installed versions, not a validated end-to-end workflow. Claude reported CUDA-Q was unbuilt/uninstalled in its inspected environments. Verify the particular environment before using it; a build need not hold up the analytical work.

Tool compatibility is an implementation question. Agreement between packages using the same assumptions is a consistency check; validation of a physical model requires separate evidence.

**Research sequence and evidence needed**

1. **Agree on the physical question and comparison conventions.** Fix the baseline observable and estimator, noise reference case, angle convention, and resource boundary. Select the size/time families that distinguish finite-time local physics from dynamics on the scale of the system.
2. **Derive and populate the reference models.** Implement both STAR generations, a specified physical-qubit mitigation menu, and an FT reference. Use published constants with explicit assumptions. Identify candidate crossover mechanisms and parameters to which the conclusion is most sensitive.
3. **Construct the common workloads and alternative mappings.** Verify Hamiltonian and state conventions, observable dynamics, angle totals, term multiplicities, and error-budget allocation. Produce fixed-circuit and optimized-step comparisons.
4. **Resolve consequential uncertainties.** Calibrate the mitigation bias–variance tradeoff, representative routed schedules, syndrome timing, STAR state preparation/growth, and FT state supply. Choose analytical derivations, simulation, compilation, or experimental comparison according to the uncertainty being resolved.
5. **Optimize and challenge the apparent regimes.** Vary hardware quality, precision, observable support, evolution time, and algorithm. Test a promising hybrid if the reference analysis identifies a specific opportunity. Revisit any assumed law that fails to organize the data.
6. **Extract the general result.** State the controlling variables, applicable assumptions, and improvements supported by the calculations. Explain both stable conclusions and boundaries that remain sensitive to an unresolved model choice.

These activities can proceed concurrently where their assumptions are compatible. The application, NISQ, STAR, and FT experts should review the interfaces where their models are combined. Avoid expanding every uncertainty into a separate project; resolve it to the level needed to assess the scientific conclusion.

**Outputs that would make the synthesis useful**

- A resource frontier: minimum physical qubits versus Hubbard size under the one-week cutoff, at selected times and precisions.
- Reachable `(L, tau)` regions at fixed `Q`, with implementation and algorithm choices identified.
- Scaling plots against observable-relevant noise exposure, total rotation angle, accumulated floors, and protected volume, showing where a common description succeeds or fails.
- A decomposition of the resource cost near each crossover: sampling, rejection, routing, protection, synthesis, and state supply. Include sensitivity to the assumptions that can change the winner.
- At least one concrete explanation of how changing a protocol or algorithm moves a boundary, if the analysis supports such an improvement.
- Reproducible workload definitions, model parameters, source versions, and representative validation cases sufficient for another expert to reconstruct the claims.

A scientifically useful result may show that a familiar three-regime picture fails for a local task, that a different observable restores it, or that an intermediate strategy improves on the reference choices. Resource estimates establish performance under their declared models. Claims about superiority over classical simulation require a separate task-matched classical comparison.

**References used as inputs**

- **R1:** [Akahoshi et al., original STAR proposal](../../refs/PRXQuantum.5.010337.pdf), PRX Quantum 5, 010337 (2024). Injection, RUS, layouts, and circuit-level surface-code fits.
- **R2:** [Toshio et al., improved STAR protocols](../../refs/PhysRevX.15.021057.pdf), Physical Review X 15, 021057 (2025). Small-angle errors and total-angle mitigation cost.
- **R3:** [Chung et al., partial-FT resource study](../../refs/2603.13093v1.pdf), arXiv:2603.13093. Preparation/growth costs and TopQAD comparisons. Its Hubbard task is ground-state energy estimation; reconstruct the appropriate dynamics costs here.
- **R4:** [CUDA-Q Logical](<refs/CUDA-Q Logical Retargetable Compilation for Fault Tolerant Quantum Computing.pdf>). Compilation-linked resource accounting and application examples.
- **R5:** [Phasecraft/Quantinuum Hubbard dynamics experiment](../../refs/2510.26300v2.pdf), arXiv:2510.26300v2. Executed dynamics, optimized circuits, and mitigation/validation methods.
- [Locality and Error Mitigation of Quantum Circuits](https://arxiv.org/abs/2303.06496). Observable-dependent mitigation costs.
- [A Theory of Trotter Error](https://arxiv.org/abs/1912.08854). Product-formula error analysis and local-observable simulation.
- [Fundamental limits of quantum error mitigation](https://www.nature.com/articles/s41534-022-00618-z). Lower bounds under specified noise and mitigation assumptions.
