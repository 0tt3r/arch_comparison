**QUOPS specification — Proctor et al., arXiv:2609.12146v1.** References below use printed/PDF pages, numbered from one. The [numeric specification](../results/quops_spec_v1/spec.json) preserves source references, PDF/code hashes, all tabulated results, and unresolved fields as `null`. Measurements, simulations, projections and historical conversions are distinguished.

**Circuit distribution.** On a width-$w$ computational register, each benchmark layer applies **CNOTs first, rotations second**. Independently sample a uniform matching of the all-to-all qubit set: $\lfloor w/2\rfloor$ disjoint CNOTs, leaving a uniformly random qubit unmatched for odd $w$. Each qubit receives $R_P(\theta)=e^{-i\theta P/2}$, with independent uniform $P\in\{X,Y,Z\}$ and $\theta\in[0,2\pi)$. SI VI.A.4 describes ordered control–target pairs. These are benchmark connections; physical routing costs remain architecture-dependent. [SI II.A–B, pp.12–13, Eq.8, Fig.4](../refs/2609.12146v1.pdf#page=13).

There are $m$ benchmark layers, individual-layer depth $d=2m$. The final layer acts on a uniformly selected subset of $k\in\{1,\ldots,w\}$ qubits, using the same recipe. With $B(n)=n+2\lfloor n/2\rfloor$,

$$
N_R=w(m-1)+k,\quad N_C=\lfloor w/2\rfloor(m-1)+\lfloor k/2\rfloor,
\qquad s=N_R+2N_C=B(w)(m-1)+B(k).
$$

Thus sizes are discrete; idles contribute zero. For dense even-width circuits, $s=wd=2wm$. The counts $N_R=s/2$, $N_C=s/4$ also hold for partial layers when both $w$ and $k$ are even. The inclusive score cone is $w^2\le s\le w^3$, equivalently $w\le d\le w^2$ only in that dense case. Approximate compilation is allowed but its error counts against the original unitary; hardware depth can differ. Data-qubit reuse is prohibited. [SI I.A, II.B–C, III.D, VI.B.4; Eqs.4,9–10,19–21; pp.12–16,31](../refs/2609.12146v1.pdf#page=13).

**Success and measurement.** Define entanglement/process fidelity $F_e=\langle\Phi|[(\Lambda\circ\mathcal U^\dagger)\otimes I](|\Phi\rangle\langle\Phi|)|\Phi\rangle$, with normalized maximally entangled $|\Phi\rangle$. Then
$\gamma=(4^wF_e-1)/(4^w-1)$, and the target is $\bar\gamma_{w,s}=\mathbb E_C\gamma(U(C),\Lambda(g(C)))$. It excludes initial preparation/final readout of computational qubits, but includes intermediate measurements and QEC. A global depolarizing channel's polarization equals its mixing parameter, not its process fidelity. [SI I.B, III.B; Eqs.5–6,15; pp.12,15](../refs/2609.12146v1.pdf#page=12).

Mirror-circuit fidelity estimation uses, in right-to-left execution notation,
$M_1=L'g_{ref}(c^\dagger)g(c)L$,
$M_2=L'g_{ref}(c^\dagger)g_{ref}(c)L$, and $M_3=L'L$.
Random local Clifford wrappers and Pauli frames produce known output strings. From Hamming-distance frequencies $h_j$,
$\hat\lambda=[4^w\sum_j(-1/2)^jh_j-1]/(4^w-1)$.
Pooling gives $\hat\gamma=\bar\lambda_1/\sqrt{\bar\lambda_2\bar\lambda_3}$; if $g=g_{ref}$, one mirror ensemble suffices and $\hat\gamma=\sqrt{\bar\lambda_1/\bar\lambda_3}$. Forward-plus-inverse approximately doubles circuit size; wrappers/reference costs remain. The estimator assumes suitable stochastic reference errors and factorization, and can have small optimistic bias. [SI IV.A–C, pp.17–22, Eqs.34–41,53–59](../refs/2609.12146v1.pdf#page=18).

There is **no universal circuit/shot count**. Willow's region used 100 QUOPS circuits → 200 mirrors × 1000 shots, plus 200 SPAM circuits/width; its rate run used 25 QUOPS circuits × 500,000 shots. Boston used 100 QUOPS circuits/shape, 1000 shots and 200 references/width. H2 generally used 10 mirrors × 200 shots; physical Helios used 20 × 200. Both used 2000 SPAM shots. Logical Helios varied from 100 to 1000 total shots; JSON records every case. [SI VII.D–IX, pp.39–43](../refs/2609.12146v1.pdf#page=39).

Pass when the **one-sided 95% lower bound** satisfies
$\hat\gamma-1.6448536\sigma\ge\alpha$, normally $\alpha=e^{-1/2}=0.60653066$. Bootstrap circuits, shots and SPAM references; the fifth bootstrap percentile is an alternative. The paper uses a preselected increasing-size sequence, stopping at the first failure; any procedure controlling score overestimation at 5% is allowed. Point-estimate scans do not certify scores. Additional region tests use multiplicity control; the combined region guarantee is 90%, with monotonicity required for downward closure. [SI IV.C–D, pp.22–24, Eqs.60–71](../refs/2609.12146v1.pdf#page=23).

**Score and rate.** $Q=\max\{s:\bar\gamma\ge\alpha,\ w^2\le s\le w^3\}$, restricted to admissible shapes. Report its width and rate. With accepted-forward time
$\tau=\operatorname{median}[\tau_{accept}+(p_{reject}/p_{accept})\tau_{reject}]$,
$\Omega_{raw}=s/\tau$, $\Omega=s\bar\gamma^2/\tau$.
Operationally,

$$
\Omega=\frac{2s\hat\gamma^2\sqrt{N_{total}N_{kept}}}{\tau_{wall}}.
$$

$N_{total}$ and $N_{kept}$ count attempted and accepted mirror executions. Include unavoidable circuit switching, calibration and reset; separable queue/transfer delays and SPAM-reference time may be excluded. Postselection's square-root correction assumes comparable rejection in both mirror halves. Lowering $\alpha$ defines error-mitigated capability with approximate overhead $1/\alpha^2$; it does not prescribe extra benchmark shots or demonstrate successful mitigation. [SI III.E–F, IV.E, pp.16–17,25–26, Eqs.22–31,72–74](../refs/2609.12146v1.pdf#page=25).

**Measured hardware.** All headline rows from [Table I, p.3](../refs/2609.12146v1.pdf#page=3):

| Platform / implementation | Q | QUOPS/s | Width |
|---|---:|---:|---:|
| Willow physical | 216 | 20,000,000 | 6 |
| Boston physical | 204 | 310,000 | 6 |
| H2-1 physical | 1320 | 353 | 12 |
| H2-1 leakage-postselected | 1392 | 231 | 12 |
| Helios-1 physical | 1504 | 303 | 16 |
| Helios-1 leakage-postselected | 1824 | 247 | 16 |
| Helios-1 logical [[7,1,3]] | 40 | 4.9 | 4 |

JSON includes all 28 rows of Tables IV–VI. Logical v3 has $\hat\gamma=0.730$, lower bound 0.699, 500 shots and 73 minutes. Three-bit nearest-angle synthesis averages 0.5 T/rotation; injected depolarizing errors charge synthesis infidelity. Only v3 ensures composable distance-three fault tolerance. [SI VII.C, pp.34–37; Table V, p.41](../refs/2609.12146v1.pdf#page=34).

At $\alpha=.01$, measured-plus-extrapolated $(Q,\Omega)$ are Willow (651,12000), Boston (899,91), H2 (7796,0.09), Helios (20442,0.29); these are not all directly demonstrated shapes. Willow/Boston per-shape polarizations appear graphically in Fig.17, without numeric tables. JSON separately includes 32 historical QV conversions. [pp.3,43–46](../refs/2609.12146v1.pdf#page=43).

**FT projections and challenge targets.** The planar model uses rotated surface codes, lattice surgery, cultivated T states, uniform depolarizing $p=10^{-3}$, 1 μs syndrome cycles and 10 μs reaction time. It scans odd code distances 3–21. GridSynth precision is
$\epsilon=\sqrt{(1-e^{-1/2})/(3N_R)}$, with quadratic synthesis charge and equal synthesis/Clifford/magic budgets. It assumes spacetime-limited $L=V/A$ logical steps. Its success proxy is
$(1-\epsilon^2)^{N_R}(1-p_{cyc})^{d_{code}V_{cliff}}(1-p_{mag})^{N_T}$.
Table III gives physical-qubit/score/width triples (5000,1950,13), (10000,9240,21), (25000,68470,41), (50000,273632,68). [SI VI.A, pp.27–29, Eqs.84–92](../refs/2609.12146v1.pdf#page=27).

The separate all-to-all forecast uses surface distances 3/5/7/9 or [[20,2,6]], $p_{2q}=p_{SPAM}=10^{-4}$, $p_{1q}=p_{memory}=10^{-5}$, 10 ms cycles, transversal/fold-transversal operations, cultivation and correlated decoding. Dynamic scheduling/fitted errors replace a single distance law; synthesis retains one-third budget. For 1000/5000/10000 physical qubits, scores are approximately 67,000/370,000/2,200,000. [SI VI.B, pp.29–33; Fig.3b](../refs/2609.12146v1.pdf#page=29).

Challenge conversion instead uses mixed-fallback synthesis
$n_T(\epsilon)=0.53\log_2(1/\epsilon)+4.86$, four T/Toffoli, and ignores Clifford/idle errors. Optimize $Q_{target}=1/(n_Tp_T+\epsilon)$, with $p_T\approx1/(2N_T)$. [SI V, Eqs.75–83; Table II, p.27](../refs/2609.12146v1.pdf#page=27).

| Problem | Width | Toffolis/execution | Target Q | Required QUOPS/s |
|---|---:|---:|---:|---:|
| RSA-2048 | 1399 | 7.1×10⁸ | 2.5×10⁸ | 5700 |
| FeMoco, 76 orbitals | 1459 | 1.0×10⁹ | 3.4×10⁸ | 800 |

RSA allows 12.07 hours/run and 9.2 expected runs; FeMoco assumes five days/run. Target rates are raw, unlike polarization-attenuated reported rates.

**Reference-code comparison.** [quops.py](../crossover_martin/crossover/benchmarks/quops.py), with full line references and 27 findings in JSON:

| Assumption | Assessment |
|---|---|
| Gate weights, default threshold/cone, maximizing size | Match. |
| Full-layer floor(w/2) CNOTs; depth 2D | Match represented layers only. |
| D=round(s/(2w)); no partial layer | Conflicts with exact sizes; (6,40) becomes size 36. |
| Random family/order | Counts omit sampling; paper starts with CNOT. |
| Architecture F≥α | Conditional polarization proxy; lacks confidence certification. |
| 30×100 shots, one-week cap | Added modeling choices. |
| Twofold mirror cost | Approximation; reference/SPAM/overhead costs omitted. |
| F² rate attenuation | Matches idealized rate; copies require justified scaling. |
| Postselection | Acceptance counts/time absent. |
| Extra 1/F² benchmark shots | Not prescribed by paper. |
| Width 3–4000; 40 sizes; 15% width steps | Restricted, potentially incomplete search. |
| First size failure, infeasible width, three dead widths | Monotonicity assumptions/heuristics; not certification. |
| Downstream synthesis/FT defaults | Different error convention, budgets, factories and logical-error law. |
| Downstream NISQ readout error in F | Conflicts with terminal-SPAM exclusion. |

**Unresolved details.** SI IV.B reverses layer order for its derivation; depth notation varies. Planar projections approximate $N_R=s/2$ at odd widths and appear to omit a code-distance factor in time. Logical-error inputs and fitted scheduler parameters are not fully tabulated. IBM rate prose conflicts with its squared-polarization formula. H2's (56,224) postselected interval is malformed, and SPAM-reuse descriptions differ. JSON preserves these issues rather than resolving them by assumption.

The paper promises a later data release ([p.8](../refs/2609.12146v1.pdf#page=8)); a bounded check of the cited [PRAQTICE tree](https://github.com/sandialabs/PRAQTICE/tree/1e73c2d13699f4b60bf133b197b1230b98f8d98d) found no identifiable QUOPS data package.
