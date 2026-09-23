# What the NISQ qubit estimate measures

This is the historical four-step compiled ledger. The current comparison is
in [INTEGRATED_COMPARISON.md](INTEGRATED_COMPARISON.md), and the direct
experimental data checks are in [PUBLISHED_BENCHMARKS.md](PUBLISHED_BENCHMARKS.md).
Later dynamics calculations qualify the required step count; the numbers
below retain their original fixed-step scope.

The curve discussed here is a sufficient sampling budget for gatewise probabilistic
error cancellation (PEC), using the compiled observable dependency cone. It
does not establish the resource frontier of physical quantum computation.
The Phasecraft experiments expose an important missing comparison: efficient
hardware compilation together with application-specific, trained mitigation.

## Reproducing the 4 × 4 point

The saved circuit has 32 active qubits and a 36-qubit square footprint. It
uses four plaquette steps at τ = 1 and U/J = 4, with 2,648 native Rzz gates,
964 scheduling layers, and a 9.84 microsecond attempt time. The reference
targets 0.001 absolute error (1% relative at D* = 0.1), with 0.0005 allocated
to statistical error at 95% confidence. Its ideal Hoeffding shot bound is
7,377,759.

The exponent in the PEC variance bound is

\[
\chi = 10.5927_{2q}+4.1201_{1q}+5.1681_{\mathrm{idle}}+0.1361_{\mathrm{SPAM}}.
\]

Thus Γ² = exp(χ) = 4.935 × 10⁸. The sufficient shot budget becomes
3.641 × 10¹⁵, about 1,135 years on one copy at the reference clock speed.
The week cutoff buys 59,239 replicas, each occupying 36 qubits:
Q = 2,132,604. No individual circuit occupies a million qubits.

| Same circuit, clock and precision | Sufficient qubit budget within one week |
|---|---:|
| Reference noise, gatewise PEC | 2,132,604 |
| Remove idle noise | 12,168 |
| Retain only two-qubit noise | 180 |
| Set the mitigation multiplier to one | 36 |

The last three rows are diagnostics. They do not model a demonstrated new
compiler or estimator. The last row leaves hardware bias uncorrected; it
only illustrates that the ideal sampling budget takes about 73 seconds.
At 0.01 absolute error, the full reference gives 21,348 qubits.

These counts follow a sufficient Hoeffding bound and factorized channel
inverse, not a necessary shot count. Exact channel knowledge is assumed.
Even that conditional statistical bound does not certify dynamics accuracy:
four steps fail the small-patch Trotter checks, and the 4 × 4 Trotter error
has not been calibrated.

Run `python python/plot_capacity.py` to regenerate
[the numerical decomposition](../results/presentation_v1/nisq_sampling_audit.json).
The presentation audit recomputes it against the saved circuit and current
assumptions. The underlying scientific artifacts have not been changed.

## The experimental comparison

The [Quantinuum paper](https://arxiv.org/abs/2510.26300v2) uses 56 qubits for
7 × 4 sites, four second-order steps, and TFLO plus time-series smoothing.
The [Willow paper](https://arxiv.org/abs/2510.26845v1) uses three steps and
up to 72 fermionic-mode qubits; its layout table includes two routing qubits
for the 6 × 6 circuit. Its mitigation includes TFLO, occupation-sector
selection, and smoothing. These are different estimators from gatewise PEC.

The Willow paper's Table I gives 400r + 20 CZ gates for 4 × 4 and
1440r + 52 for 6 × 6. At r = 4 these are 1,620 and 5,812, compared with
our 2,648 and 7,952 Rzz gates. These four-step experimental-compiler counts
are evaluations of the published formulas, not claims about circuits run in
that experiment. A comparison of errors and durations must respect the
different native gate implementations.

Different observables, states, precision criteria, and noise channels also
matter. Demonstrating useful dynamics does not establish the specific
0.001 absolute-error contract for our central doublon at every interacting
point. Conversely, that distinction does not justify treating the current
PEC curve as an estimate of the best physical strategy. In particular,
the Willow comparison cannot be dismissed as an all-to-all ion-trap advantage.

## What needs to change in the NISQ comparison

Retain gatewise PEC as one specified strategy. Benchmark compilation against
the published native circuits, including the compression of single-qubit
operations, scheduling, and the treatment of idle intervals. The current
snake routing and separately charged one-qubit primitives are not a
validated approximation to those compilers.

Add an observable estimator based on the experimental strategy and measure
both its variance and residual bias. For a fixed affine correction
\(\widehat D=a+b\overline X\), conditional on its fitted coefficients,

\[
\operatorname{Var}(\widehat D\mid a,b)=b^2\operatorname{Var}(X)/M.
\]

The relevant multiplier is then the fitted response and raw variance,
not automatically a product of every inverse gate channel. Training adds
cost and uncertainty; transferring a fit from U = 0 to U = 4 can leave bias.
Validation must use held-out data and include that bias in the same error
budget as the other arms. Postselection acceptance alone cannot do this.

The next useful comparison is therefore a matched experimental benchmark
and a bias–variance curve for physical estimators. No optimized NISQ
capacity or revised architecture crossover is established by this audit.
