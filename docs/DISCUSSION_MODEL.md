# Back-of-the-envelope+ for the team discussion

Start with the [capacity figure](../figures/discussion_capacity.pdf), the
[generated numerical summary](../results/discussion_v3/summary.md), and the
[time-calibration figure](../figures/time_calibration.pdf). The
[0.01 absolute-error figure](../figures/discussion_capacity_loose.pdf) is the
precision sensitivity. The current model optimizes each arm's error allocation
and resource-state supply using measured finite-patch Trotter coefficients.

The hardware and target are unchanged: p = 0.001, 10 ns gates, 100 ns
measurement/reset, 240 ns syndrome cycles, U/J = 4, and one week. Strict
precision is 0.001 absolute, or 1% relative at the reference doublon D* = 0.1.
The comparison concerns a single local observable after evolving the full
lattice. Both [earlier](../results/discussion_v1/summary.md)
[discussion versions](../results/discussion_v2/summary.md) retain their data
and original figures.

The [delivery sensitivity](DELIVERY_AND_LARGE_PATCH.md) adds explicit compact,
intermediate, and fast data blocks, reoptimizes source supply, and varies
analog parallelism. It preserves the curves above as the reference.

## Time dependence and algorithm choice

We measured the local error for **both termwise and plaquette Strang orderings**
on 2 × 2 and 3 × 3 patches at 18 times through tau = 8. The sector evolution
preserves the order in the corresponding compiled circuits; independent tests
compare it with the full qubit circuit and sparse Hamiltonian evolution.
At each time, the starting count is ceil[16 max(1,tau)], followed by doubling.
The two finest estimates of |error| r² must agree within 5%; one further
doubling is used if needed. Signed errors and all step counts are retained.

On 3 × 3, the coefficients are:

| tau | Plaquette | Termwise |
|---:|---:|---:|
| 1 | 0.218 | 0.0767 |
| 2 | 0.926 | 0.381 |
| 4 | 5.27 | 1.96 |
| 8 | 7.67 | 23.74 |

The r⁻² law persists, while the coefficient oscillates and depends on ordering.
Additional checks on both small squares at r = ceil(8 tau) and twice that value capture
finite-step corrections at the optimizer's lower limit. The coefficient used
for planning takes the maximum across these checks and the asymptotic pair.
This was needed when a direct loose-precision check exceeded the asymptotic
estimate plus its 10% margin.
The planning model uses

\[
C_{\mathcal A}^{\rm plan}(\tau)
=1.1\max_{s\leq\tau,\ s\in\mathcal T}C_{\mathcal A}^{3\times3}(s),
\qquad \epsilon_T=C_{\mathcal A}^{\rm plan}(\tau)/r^2.
\]

The running maximum avoids exploiting isolated small-patch cancellations.
It is a planning choice for the final-time observable; it does not change the
target to accuracy over an entire curve. The 10% margin is an assumption, not
a bound on finite-size effects. Pointwise costs are retained in the
[ablation table](../results/discussion_v3/ablation.csv).

The 4 × 3 checks at tau = 2 and 4 quantify patch sensitivity. Already at
tau = 2, the plaquette and termwise coefficients differ from 3 × 3 by about
23% and 84%, respectively. Larger-lattice transfer remains uncalibrated.
The main model uses 3 × 3 consistently for L ≥ 4 rather than switching patch
size at a few times. L = 2 uses its own measured 2 × 2 coefficients. STAR uses the termwise coefficient; FT and hybrid use plaquettes.
The experimental NISQ gate proxy transfers the plaquette coefficient, without
claiming its precise ordering matches the experiment.

Only measured time nodes are used. In particular, the tau = L/4 family stops
at L = 32, where tau reaches 8; fixed-time curves still sweep through L = 48.
There is no new fitted universal time exponent. The exposed config supports
pointwise coefficients and a different measured patch as sensitivities.

## Joint error allocation

For each candidate, statistics receive what remains of the same total error:

\[
\epsilon_{\rm stat}=\epsilon-\epsilon_T-\epsilon_L-\epsilon_{\rm syn}
-\Gamma N_T q_T-\epsilon_{\rm res}>0.
\]

The optimizer searches Trotter allowances, logical reserves, synthesis
precision, and neighboring integer step counts. It enforces r ≥ ceil(8 tau)
as a small-step guard. Unused categories cost zero: NISQ has no encoded or
synthesis error, termwise STAR has no supplied T states or synthesis, and
the hybrid does not synthesize its analog rotations.

Mitigated arms reserve 10% of the target for residual estimator bias, including
mismatch and unmodeled SPAM. **Affordability of this reserve is assumed.**
Logical and magic-state faults retain sqrt(F) amplification. FT's modeled
faults are charged to logical/synthesis/magic terms; it has no mitigation
reserve. Each selected allocation is recorded in
[budgets.csv](../results/discussion_v3/budgets.csv).

FT and hybrid both compare the same raw cultivated-state source with zero or
one 15-to-1 cleanup stage and 1–16 sources. Magic error is charged at the model
value, rather than being restricted to a fixed 7.5% share. Other source quality
versus throughput tradeoffs are not invented: optimization is over this
explicit menu. The reported solution is the best candidate in the declared
grid, not a global optimum over architectures. A finer-grid sensitivity improves the tau = 1 examples by at most 18%.
Direct physics checks at selected optimized counts fit the planning envelope
on 3 × 3; neither check validates larger-lattice transfer.

Effective-mitigation and encoded curves use Gaussian planning counts,

\[
M_0=\left\lceil\frac{1.96^2D_*(1-D_*)}{\epsilon_{\rm stat}^2}\right\rceil,
\quad Q_{\rm week}=q\max\left(1,\left\lceil M_0Ft_{\rm shot}/B\right\rceil\right).
\]

PEC references retain the sufficient Hoeffding count. Thus their difference
from the effective line includes statistical convention as well as estimator
choice. The larger difference comes from mitigation exposure.

## Resource laws and source assumptions

| Arm | Main cost law | Declared implementation |
|---|---|---|
| NISQ effective | F = exp(kappa p G2) | Willow gate scale; assumed total mitigation penalty |
| NISQ PEC | log F ≈ 4pG2, or saved full-noise exposure | Specified inverse channels and Hoeffding sampling |
| STAR | log F = 4 alpha p Theta + 4 f R | Termwise Strang; analog rotations |
| FTQC | F = 1; T supply and consumption | Plaquette Strang; synthesized variable angles |
| Hybrid | Analog exponent plus fixed-T supply | Plaquette Strang; analog variable angles |

NISQ uses G2 = 4372 (r/3) (L/6)³ and depth = 294 (r/3) (L/6), anchored to
[Willow Table I](https://arxiv.org/abs/2510.26845v1). The chosen effective
sampling multiplier at the anchor is 100, with 10–1000 sensitivity. It is
not fitted to experimental precision. The PEC alternatives use the exact
2q inverse-channel factor, or saved compiled exposure with bulk terms linear
in r and SPAM kept separate. The compiled proxy recovers the archived r = 4
reference when evaluated with its original statistical allowance.

STAR uses Theta = (5L² − 4L) tau and R = [16L(L−1) + 2L² − 2]r + 1, alpha = 1.5,
f = 1e-5, and 18 code clocks per rotation, from the
[STAR](https://doi.org/10.1103/PRXQuantum.5.010337) and
[small-angle](https://doi.org/10.1103/PhysRevX.15.021057) reference laws with
chosen delivery assumptions. The central lane count is four; scenarios use
8–1 lanes and f = 0–1e-4. Angle exposure is independent of step count; floor
exposure is not. At optimized strict 4 × 4, tau = 8, the floor and angle
terms are 2.28 and 3.07. The earlier floor-dominated example depended on the
cubic-time step prescription.

FT and hybrid share a 1,024-qubit [cultivation](https://arxiv.org/abs/2409.17595)
source with error 2e-9, 100 syndrome cycles per state plus 2d handoff cycles,
and optional cleanup. Supply overlaps consumption. Nine code clocks are
charged per T, and the hybrid adds analog delivery. The same surface-code
fit and per-attempt logical criterion apply to all encoded arms. Data/operation
area is ceil[1.5(2L²+6)] + 6b protected tiles, with source and buffer area added.
Sustained lanes are an assumption rather than a dependency-resolved schedule.

## What changed

At strict precision and tau = 1, 4 × 4 FT selects r = 31 and roughly 36% of the
error for magic states. It avoids cleanup and drops from 5.66 million qubits
to about 1.43 million. Strict 6 × 6 FT still chooses cleanup. The hybrid is
slightly more expensive than termwise STAR at 4 × 4 and less expensive at
6 × 6; the termwise algorithm already avoids fixed T gates.

Optimized effective NISQ uses r = 20 at 4 × 4, assigning about 61% to Trotter
error. Its 36-qubit result is a conditional effective-estimator scenario;
the explicit PEC alternatives remain much more costly. This is not a new
experimental performance claim. Every arm gets the same opportunity to
trade systematic allowance against sampling.

The [ablation](../results/discussion_v3/ablation.csv) first optimizes the
budget under the previous time law, then replaces the law with the measured
envelope. It also compares pointwise coefficients. The capacity bands show
hardware/delivery/mitigation scenarios, not uncertainty in bulk accuracy.

## Reproduction

Edit [discussion_model.json](../configs/discussion_model.json) for delivery
assumptions and [budget_optimization.json](../configs/budget_optimization.json)
for the coefficient policy and optimization grid, then run:

```bash
make PYTHON=/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python discussion
```

This reads saved calibration and regenerates `results/discussion_v3/` costs,
budgets, capacity and ablation tables, and the three current figures. It does
not rerun dynamics. See [REPRODUCING.md](REPRODUCING.md) for calibration and
independent checks. The complete narrative is in [main.pdf](../tex/main.pdf)
with a dedicated bibliography.
