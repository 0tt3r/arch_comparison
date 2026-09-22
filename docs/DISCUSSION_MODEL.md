# Back-of-the-envelope+ for the team discussion

Start with the [figure](../figures/discussion_capacity.pdf) and the
[generated numerical summary](../results/discussion_v2/summary.md).
The [0.01 absolute-error figure](../figures/discussion_capacity_loose.pdf)
shows the precision sensitivity. The four panels cover tau = 1, 4, 8 and
**tau = L/4**. The [original fixed-step example](../results/discussion_v1/summary.md)
is archived, including its figures.

We retain p = 0.001, 10 ns gates, 100 ns measurement/reset, 240 ns syndrome
cycles, U/J = 4, and a week of runtime. The main precision is 0.001 absolute
(1% relative at D* = 0.1). This is an analytical discussion model with editable
assumptions; the refinement uses saved dynamics and compilation data.

## Step-count scenario

The saved 2 × 2 and 3 × 3 plaquette calculations at tau = 1 and r = 16, 32, 64
give error × r² between 0.2178 and 0.2231. We round to 0.22 and choose

\[
\delta D_{\rm model}=0.22\tau^3/r^2,\qquad
r=\left\lceil\sqrt{0.22\tau^3/(\epsilon/4)}\right\rceil.
\]

This gives **r = 10 at epsilon = 0.01 and r = 30 at epsilon = 0.001** for
tau = 1. The full-lattice example allocates the entire 25% systematic reserve
to Trotter error. Reserving half for spatial truncation would instead give
r = 14 and 42: using 10 and 42 mixes budget conventions.

The [fit record](../results/discussion_v2/trotter_fit.json) retains all six
input points. Transfer to larger lattices and other product-formula orderings
is an assumption. The tau³ time dependence is also assumed; it is not measured
by these tau = 1 data. The `trotter_time_power` and `trotter_budget_fraction`
knobs expose both choices. This is a small-patch-informed scenario, not an
accuracy certificate for the larger lattices.

## Resource laws

The effective-mitigation and encoded curves use a Gaussian planning count:

\[
M_0=\left\lceil\frac{1.96^2 D_*(1-D_*)}{(\epsilon/2)^2}\right\rceil,
\quad
Q_{\rm week}=q\max\left(1,\left\lceil M_0 F t_{\rm shot}/B\right\rceil\right).
\]

This gives about 1.38 million ideal samples at the strict target. The explicit
PEC references use the sufficient Hoeffding count, about 7.38 million, to
retain the earlier bound and reproduce the compiled reference. Thus their
difference from the effective curve includes the statistical convention as
well as the mitigation model; the exponential overhead is the larger lever.

| Arm | Main cost law | What sets the constants |
|---|---|---|
| NISQ effective | F = exp(kappa p G2) | Willow circuit counts; chosen total mitigation penalty |
| NISQ PEC | log F ≈ 4pG2, or saved full-noise exposure | Specified channel inversion; Hoeffding sampling |
| STAR | F = exp(4 alpha p Theta + 4 f R) | Angle law, injection floor, mean delivery cost |
| FTQC | F = 1; time set by T supply and consumption | Synthesized T count and rounded cultivation source |
| Hybrid | Analog exponent for variable angles; FT supply for fixed T states | Same plaquette workload and source as FTQC |

For NISQ, use G2 = 4372 (r/3) (L/6)³ and depth = 294 (r/3) (L/6).
These interpolate the 6 × 6, three-step counts from
[Willow Table I](https://arxiv.org/abs/2510.26845v1). Choose F = 100 at that
gate count and p = 0.001; vary it from 10 to 1000. This effective *total*
sampling overhead stands in for observable recovery, selection, and training.
It is **not fitted to experimental accuracy or shot data**, and assumes
remaining mitigation bias is affordable. We use our reference clock.

The NISQ band now also contains two explicit PEC references:

- **Two-qubit noise:** exact depolarizing inverse-channel overhead on the
  rounded gate proxy, log F ≈ 4pG2. This excludes other noise channels.
- **Compiled-noise proxy:** saved full-noise exposure, area, and timing for
  L = 2, 4, 6, 8. Bulk exposure scales linearly with step count; SPAM is kept
  separately. Beyond L = 8, the same cubic routing proxy is used. At L = 4,
  r = 4, this recovers the earlier 2,132,604-qubit estimate. Other step counts
  and larger sizes are extrapolations, not new compiled schedules.

These references do not bound the performance of every mitigation strategy.
At strict precision the refined central effective curve reaches 4 × 4 with
about 1,400 qubits; the PEC references remain at 2 × 2 over the displayed
budgets. Their separation is a useful question for the group, not a claim
that the effective curve achieves the target bias.

For STAR, use termwise second-order evolution. At U/J = 4,
Theta = (5L² − 4L) tau and R = [16L(L−1) + 2L² − 2] r + 1.
Set alpha = 1.5, f = 1e-5, and 18 code clocks per completed rotation, divided
by four effective parallel lanes. These combine the
[STAR](https://doi.org/10.1103/PRXQuantum.5.010337) and
[small-angle](https://doi.org/10.1103/PhysRevX.15.021057) reference laws with
chosen delivery assumptions. The scenario range uses f = 0–1e-4 and 8–1 lanes.

FTQC uses the existing plaquette counts and logarithmic synthesis fit. The
hybrid uses those same plaquette counts, but injects the variable rotations
and supplies only the fixed T states. Both use a 1,024-qubit source with state
error 2e-9, roughly 100 syndrome cycles per state plus 2d handoff cycles. This
rounds the existing [cultivation](https://arxiv.org/abs/2409.17595) reference.
Add one 15-to-1 cleanup stage when the total magic error exceeds its budget,
including the hybrid's sqrt(F) amplification. Charge nine code clocks per
consumed T state, optimize 1–16 sources, and overlap supply with the combined
analog/T consumption time. Vary source service time by 0.5–4 and lanes by 8–1.

Encoded arms retain the surface-code fit and choose distance from the
per-attempt fault budget, including sqrt(F) amplification where F > 1.
Allocate ceil[1.5(2L²+6)] + 6b protected tiles for data and b operation lanes,
plus source/buffer tiles when needed. Lanes represent sustained delivery,
without an explicit dependency schedule. The arms use different second-order
implementations suited to their gate costs; the shared step prescription is
an additional approximation. All evolve the full lattice.

## What changes the comparison

At 4 × 4, tau = 1, changing r = 4 to 30 raises the effective NISQ exponent
from 1.82 to 13.64. STAR's exponent goes from 0.420 to 0.650, entirely through
the floor term. FT's T count grows from 20,736 to 166,744; crossing the chosen
cleanup threshold contributes to its area jump from about 79,000 to 5.66
million qubits. [Step sensitivity data](../results/discussion_v2/step_sensitivity.csv)
separate these effects. The hybrid is slightly more expensive than termwise
STAR at 4 × 4 and about 2.5 times cheaper at 6 × 6 under the strict central
assumptions. Termwise STAR already avoids the fixed T gates that motivate
the hybrid; a hybrid advantage is not automatic.

At fixed tau, refinement enters the NISQ exponent through rL³, while STAR's
angle term is independent of r and its floor scales as frL². With the chosen
tau³/r² error law, tau ∝ L gives r ∝ L^(3/2): NISQ's exponent scales as
L^(9/2), STAR's angle term as L³, and its floor as L^(7/2). An L⁴ NISQ law
would instead assume r ∝ tau; choose time power two to explore that scenario.
These are laws of the chosen implementations. Routing, observable support,
and lane scaling can change them. Bands show scenario ranges, not confidence
intervals or bounds on unmodeled effects.

To change an assumption, edit
[discussion_model.json](../configs/discussion_model.json), then run:

```bash
make PYTHON=/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python discussion
```

The formulas are in [discussion_model.py](../python/discussion_model.py).
Outputs include the assumption snapshot, small-patch fit inputs, 2,688 scenario
points, CSV capacity/step-sensitivity tables, and a short Markdown summary in
`results/discussion_v2/`. `make paper` also rebuilds the figures and incorporates
this discussion model into the main manuscript.
