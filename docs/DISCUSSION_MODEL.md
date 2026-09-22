# Back-of-the-envelope+ for the team discussion

Start with the [figure](../figures/discussion_capacity.pdf) and the
[generated numerical summary](../results/discussion_v1/summary.md).
The [0.01 absolute-error figure](../figures/discussion_capacity_loose.pdf)
shows the precision sensitivity. This model is intended for a discussion of
mechanisms and plausible scales. It does not require new compiler runs,
exact dynamics, or experimental fitting.

We retain p = 0.001, 10 ns gates, 100 ns measurement/reset, 240 ns syndrome
cycles, U/J = 4, and a week of runtime. The main precision is 0.001 absolute
(1% relative at D* = 0.1). Set r = ceil(4 tau). We assume Trotter and residual
mitigation bias fit their budgets for this illustration; r is not calibrated
to that precision. All arms use a Gaussian planning count

\[
M_0=\left\lceil\frac{1.96^2 D_*(1-D_*)}{(\epsilon/2)^2}\right\rceil,
\quad
Q_{\rm week}=q\max\left(1,\left\lceil M_0 F t_{\rm shot}/B\right\rceil\right).
\]

This gives about 1.38 million ideal samples at the strict target. The
Gaussian count is a planning approximation; the older Hoeffding bound and
compiled-reference curves remain available separately.

| Arm | Main cost law | What sets the constants |
|---|---|---|
| NISQ | F = exp(kappa p G2) | Willow circuit counts; chosen effective mitigation penalty |
| STAR | F = exp(4 alpha p Theta + 4 f R) | Angle law, injection floor, and mean delivery cost |
| FTQC | F = 1; time set by T supply and consumption | Synthesized T count and rounded cultivation source |

For NISQ, use G2 = 4372 (r/3) (L/6)^3 and depth = 294 (r/3) (L/6).
These interpolate the 6 × 6, three-step counts from
[Willow Table I](https://arxiv.org/abs/2510.26845v1). Choose F = 100 at that
gate count and p = 0.001; vary it from 10 to 1000. This effective *total*
sampling overhead stands in for observable recovery, selection, and training.
It is **not fitted to experimental accuracy or shot data**. We use our fast
reference clock, not the experiment's measured wall-clock runtime.

For STAR, use termwise second-order evolution. At U/J = 4,
Theta = (5L² − 4L) tau and R = [16L(L−1) + 2L² − 2] r + 1.
Set alpha = 1.5, f = 1e-5, and 18 code clocks per completed rotation, divided
by four effective parallel lanes. These are the existing
[STAR](https://doi.org/10.1103/PRXQuantum.5.010337) and
[small-angle](https://doi.org/10.1103/PhysRevX.15.021057) reference laws and
chosen delivery assumptions. The scenario range uses f = 0–1e-4 and 8–1 lanes.

For FTQC, use the existing plaquette counts and logarithmic synthesis fit.
Supply states at error 2e-9 from a 1,024-qubit source, roughly 100 syndrome
cycles per state plus 2d handoff cycles. This rounds the existing
[cultivation](https://arxiv.org/abs/2409.17595) reference, rather than adding a
new factory calculation. Add one 15-to-1 cleanup stage when the total magic
error exceeds its per-shot budget. Charge nine code clocks per consumed T
state divided by four effective lanes, and optimize 1–16 supply sources.
Vary source service time by 0.5–4 and lanes by 8–1.

Both encoded arms retain the existing surface-code fit and choose distance
from the per-attempt fault budget. Allocate ceil[1.5(2L²+6)] + 6b protected
tiles for data and b operation lanes, plus FT source/buffer tiles. Parallel
lanes represent sustained delivery, without an explicit dependency schedule.
The arms use different second-order implementations suited to their gate
costs; this is not a matched decomposition comparison. All costs evolve the
full lattice, without physical light-cone truncation.

The useful scaling contrast is immediate: the chosen NISQ routing model
has log F proportional to rL³, STAR has log F proportional to pL²tau + frL²,
and FTQC has polynomial cost with logarithmic synthesis/protection factors.
Changing the routing, observable support, or lane scaling changes these
particular powers. The band shows alternative constants, not a confidence
interval or a bound on unmodeled effects.

To change an assumption, edit
[discussion_model.json](../configs/discussion_model.json), then run:

```bash
make PYTHON=/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python discussion
```

The formulas are in [discussion_model.py](../python/discussion_model.py).
Outputs include the full assumption snapshot, 1,296 scenario points, CSV
capacity tables, and the short Markdown summary. `make paper` also rebuilds
the figures and incorporates this discussion model into the main manuscript.
