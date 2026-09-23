# Current comparison: component scenarios and NISQ requirements

The current results are in [integrated_v1](../results/integrated_v1/summary.md).
They supersede the old scalar-coefficient capacity plot as the main discussion
figure. Old calculations and their inputs remain available for replay.
This pass integrates existing evidence; it does not reconstruct an experimental
mitigation pipeline or launch larger noisy-dynamics simulations.

## Shared task and what a band means

The workload remains the half-filled Neel quench at U/J=4 with an open square
lattice and a central doublon observable. We evaluate tau=1,2, absolute error
0.01 and 0.001, and p=0.001 and 0.0001. The clock and the one-week cutoff are
unchanged. Each encoded arm reselects steps, error allocations and source counts.
The Trotter coefficient is the existing running finite-patch envelope, including
4x3 and 4x4; its transfer to larger sizes remains an assumption. All angles use
exp(-i theta P).

The three scenarios are finite, declared alternatives rather than confidence
intervals. The central component scenario uses the measured growth-campaign
upper limits in the injection submodel. The zero-growth scenario sets the
distance-9 growth component to zero, an optimistic sensitivity consistent with
its zero-failure count but not certified by that count. The slow scenario
increases a declared detector-rejection coefficient from one to four and slows
the cultivated-state source by four. Neither change is a measured full-device
error bar. Source error remains 2e-9 at both physical error rates.

A robust encoded capacity leader has larger capacity even in its least
favorable scenario than both other encoded arms in their most favorable
scenarios. This criterion is relative to this menu. Sparse lattice sizes,
integer replicas and the size ceiling are retained in the tables. No
unconditional NISQ boundary is inferred.

## STAR and hybrid: keep error and acceptance together

For each actual workload angle we evaluate the existing generator-Pauli
instrument in `star_channel_study.py`. It includes perfect syndrome projection,
independent generator faults, noisy probabilistic coherent error cancellation,
and repeat-until-success (RUS) branches. Preparation rejection is modeled with
the existing declared detector factor. Generic syndrome/preparation faults are
not newly simulated. Thus this is a specified component composition rather than
end-to-end noisy injection validation.

The available injection seeds are 5, 9 and 25. Growth ends at distance 13 or 25;
only seed/target pairs present in the component data are used.
Seeds 5 and 9 use their own simultaneous X/Z campaign upper limits from
`results/frontier_v1/growth_summary.json`; no distance law is fitted through
unmeasured seeds. Seed 25 does not grow. The zero-growth scenario excludes seed
5, whose growth experiment measured failures at p=1e-3. All seeds pay their
angle-dependent acceptance and retry time, so selecting a larger seed is not a
free reduction of a floor.

The growth measurements start with an ideally encoded state and include growth,
memory and readout. In this scenario they are interpreted as an independent
dephasing contribution at each growth inside the RUS and correction processes.
Their sum as a Pauli-channel bound requires the common-Pauli-channel assumption;
placing that bound in the injection instrument is an additional composition
assumption. Some protected faults are consequently charged conservatively in
both the growth component and the logical volume. These upper limits are not
measurements of arbitrary-angle injected states.

The old helper parameterizes growth by a distance-9 anchor. We invert that
parameter conversion so that its growth calls receive exactly the specified
seed's saved upper limit. The old distance/power-law growth extrapolation is
therefore not used to predict an unmeasured seed.

For each rotation channel, let c be its complex transverse multiplier and
q=(1-Re(c))/2. The campaign uses

\[
\chi=-2\sum_j\log(1-2q_j),\qquad
b_{\rm coh}\leq\tfrac12\left[\prod_j
 (1+|\operatorname{Im}c_j|/(1-2q_j))-1\right].
\]

Multiplicity is included in the sums/products. The coherent bound and the RUS
tail bound consume explicit additional error allowance. We retain the previous
10% residual-error reserve, separately, for unmodeled effects. The old alpha=1.5
law is shown only as a historical reference; alpha=2 to 5 is not silently treated
as a universal full-gadget fit.

The receiver distance is at least the selected growth/channel anchor d_f.
For a larger selected data distance d, the reference service time is stretched
by d/d_f and an additional 2d-cycle protected handoff is charged per rotation.
Those extra operations enter the protected tile-time volume. This is a declared
continuation of the layout model; no new growth-channel calibration at that
larger distance is claimed. The distance-13 alternative avoids imposing the
distance-25 footprint at low p, and uses its own channel and growth data.
Hybrid additionally respects the cultivated source's distance-15 handoff floor.

The analog service time includes state-preparation rejection, teleportation and
noisy correction retries. Four analog lanes and six support tiles per lane are
the retained aggregate parallelism assumption. Clifford scheduling outside
these macros remains approximate. Hybrid supplies the plaquette circuit's fixed
T rotations with the same source menu as FT, while its dynamic rotations use
this analog instrument. Supply overlaps consumption; serial analog and digital
consumption are added before taking the maximum with source service.

## FT synthesis and delivery

FT selects individual synthesis or equal-angle Hamming-weight phasing in batches
4, 8 or 16. Arithmetic counts and the measurement-assisted uncomputation have
already been checked. Compact and fast data blocks each have one T-consumption
port; a fast block does not acquire four T ports by dividing its consumption
rate by the analog lane count. Workspace, source copies, cleanup, waiting and
logical faults are charged. The source menu and quality are shared with hybrid.

These are aggregate delivery scenarios. The separately measured 5-16x slowdown
of one routed onsite construction is not multiplied into the whole program as
though it were a full-program measurement. Full routed schedules could change
the ranking of Hamming-weight batches and the architecture ordering.

## NISQ: solve for the required accuracy and sampling cost

The physical circuit scale uses the existing experimentally anchored native
proxy. This does not fit experimental estimator performance. For a selected
step count r, physical footprint a and attempt duration t, the number of
evaluation attempts available at Q qubits and deadline T is

\[
 M=\max\{0,\lfloor Q/a\rfloor\lfloor T/t\rfloor-M_{\rm train}\}.
\]

Training costs are either zero or six million attempts of the same duration.
These are explicit sensitivities, not claims about an experimental training
protocol. F is the estimator variance per attempted circuit relative to the
reference Bernoulli variance v*=0.1(1-0.1). Postselection rejection is already
included in F and must not be charged again. Training-fit uncertainty must be
included by an actual estimator in its F or remaining-offset allowance.

With Trotter allowance e_T=C/r^2, remaining systematic offset b, and Gaussian
planning factor z=1.96, the accuracy condition is

\[
 e_T+b+z\sqrt{Fv_*/M}\leq\epsilon.
\]

Consequently,

\[
 F_{\max}(b,r)=\frac{M(\epsilon-e_T-b)^2}{z^2v_*},\qquad
 b_{\max}(F,r)=\epsilon-e_T-z\sqrt{Fv_*/M}.
\]

An exhausted error budget or insufficient attempts is infeasible, rather than a
zero-offset certificate. The step menu is scanned directly and each selected
r is saved. An estimator must achieve its F and b at that r; the maximization
does not assume that its overhead stays fixed when steps change. The optional
kappa column translates that point to F=exp(kappa p G2), without fitting that
law to trained mitigation.

The deadline is the smaller of one week and the selected encoded runtime divided
by the requested speedup (one or ten). At each Q the comparator is the fastest
fitting member among the three selected encoded implementations for that
scenario. Those implementations minimize one-week area; this does not claim
fresh fixed-Q runtime optimization over every factory/allocation choice.
If no encoded copy fits, the table states `none_fits` and evaluates only the
one-week NISQ requirement. It does not claim a measured NISQ win.

These are Gaussian planning requirements. They are not a confidence certificate
for an unspecified trained estimator, and make no claim about whether the
displayed offset is achievable. Their purpose is to identify when a targeted
accuracy test could change an architecture conclusion.

## Reproduction

Run `make integrated` with the qb_gsee interpreter. The scalar channel cache
and completed resource rows are restartable; input changes fail closed rather
than mixing contracts. The run reads historical results without changing them.
`make integrated-figures` redraws the saved result; `make integrated-check`
checks input hashes, error budgets, replicas and the threshold identities.
The main manuscript and presentation manifest include this result version.
